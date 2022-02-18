/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "FERVEProbe.h"
#include "FEMicroMaterial.h"
#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/FEModel.h>
#include <FECore/FEDomain.h>
#include <FECore/log.h>

FERVEProbe::FERVEProbe(FEModel* fem) : FECallBack(fem, CB_ALWAYS)
{
	m_rve = nullptr;
	m_file = "rve.xplt";
	m_bdebug = false;
}

bool FERVEProbe::Init()
{
	// make sure we have an RVE
	if (m_rve == nullptr) return false;
	return FECallBack::Init();
}

bool FERVEProbe::Execute(FEModel& fem, int nwhen)
{
	if (nwhen == CB_INIT)	// initialize the plot file
	{
		assert(m_rve);
		if (m_rve == nullptr) return false;

		// create a plot file
		m_xplt = new FEBioPlotFile(m_rve);
		if (m_xplt->Open(m_file.c_str()) == false)
		{
			feLog("Failed creating probe.\n\n");
			delete m_xplt; m_xplt = 0;
		}

		// write the initial state
		Save();
	}
	else if (nwhen == CB_MINOR_ITERS)
	{
		if (m_bdebug) Save();
	}
	else if (nwhen == CB_MAJOR_ITERS)	// store the current state
	{
		Save();
	}
	else if (nwhen == CB_SOLVED)	// clean up
	{
		if (m_xplt) delete m_xplt;
		m_xplt = 0;
	}

	return true;
}

void FERVEProbe::Save()
{
	assert(m_rve);
	if (m_rve)
	{
		if (m_xplt) m_xplt->Write((float)m_rve->GetCurrentTime());
	}
}

void FERVEProbe::SetRVEModel(FEModel* rve)
{
	m_rve = rve;
}

//=============================================================================
BEGIN_FECORE_CLASS(FEMicroProbe, FERVEProbe)
	ADD_PARAMETER(m_neid  , "element_id");
	ADD_PARAMETER(m_ngp   , "gausspt"   );
	ADD_PARAMETER(m_file  , "file"      );
	ADD_PARAMETER(m_bdebug, "debug"     );
END_FECORE_CLASS();

FEMicroProbe::FEMicroProbe(FEModel* fem) : FERVEProbe(fem)
{
	m_neid = -1;	// invalid element - this must be defined by user
	m_ngp = -1;		// invalid gauss point
}

bool FEMicroProbe::Init()
{
	// get the (parent) model
	FEModel& fem = *GetFEModel();

	// get the micro-material
	FEMicroMaterial* mat = dynamic_cast<FEMicroMaterial*>(GetParent());
	if (mat == nullptr)
	{
		feLogError("The RVE probe must be a property of a micro-material.");
		return false;
	}

	// find the element from the ID
	FEMesh& mesh = fem.GetMesh();
	FEElement* pel = mesh.FindElementFromID(m_neid);
	if (pel == nullptr)
	{
		feLogError("Invalid Element ID for micro probe %d in material %d (%s)", m_neid, mat->GetID(), mat->GetName().c_str());
		return false;
	}

	// make sure the element's parent domain has this material
	FEDomain* dom = dynamic_cast<FEDomain*>(pel->GetMeshPartition());
	if ((dom == nullptr) || (dom->GetMaterial() != mat))
	{
		feLogError("The specified element does not belong to this material's domain.");
		return false;
	}

	// get the gauss-point
	int nint = pel->GaussPoints();
	if ((m_ngp >= 0) && (m_ngp < nint))
	{
		FEMaterialPoint* mp = pel->GetMaterialPoint(m_ngp);
		FEMicroMaterialPoint* mmp = mp->ExtractData<FEMicroMaterialPoint>();
		if (mmp == nullptr) return false;
		SetRVEModel(&mmp->m_rve);
	}
	else
	{
		feLogError("Invalid gausspt number for micro-probe %d in material %d (%s)", m_ngp, mat->GetID(), mat->GetName().c_str());
		return false;
	}


	return FERVEProbe::Init();
}
