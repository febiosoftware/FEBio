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
#include "FEElasticMultiscaleDomain2O.h"
#include "FEMicroMaterial2O.h"
#include "FERVEProbe.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"
#include <FECore/log.h>


//-----------------------------------------------------------------------------
//! constructor
FEElasticMultiscaleDomain2O::FEElasticMultiscaleDomain2O(FEModel* pfem) : FEElasticSolidDomain2O(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialize element data
bool FEElasticMultiscaleDomain2O::Init()
{
	// initialize base class first
	if (FEElasticSolidDomain2O::Init() == false) return false;

	FEModel& fem = *GetFEModel();

	// initialze RVEs
	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
	FEMesh& m = *GetMesh();
		
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);
	FERVEModel2O& rve = pmat->m_mrve;
			
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j)
		{
			x0[j] = m.Node(el.m_node[j]).m_r0;
			xt[j] = m.Node(el.m_node[j]).m_rt;
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
			mmpt2O.m_elem_id = el.GetID();
			mmpt2O.m_gpt_id = j;

			// initialize the material point RVE
			// This essentially copies the parent RVE to the material point RVE
			if (mmpt2O.m_rve.Init(rve) == false) return false;
		}
	}

	// initialize surface RVEs
	int nnf = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& face = m_surf.Element(i);
		int nint = face.GaussPoints();
		for (int n=0; n<nint; ++n, ++nnf)
		{
			for (int k=0; k<2; ++k)
			{
				FEMaterialPoint& mp = *m_surf.GetData(nnf).m_pt[k];
				FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
				mmpt2O.m_elem_id = -i-1;
				mmpt2O.m_gpt_id = n + k*nint;

				// Initialize the material point RVE
				// This essentially copies the parent RVE model to the material points
				mmpt2O.m_rve.Init(rve);
			}
		}
	}

		// TODO: I need to move the code below to FERVEProbe::Init.
/*
		if (p.m_neid > 0)
		{
			FEElement* pel = FindElementFromID(p.m_neid);
			if (pel == 0)
			{
				feLogError("Invalid Element ID for micro probe %d in material %d (%s)", i + 1, m_pMat->GetID(), m_pMat->GetName().c_str());
				return false;
			}

			int nint = pel->GaussPoints();
			int ngp = p.m_ngp - 1;
			if ((ngp>=0)&&(ngp<nint))
			{
				FEMaterialPoint& mp = *pel->GetMaterialPoint(ngp);
				FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
				FERVEProbe* prve = new FERVEProbe(fem, mmpt.m_rve, p.m_szfile.c_str());
				p.m_probe = prve;
				prve->SetDebugFlag(p.m_bdebug);
			}
			else
			{
				feLogError("Invalid gausspt number for micro-probe %d in material %d (%s)", i + 1, m_pMat->GetID(), m_pMat->GetName().c_str());
				return false;
			}
		}
		else 
		{
			int fid = -p.m_neid-1;
			if ((fid < 0) || (fid >= m_surf.Elements()))
			{
				feLogError("Invalid surface ID for micro-probe");
				return false;
			}
			int nnf = 0;
			for (int j=0; j<fid; ++j)
			{
				FESurfaceElement& el = m_surf.Element(j);
				int nint = el.GaussPoints();
				nnf += nint;
			}

			FESurfaceElement& face = m_surf.Element(fid);
			int nint = face.GaussPoints();
			int k = (p.m_ngp-1)/nint;
			int gpt = (p.m_ngp-1)%nint;
			nnf += gpt;

			FEMaterialPoint& mp = *m_surf.GetData(nnf).m_pt[k];
			FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
			if (mmpt2O.m_elem_id != p.m_neid ) return false;
			if (mmpt2O.m_gpt_id  != p.m_ngp-1) return false;

			FERVEProbe* prve = new FERVEProbe(fem, mmpt2O.m_rve, p.m_szfile.c_str());
			p.m_probe = prve;
			prve->SetDebugFlag(p.m_bdebug);
		}
*/

	return true;
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::Update(const FETimeInfo& timeInfo)
{
	try
	{
		// call base class
		FEElasticSolidDomain2O::Update(timeInfo);
	}
	catch (FEMultiScaleException)
	{
		// store all the probes
		FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);
		int NP = pmat->Probes();
		for (int i=0; i<NP; ++i)
		{
			FERVEProbe& p = pmat->Probe(i);
			if (p.GetDebugFlag()) p.Save();
		}

		// retrhow
		throw;
	}
}
