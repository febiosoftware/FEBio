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
#include "FESBMPointSource.h"
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <algorithm>
#include <iostream>

BEGIN_FECORE_CLASS(FESBMPointSource, FEBodyLoad)
	ADD_PARAMETER(m_sbm, "sbm");
	ADD_PARAMETER(m_pos.x, "x");
	ADD_PARAMETER(m_pos.y, "y");
	ADD_PARAMETER(m_pos.z, "z");
	ADD_PARAMETER(m_val, "value");
	ADD_PARAMETER(m_weighVolume, "weigh_volume");
END_FECORE_CLASS();

FESBMPointSource::FESBMPointSource(FEModel* fem) : FEBodyLoad(fem), m_search(&fem->GetMesh())
{
	static bool bfirst = true;
	m_sbm = -1;
	m_pos = vec3d(0,0,0);
	m_val = 0.0;
	m_reset = bfirst;
	m_doReset = true;
	m_weighVolume = true;
	bfirst = false;
}

bool FESBMPointSource::Init()
{
	if (m_sbm == -1) return false;
	if (m_search.Init() == false) return false;
	return FEBodyLoad::Init();
}

// allow species to accumulate at the point source
void FESBMPointSource::Accumulate(double dc) {
	double rt[3] = { 0, 0, 0 };
	FEElement* el = m_search.FindElement(m_pos, rt);
	if (el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;

	int sbmid = -1;
	int sbms = mat->SBMs();
	for (int j = 0; j < sbms; ++j)
	{
		int sbmj = mat->GetSBM(j)->GetSBMID();
		if (sbmj == m_sbm)
		{
			sbmid = j;
			break;
		}
	}
	if (sbmid == -1) return;

	m_val = dc + m_val; // prevent negative concentrations
	m_accumulate = true;
}

void FESBMPointSource::Update()
{
	if (m_reset && m_doReset) ResetSBM();

	// find the element in which the point lies
	double rt[3] = { 0, 0, 0 };
	FEElement* el = m_search.FindElement(m_pos, rt);
	if (el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;

	// calculate the element volume
	FEMesh* mesh = dom->GetMesh();
	double Ve = mesh->ElementVolume(*el);

	// we prescribe the element average to the integration points
	const int nint = el->GaussPoints();
	double val = (m_weighVolume ? m_val / Ve : m_val);

	// Make sure the material has the correct sbm
	int sbmid = -1;
	int sbms = mat->SBMs();
	for (int j = 0; j<sbms; ++j)
	{
		int sbmj = mat->GetSBM(j)->GetSBMID();
		if (sbmj == m_sbm)
		{
			sbmid = j;
			break;
		}
	}
	if (sbmid == -1) return;

	// set the concentration of all the integration points
	double H[FEElement::MAX_NODES];
	double m_q[3];
	m_q[0] = m_q[1] = m_q[2] = 0.0;
	FESolidElement* m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	if (m_el == nullptr) return;
	m_el->shape_fnc(H, m_q[0], m_q[1], m_q[2]);
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint* mp = el->GetMaterialPoint(i);
		FESolutesMaterialPoint& pd = *(mp->ExtractData<FESolutesMaterialPoint>());
		// if this point source has not yet been added to the integration points do it then turn off the flag so we don't double count
		if (m_accumulate) {

			pd.m_sbmr[sbmid] = std::max(0.0, H[i] * val + pd.m_sbmrp[sbmid]);
			pd.m_sbmrp[sbmid] = pd.m_sbmr[sbmid];
		}
		else {
			pd.m_sbmr[sbmid] = pd.m_sbmr[sbmid];
			pd.m_sbmrp[sbmid] = pd.m_sbmrp[sbmid];
		}
	}
	m_accumulate = false; // don't double count a point source
	m_val = 0;
}

void FESBMPointSource::SetPosition(const vec3d& pos)
{
	m_pos = pos;
}

vec3d FESBMPointSource::GetPosition() const
{
	return m_pos;
}

void FESBMPointSource::SetSBMID(int id)
{
	m_sbm = id;	
}

int FESBMPointSource::GetSBMID() const
{
	return m_sbm;
}

void FESBMPointSource::SetValue(double val)
{
	m_val = val;
}

double FESBMPointSource::GetValue() const
{
	return m_val;
}

void FESBMPointSource::SetWeighVolume(bool b)
{
	m_weighVolume = b;
}

void FESBMPointSource::SetResetFlag(bool b)
{
	m_doReset = b;
}

void FESBMPointSource::SetAccumulateFlag(bool b)
{
	m_accumulate = b;
}

void FESBMPointSource::ResetSBM()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();

		FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom.GetMaterial());
		if (mat)
		{
			// Make sure the material has the correct sbm
			int sbmid = -1;
			int sbms = mat->SBMs();
			for (int j = 0; j<sbms; ++j)
			{
				int sbmj = mat->GetSBM(j)->GetSBMID();
				if (sbmj == m_sbm)
				{
					sbmid = j;
					break;
				}
			}

			if (sbmid != -1)
			{
				for (int j = 0; j < NE; ++j)
				{
					FEElement& el = dom.ElementRef(j);

					// set the concentration of all the integration points
					int nint = el.GaussPoints();
					for (int k = 0; k < nint; ++k)
					{
						FEMaterialPoint* mp = el.GetMaterialPoint(k);
						FESolutesMaterialPoint& pd = *(mp->ExtractData<FESolutesMaterialPoint>());
						pd.m_sbmr[sbmid] = 0.0;
						pd.m_sbmrp[sbmid] = 0.0;
					}
				}
			}
		}
	}
}
