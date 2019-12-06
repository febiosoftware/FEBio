/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

BEGIN_FECORE_CLASS(FESBMPointSource, FEBodyLoad)
	ADD_PARAMETER(m_sbm, "sbm");
	ADD_PARAMETER(m_pos.x, "x");
	ADD_PARAMETER(m_pos.y, "y");
	ADD_PARAMETER(m_pos.z, "z");
	ADD_PARAMETER(m_val, "value");
END_FECORE_CLASS();

FESBMPointSource::FESBMPointSource(FEModel* fem) : FEBodyLoad(fem)
{
	m_sbm = -1;
	m_pos = vec3d(0,0,0);
	m_val = 0.0;
	m_valp = 0.0;

	m_closestPoint = nullptr;
	m_local_sbm = -1;
}

bool FESBMPointSource::Init()
{
	if (m_sbm == -1) return false;
	return FEBodyLoad::Init();
}

void FESBMPointSource::Update()
{
	int localId = 0;
	FEMaterialPoint* mp = FindClosestMaterialPoint(localId);

	if (mp != m_closestPoint)
	{
		if (m_closestPoint == nullptr)
		{
			FESolutesMaterialPoint& pd = *(mp->ExtractData<FESolutesMaterialPoint>());
			pd.m_sbmr[localId] = m_val;
			pd.m_sbmrp[localId] = m_val;
		}
		else
		{
			FESolutesMaterialPoint& ps = *(m_closestPoint->ExtractData<FESolutesMaterialPoint>());
			FESolutesMaterialPoint& pd = *(mp->ExtractData<FESolutesMaterialPoint>());

			pd.m_sbmr[localId] = ps.m_sbmr[m_local_sbm];
			pd.m_sbmrp[localId] = ps.m_sbmrp[m_local_sbm];

			ps.m_sbmr[m_local_sbm] = 0.0;
			ps.m_sbmrp[m_local_sbm] = 0.0;
		}

		m_closestPoint = mp;
		m_local_sbm = localId;
	}
}

/*void FESBMPointSource::Update()
{
	int localId = 0;
	FEMaterialPoint* mp = FindClosestMaterialPoint(localId);

	// reset current closest point
	if ((mp != m_closestPoint) && (m_closestPoint!=nullptr))
	{
		FESolutesMaterialPoint& ps = *(m_closestPoint->ExtractData<FESolutesMaterialPoint>());
		ps.m_sbmr[m_local_sbm] -= m_valp;
		ps.m_sbmrp[m_local_sbm] -= m_valp;
	}

	// set the sbm concentration to the user-specified value
	if (mp != m_closestPoint)
	{
		m_closestPoint = mp;
		m_local_sbm = localId;

		FESolutesMaterialPoint& ps = *(mp->ExtractData<FESolutesMaterialPoint>());
		m_valp = m_val;
		ps.m_sbmr[localId] += m_val;
		ps.m_sbmrp[localId] += m_val;
	}
}
*/

void FESBMPointSource::UpdatePos(vec3d pos)
{
	if (this) { m_pos = pos; }
}

void FESBMPointSource::UpdateSBM(int id, double val)
{
	if (this) { m_sbm = id;	m_val = val; }
}

FEMaterialPoint* FESBMPointSource::FindClosestMaterialPoint(int& localID)
{
	double minDist = 1e99;
	FEMaterialPoint* minPt = nullptr;
	localID = -1;

	// find the integration point that is closest to m_pos
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		// get the domain and make sure it has a multiphasic material
		FEDomain* dom = &mesh.Domain(i);
		FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
		if (mat)
		{
			// does this material have the correct sbm
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

			if (sbmid >= 0)
			{
				int NE = dom->Elements();
				for (int j = 0; j<NE; ++j)
				{
					FEElement& el = dom->ElementRef(j);
					int nint = el.GaussPoints();
					int neln = el.Nodes();

					vec3d r0[FEElement::MAX_NODES];
					for (int k=0; k<neln; ++k) r0[k] = mesh.Node(el.m_node[k]).m_r0;

					for (int k = 0; k<nint; ++k)
					{
						FEMaterialPoint* mp = el.GetMaterialPoint(k);
						
						vec3d rk = el.Evaluate(r0, k);

						double dk2 = (rk - m_pos)*(rk - m_pos);
						if (dk2 < minDist)
						{
							minDist = dk2;
							minPt = mp;
							localID = sbmid;
						}
					}
				}
			}
		}
	}

	return minPt;
}
