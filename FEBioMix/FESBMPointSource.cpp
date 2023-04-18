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
#include "FESolute.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEElemElemList.h>
#include <algorithm>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <FECore/FEAnalysis.h>

BEGIN_FECORE_CLASS(FESBMPointSource, FEBodyLoad)
	ADD_PARAMETER(m_sbmId, "sbm");
	ADD_PARAMETER(m_rate, "rate");
	ADD_PARAMETER(m_pos.x, "x");
	ADD_PARAMETER(m_pos.y, "y");
	ADD_PARAMETER(m_pos.z, "z");
	ADD_PARAMETER(m_weighVolume, "weigh_volume");
END_FECORE_CLASS();

FESBMPointSource::FESBMPointSource(FEModel* fem) : FEBodyLoad(fem), m_search(&fem->GetMesh())
{
	//static bool bfirst = true;
	m_sbmId = -1;
	m_pos = vec3d(0, 0, 0);
	m_rate = 0.0;
	//m_reset = bfirst;
	//m_doReset = true;
	m_weighVolume = true;
	//bfirst = false;
}

vec3d FESBMPointSource::GetPosition() const
{
	return m_pos;
}

void FESBMPointSource::SetPosition(const vec3d& pos)
{
	m_pos = pos;
}

int FESBMPointSource::GetSBMID() const
{
	return m_sbmId;
}


void FESBMPointSource::SetSBMID(int sbmID)
{
	m_sbmId = sbmID;
}

double FESBMPointSource::GetRate() const
{
	return m_rate;
}

void FESBMPointSource::SetRate(double rate)
{
	m_rate = rate;
}

double FESBMPointSource::GetdC() const
{
	return m_dC;
}

double FESBMPointSource::GetdCp() const
{
	return m_dCp;
}

void FESBMPointSource::SetdC(double dC)
{
	m_dC = dC;
}


void FESBMPointSource::SetRadius(double radius)
{
	m_radius = radius;
	m_Vc = (4.0 / 3.0) * PI * pow(m_radius, 3.0);
}

void FESBMPointSource::SetAccumulateFlag(bool b)
{
	m_accumulate = b;
}

//void FESBMPointSource::SetWeighVolume(bool b)
//{
//	m_weighVolume = b;
//}

//void FESBMPointSource::SetResetFlag(bool b)
//{
//	m_doReset = b;
//}



bool FESBMPointSource::Init()
{
	if (m_sbmId == -1) return false;
	if (m_search.Init() == false) return false;
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	return FEBodyLoad::Init();
}

// allow species to accumulate at the point source
void FESBMPointSource::Accumulate(double dc) {
	// find the element in which the point lies
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	if (m_el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;

	// Make sure the material has the correct solute
	int sbmid = -1;
	int sbms = mat->SBMs();
	for (int j = 0; j < sbms; ++j)
	{
		int sbmj = mat->GetSBM(j)->GetSBMID();
		if (sbmj == m_sbmId)
		{
			sbmid = j;
			break;
		}
	}
	if (sbmid == -1) return;

	m_rate = dc + m_rate;
	m_accumulate = true;
}

void FESBMPointSource::Update()
{
	//if (m_reset && m_doReset) ResetSBM();
	if (m_accumulate) {
		m_dCp = m_dC;
		// find the element in which the point lies
		m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
		if (m_el == nullptr) return;

		// make sure this element is part of a multiphasic domain
		FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
		FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
		if (mat == nullptr) return;

		// calculate the element volume
		FEMesh* mesh = dom->GetMesh();
		double Ve = mesh->ElementVolume(*m_el);

		// we prescribe the element average to the integration points
		const int nint = m_el->GaussPoints();
		// Make sure the material has the correct sbm
		int sbmid = -1;
		int sbms = mat->SBMs();
		for (int j = 0; j < sbms; ++j)
		{
			int sbmj = mat->GetSBM(j)->GetSBMID();
			if (sbmj == m_sbmId)
			{
				sbmid = j;
				break;
			}
		}
		if (sbmid == -1) return;
		
		double m_dt = mesh->GetFEModel()->GetCurrentStep()->m_dt;
		std::vector<FEMaterialPoint*> possible_ints;
		double total_elem = 0;
		FindIntInRadius(possible_ints, total_elem);
		double total_change = 0.0;
		// if the cell is not projecting on top of an integration point (i.e. too small) 
		// then project it's concentration to the nodes via the shape functions
		//SL TODO: Currently assigns just to the current element. Should instead find the ~8 closest integration points
		//		   regardless of element.
		if (possible_ints.size() == 0)
		{
			FindNodesInRadius(possible_ints, total_elem);
		}
		if (possible_ints.size() == 0) 
		{
			// set the concentration of all the integration points
			double H[FEElement::MAX_NODES];
			m_el->shape_fnc(H, m_q[0], m_q[1], m_q[2]);
			for (int i = 0; i < nint; ++i)
			{
				double H[FEElement::MAX_NODES];
				m_el->shape_fnc(H, m_q[0], m_q[1], m_q[2]);
				FEMaterialPoint& mp = *m_el->GetMaterialPoint(i);
				FESolutesMaterialPoint& pd = *(mp.ExtractData<FESolutesMaterialPoint>());
				double new_r = H[i] * m_rate * nint / Ve + pd.m_sbmr[sbmid];
				pd.m_sbmr[sbmid] = new_r;
				pd.m_sbmrp[sbmid] = pd.m_sbmr[sbmid];
				total_change += m_rate / nint;
			}
		}
		// else evenly distribute it among the integration points that the cell is on top of.
		//SL TODO: currently may lead to a little bias when neighboring elements are smaller resulting in
		//		   higher density of integration points
		else {
			int nint_in = possible_ints.size();
			for (auto iter = possible_ints.begin(); iter != possible_ints.end(); ++iter)
			{
				FEMaterialPoint& mp = **iter;
				FESolutesMaterialPoint& pd = *(mp.ExtractData<FESolutesMaterialPoint>());
				// scale by H so that the total integral over each element is consistent.
				double H = double(nint) / double(nint_in);
				double new_r = H * m_rate / Ve + pd.m_sbmr[sbmid];
				pd.m_sbmr[sbmid] = new_r;
				pd.m_sbmrp[sbmid] = pd.m_sbmr[sbmid];
				total_change += m_rate / (nint_in);
			}
		}	
		m_dC = -total_change / m_Vc;
		m_accumulate = false; // don't double count a point source
	}
}

void FESBMPointSource::FindIntInRadius(std::vector<FEMaterialPoint*> &possible_ints, double &total_elem) {
	
	// get element and set up buffers
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	std::unordered_set<FESolidElement*> visited;
	std::set<FESolidElement*> next;
	//std::vector<FEMaterialPoint*> possible_ints;
	visited.reserve(1000);
	possible_ints.reserve(500);

	//we will need to check the current element first
	next.insert(m_el);
	// create the element adjacency list.
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;
	// calculate the element volume
	auto mesh = dom->GetMesh();
	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(mesh);

	//while there are still elements to evaluate
	while (next.size()) {
		// get the element to be evaluated
		FESolidElement* cur = *next.begin();
		// remove the current element from the next buffer and add it to the visited buffer
		next.erase(next.begin());
		visited.insert(cur);
		// get the current element bounds
		std::vector<vec3d> cur_element_bounds;
		// add integration points within the radius
		bool int_flag = false;
		for (int i = 0; i < cur->GaussPoints(); i++)
		{
			FEMaterialPoint* mp = cur->GetMaterialPoint(i);
			vec3d disp = mp->m_r0 - m_pos;
			if (disp.norm() <= m_radius) {
				possible_ints.push_back(mp);
				int_flag = true;
			}
		}
		if (int_flag) 
		{
			total_elem++;
		}

		// Add neighboring element to the next buffer as long as they haven't been visited.
		// get the global ID of the current element
		int cur_id = cur->GetID()-1;
		// for each neighboring element
		for (int i = 0; i < EEL.NeighborSize(); i++)
		{
			if (EEL.Neighbor(cur_id, i))
			{
				// if that element has not been visited yet add it to the next list
				if (!visited.count(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i))))
				{
					next.insert(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i)));
				}
			}
		}
	}
}

void FESBMPointSource::FindNodesInRadius(std::vector<FEMaterialPoint*>& possible_ints, double& total_elem) {

	// get element and set up buffers
	m_el = dynamic_cast<FESolidElement*>(m_search.FindElement(m_pos, m_q));
	std::unordered_set<FESolidElement*> visited;
	std::set<FESolidElement*> next;
	//std::vector<FEMaterialPoint*> possible_ints;
	visited.reserve(1000);
	std::vector<FENode*> possible_nodes;
	possible_nodes.reserve(500);

	//we will need to check the current element first
	next.insert(m_el);
	// create the element adjacency list.
	FEDomain* dom = dynamic_cast<FEDomain*>(m_el->GetMeshPartition());
	FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom->GetMaterial());
	if (mat == nullptr) return;
	// calculate the element volume
	auto mesh = dom->GetMesh();
	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(mesh);

	//while there are still elements to evaluate
	while (next.size()) {
		// get the element to be evaluated
		FESolidElement* cur = *next.begin();
		// remove the current element from the next buffer and add it to the visited buffer
		next.erase(next.begin());
		visited.insert(cur);
		// get the current element bounds
		std::vector<vec3d> cur_element_bounds;
		// add integration points within the radius
		bool int_flag = false;
		for (int i = 0; i < cur->Nodes(); i++)
		{
			FENode* mn = &(mesh->Node(cur->m_node[i]));
			vec3d disp = mn->m_rt - m_pos;
			FEMaterialPoint* closest;
			if (disp.norm() <= m_radius) {
				// find the closest mp in the element
				double min_d = std::numeric_limits<double>::max();
				for (int j = 0; j < cur->GaussPoints(); j++)
				{
					double disp2 = (mn->m_rt - cur->GetMaterialPoint(j)->m_rt).norm();
					if (disp2 < min_d)
					{
						min_d = disp2;
						closest = cur->GetMaterialPoint(j);
					}
				}
				possible_ints.push_back(closest);
				int_flag = true;
			}
		}
		if (int_flag)
		{
			total_elem++;
		}

		// Add neighboring element to the next buffer as long as they haven't been visited.
		// get the global ID of the current element
		int cur_id = cur->GetID() - 1;
		// for each neighboring element
		for (int i = 0; i < EEL.NeighborSize(); i++)
		{
			if (EEL.Neighbor(cur_id, i))
			{
				// if that element has not been visited yet add it to the next list
				if (!visited.count(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i))))
				{
					next.insert(dynamic_cast<FESolidElement*>(EEL.Neighbor(cur_id, i)));
				}
			}
		}
	}
}

//void FESBMPointSource::ResetSBM()
//{
//	FEModel& fem = *GetFEModel();
//	FEMesh& mesh = fem.GetMesh();
//
//	for (int i = 0; i < mesh.Domains(); ++i)
//	{
//		FEDomain& dom = mesh.Domain(i);
//		int NE = dom.Elements();
//
//		FEMultiphasic* mat = dynamic_cast<FEMultiphasic*>(dom.GetMaterial());
//		if (mat)
//		{
//			// Make sure the material has the correct sbm
//			int sbmid = -1;
//			int sbms = mat->SBMs();
//			for (int j = 0; j<sbms; ++j)
//			{
//				int sbmj = mat->GetSBM(j)->GetSBMID();
//				if (sbmj == m_sbmId)
//				{
//					sbmid = j;
//					break;
//				}
//			}
//
//			if (sbmid != -1)
//			{
//				for (int j = 0; j < NE; ++j)
//				{
//					FEElement& el = dom.ElementRef(j);
//
//					// set the concentration of all the integration points
//					int nint = el.GaussPoints();
//					for (int k = 0; k < nint; ++k)
//					{
//						FEMaterialPoint* mp = el.GetMaterialPoint(k);
//						FESolutesMaterialPoint& pd = *(mp->ExtractData<FESolutesMaterialPoint>());
//						pd.m_sbmr[sbmid] = 0.0;
//						pd.m_sbmrp[sbmid] = 0.0;
//					}
//				}
//			}
//		}
//	}
//}
