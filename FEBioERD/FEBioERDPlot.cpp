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
#include "FEBioERDPlot.h"
#include <FEBioMix/FEBioMixPlot.h>
#include <FEBioMix/FEBiphasicSolidDomain.h>
#include <FEBioMix/FEBiphasicSoluteSolidDomain.h>
#include <FEBioMix/FEMultiphasicSolidDomain.h>
#include "FEBioERD/FEElasticReactionDiffusionSolidDomain.h"
#include <FEBioMech/FEElasticSolidDomain.h>
#include <FEBioMix/FEBiphasic.h>
#include <FEBioMix/FEBiphasicSolute.h>
#include <FEBioMix/FEMultiphasic.h>
#include "FEBioERD/FEElasticReactionDiffusion.h"
#include "FEBioMech/FEElasticMixture.h"
#include <FECore/FEModel.h>
#include <FECore/writeplot.h>
#include <FECore/FEEdgeList.h>
#include <FECore/mathalg.h>

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotGrowthRatio::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* pmf = kg->GetGrowthMaterial();
	//FEVolumeGrowth* vg = dynamic_cast<FEVolumeGrowth*>(pmf);
	if (pmf == nullptr) return false;
	//if (vg == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double g = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& pt = *el.GetMaterialPoint(j);
			FEKinematicMaterialPointERD& kp = *pt.ExtractData<FEKinematicMaterialPointERD>();
			g += kp.m_theta;
		}
		g /= (double)el.GaussPoints();
		a << g;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotGrowthTensor::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* pmf = kg->GetGrowthMaterial();
	//FEVolumeGrowth* vg = dynamic_cast<FEVolumeGrowth*>(pmf);
	if (pmf == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FESolidElement& se = dynamic_cast<FESolidElement&>(dom.ElementRef(i));

		std::vector<mat3ds> SPDs_gausspts;
		// determine shape function value for the local position
		// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
		for (int j = 0; j < se.GaussPoints(); ++j) {
			FEMaterialPoint& pt = *se.GetMaterialPoint(j);
			vec3d n0 = pmf->m_fiber.unitVector(pt);
			mat3d gt = pmf->GrowthTensor(pt, n0);
			mat3ds gts = gt.sym();
			SPDs_gausspts.push_back(gts);
		}

		// array containing the SPD for each node in the element
		mat3ds SPDs_nodes[FESolidElement::MAX_NODES];
		// array for the shape function values
		double H[FESolidElement::MAX_NODES];
		// project the spds from integration points to the nodes
		se.project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
		double centerv;
		switch (se.Type())
		{
		case FE_Element_Type::FE_HEX8G8:
		case FE_Element_Type::FE_HEX8RI:
		case FE_Element_Type::FE_HEX8G1:
		case FE_Element_Type::FE_HEX20G8:
		case FE_Element_Type::FE_HEX20G27:
		case FE_Element_Type::FE_HEX27G27:
		{
			centerv = 0.0;
		}
		case FE_Element_Type::FE_TET4G1:
		case FE_Element_Type::FE_TET4G4:
		case FE_Element_Type::FE_TET5G4:
		case FE_Element_Type::FE_TET10G1:
		case FE_Element_Type::FE_TET10G4:
		case FE_Element_Type::FE_TET10G8:
		case FE_Element_Type::FE_TET10GL11:
		case FE_Element_Type::FE_TET10G4RI1:
		case FE_Element_Type::FE_TET10G8RI4:
		case FE_Element_Type::FE_TET15G4:
		case FE_Element_Type::FE_TET15G8:
		case FE_Element_Type::FE_TET15G11:
		case FE_Element_Type::FE_TET15G15:
		case FE_Element_Type::FE_TET15G15RI4:
		case FE_Element_Type::FE_TET20G15:
		case FE_Element_Type::FE_PENTA6G6:
		case FE_Element_Type::FE_PENTA15G8:
		case FE_Element_Type::FE_PENTA15G21:
		case FE_Element_Type::FE_PYRA5G8:
		case FE_Element_Type::FE_PYRA13G8:
		{
			centerv = 0.5;
		}
		default:
		{
			break;
		}
		}
		se.shape_fnc(H, centerv, centerv, centerv);
		mat3ds g = weightedAverageStructureTensor(SPDs_nodes, H, se.Nodes());
		a << g;
	}
	return true;
}

//=============================================================================
// Stress traces

bool FEPlotTraceStresses::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		mat3ds s = mat3ds(0.0);
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint* mp = el.GetMaterialPoint(j);
			const FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
			s += ep->m_s;
		}
		s /= (double)el.GaussPoints();
		a << s.tr();
	}
	return true;
};

bool FEPlotGrowthElasticDeformationGradient::Save(FEDomain& dom, FEDataStream& a)
{
	//// For each element get the growth tensor and then solve the average value.

	FEKinematicGrowthRateDependent* kgm = dom.GetMaterial()->ExtractProperty<FEKinematicGrowthRateDependent>();
	if (kgm == nullptr) return false;

	writeAverageElementValue<mat3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEKinematicMaterialPointERD& pt = *mp.ExtractData<FEKinematicMaterialPointERD>();
		return pt.m_Fe;
		});

	return true;
};

bool FEPlotGrowthDeformationGradient::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* gmat = kg->GetGrowthMaterial();
	//FEVolumeGrowth* vg = dynamic_cast<FEVolumeGrowth*>(pmf);
	if (gmat == nullptr) return false;
	//if (vg == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		mat3d g = mat3dd(1.0);
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEKinematicMaterialPointERD& kp = *mp.ExtractData<FEKinematicMaterialPointERD>();
			g += kp.m_Fg;
		}
		g /= (double)el.GaussPoints();
		a << g;
	}
	return true;
}

bool FEPlotJacobian::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double J = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
			J += ep.m_J;
		}
		J /= (double)el.GaussPoints();
		a << J;
	}
	return true;
}

bool FEPlotGrowthElasticJacobian::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* gmat = kg->GetGrowthMaterial();
	//FEVolumeGrowth* vg = dynamic_cast<FEVolumeGrowth*>(pmf);
	if (gmat == nullptr) return false;
	//if (vg == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double J_e = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEKinematicMaterialPointERD& kp = *mp.ExtractData<FEKinematicMaterialPointERD>();
			J_e += kp.m_Je;
		}
		J_e /= (double)el.GaussPoints();
		a << J_e;
	}
	return true;
}

bool FEPlotGrowthJacobian::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* gmat = kg->GetGrowthMaterial();
	if (gmat == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double J_g = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEKinematicMaterialPointERD& kp = *mp.ExtractData<FEKinematicMaterialPointERD>();
			J_g += kp.m_Jg;
		}
		J_g /= (double)el.GaussPoints();
		a << J_g;
	}
	return true;
}

bool FEPlotGrowthK::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* gmat = kg->GetGrowthMaterial();
	if (gmat == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double k_theta = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			k_theta += gmat->ActivationFunction(mp);
		}
		k_theta /= (double)el.GaussPoints();
		a << k_theta;
	}
	return true;
}

bool FEPlotGrowthPhi::Save(FEDomain& dom, FEDataStream& a)
{
	// Try to get the kinematic growth material.
	FEMaterial* mat = dom.GetMaterial();
	FEKinematicGrowthRateDependent* kg = mat->ExtractProperty<FEKinematicGrowthRateDependent>();
	// Check if we were successful.
	if (kg == nullptr) return false;
	FEGrowthTensorERD* gmat = kg->GetGrowthMaterial();
	if (gmat == nullptr) return false;
	// For each element get the growth tensor and then solve the average value.
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		double phi = 0.0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			phi += gmat->EnvironmentalFunction(mp);
		}
		phi /= (double)el.GaussPoints();
		a << phi;
	}
	return true;
}

FEPlotActualSoluteConcentrationERD::FEPlotActualSoluteConcentrationERD(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i < ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
	SetUnits(UNIT_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotActualSoluteConcentrationERD::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticReactionDiffusionInterface* pm = dynamic_cast<FEElasticReactionDiffusionInterface*>(dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i < (int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] < 0) negs++;
	}
	if (negs == nsols) return false;

	// loop over all elements
	int N = dom.Elements();
	for (int i = 0; i < N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k = 0; k < nsols; ++k)
		{
			int nsid = lid[k];
			if (nsid == -1) a << 0.f;
			else
			{
				// calculate average concentration
				double ew = 0;
				for (int j = 0; j < el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					ew += pm->GetActualSoluteConcentration(mp, nsid);
				}
				ew /= el.GaussPoints();
				a << ew;
			}
		}

	}
	return true;
}

//=================================================================================================
// FEPlotEffectiveSoluteConcentration
//=================================================================================================

FEPlotEffectiveSoluteConcentrationERD::FEPlotEffectiveSoluteConcentrationERD(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_NODE)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i < ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
	SetUnits(UNIT_CONCENTRATION);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveSoluteConcentrationERD::Save(FEDomain& dom, FEDataStream& a)
{
	// get the dof
	DOFS& dofs = GetFEModel()->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	if (nsol == -1) return false;

	// get the start index
	const int dof_C = GetFEModel()->GetDOFIndex("concentration", 0);

	FEElasticReactionDiffusionInterface* pm = dynamic_cast<FEElasticReactionDiffusionInterface*>(dom.GetMaterial());
	if (pm == 0) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int negs = 0;
	for (int i = 0; i < (int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] < 0) negs++;
	}
	if (negs == nsol) return false;

	// save the concentrations
	int N = dom.Nodes();
	for (int i = 0; i < N; ++i)
	{
		FENode& node = dom.Node(i);
		for (int j = 0; j < nsol; ++j)
		{
			double c = (lid[j] >= 0 ? node.get(dof_C + j) : 0.0);
			a << c;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSoluteFluxERD::FEPlotSoluteFluxERD(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY_VEC3F, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i < ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
			m_sol.push_back(ps->GetID());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);
	SetUnits(UNIT_MOLAR_FLUX);
}

//-----------------------------------------------------------------------------
bool FEPlotSoluteFluxERD::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticReactionDiffusionInterface* pm = dynamic_cast<FEElasticReactionDiffusionInterface*>(dom.GetMaterial());
	if ((pm == 0) || (pm->Solutes() == 0)) return false;

	// figure out the local solute IDs. This depends on the material
	int nsols = (int)m_sol.size();
	vector<int> lid(nsols, -1);
	int nsc = 0;
	for (int i = 0; i < (int)m_sol.size(); ++i)
	{
		lid[i] = pm->FindLocalSoluteID(m_sol[i]);
		if (lid[i] != -1) nsc++;
	}
	if (nsc == 0) return false;

	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k = 0; k < nsols; ++k)
		{
			int nsid = lid[k];
			if (nsid == -1) a << vec3d(0, 0, 0);
			else
			{
				// calculate average flux
				vec3d ew = vec3d(0, 0, 0);
				for (int j = 0; j < el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					ew += pm->GetSoluteFlux(mp, nsid);
				}

				ew /= el.GaussPoints();

				a << ew;
			}
		}
	}
	return true;
}