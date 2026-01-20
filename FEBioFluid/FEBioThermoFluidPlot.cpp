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
#include "FEBioThermoFluidPlot.h"
#include "FEFluidMaterial.h"
#include "FEThermoFluid.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/writeplot.h>
#include <FECore/FEModel.h>

bool FEPlotFluidTemperature::Save(FEDomain& dom, FEDataStream& a)
{
	FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		return pfluid->Temperature(const_cast<FEMaterialPoint&>(mp));
		});
	return true;
}

//! Store the nodal temperatures
bool FEPlotNodalFluidTemperature::Save(FEMesh& m, FEDataStream& a)
{
	// get the dilatation dof index
	int dof_T = GetFEModel()->GetDOFIndex("T");
	if (dof_T < 0) return false;

	// loop over all nodes
	writeNodalValues<double>(m, a, [=](const FENode& node) {
		return node.get(dof_T);
		});
	return true;
}

bool FEPlotFluidPressureTangentTemperature::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->Tangent_Temperature(mp_noconst);
		});

	return true;
}

bool FEPlotFluidHeatFlux::Save(FEDomain& dom, FEDataStream& a)
{
	FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	// write solid element data
	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEThermoFluidMaterialPoint* ppt = mp.ExtractData<FEThermoFluidMaterialPoint>();
		return (ppt ? ppt->m_q : vec3d(0.));
		});

	return true;
}

bool FEPlotFluidRelativeThermalPecletNumber::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [&pfluid](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* fpt = mp.ExtractData<FEFluidMaterialPoint>();
		//        const FEElasticMaterialPoint* ept = mp.ExtractData<FEElasticMaterialPoint>();
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		double cp = pfluid->GetElastic()->IsobaricSpecificHeatCapacity(mp_noconst);
		double K = pfluid->GetConduct()->ThermalConductivity(mp_noconst);
		double rho = pfluid->Density(mp_noconst);
		vec3d v(0, 0, 0);
		//        if (ept) v = ept->m_v;
		return (fpt->m_vft - v).Length() * rho * cp / K;
		});

	return true;
}

bool FEPlotFluidThermalConductivity::Save(FEDomain& dom, FEDataStream& a)
{
	FEFluidThermalConductivity* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidThermalConductivity>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->ThermalConductivity(mp_noconst);
		});

	return true;
}

bool FEPlotFluidIsochoricSpecificHeatCapacity::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->IsochoricSpecificHeatCapacity(mp_noconst);
		});

	return true;
}

bool FEPlotFluidIsobaricSpecificHeatCapacity::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->IsobaricSpecificHeatCapacity(mp_noconst);
		});

	return true;
}
