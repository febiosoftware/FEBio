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
#include <FEBioFluid/FEFluidMaterial.h>
#include "FEThermoFluid.h"
#include "FEThermoElasticFluid.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FECore/writeplot.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
bool FEPlotFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEThermoFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluidMaterial>();
    if (pfluid == 0) return false;

    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
        return (ppt ? ppt->m_vft : vec3d(0.));
    });

    return true;
}

bool FEPlotFluidTemperature::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluid>();
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

//-----------------------------------------------------------------------------
bool FEPlotFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
    FEThermoFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluidMaterial>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
        return (pt ? pt->m_pf : 0.0);
    });
    return true;
}

bool FEPlotFluidPressureTangentTemperature::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->Tangent_Temperature(mp_noconst);
		});

	return true;
}

bool FEPlotFluidHeatFlux::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
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
		double K = pfluid->GetHeatFlux()->Conductivity(mp_noconst);
		double rho = pfluid->Density(mp_noconst);
		vec3d v(0, 0, 0);
		//        if (ept) v = ept->m_v;
		return (fpt->m_vft - v).Length() * rho * cp / K;
		});

	return true;
}

bool FEPlotFluidThermalConductivity::Save(FEDomain& dom, FEDataStream& a)
{
    FEThermoFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluid>();
    if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->GetHeatFlux()->Conductivity(mp_noconst);
		});

	return true;
}

bool FEPlotFluidIsochoricSpecificHeatCapacity::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->IsochoricSpecificHeatCapacity(mp_noconst);
		});

	return true;
}

bool FEPlotFluidIsobaricSpecificHeatCapacity::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->IsobaricSpecificHeatCapacity(mp_noconst);
		});

	return true;
}

bool FEPlotFluidSpecificFreeEnthalpy::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == nullptr) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->SpecificFreeEnthalpy(mp_noconst);
		});

	return true;
}

bool FEPlotFluidSpecificGaugeEnthalpy::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->SpecificGaugeEnthalpy(mp_noconst);
		});

	return true;
}

bool FEPlotFluidSpecificInternalEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->SpecificInternalEnergy(mp_noconst);
		});

	return true;
}

bool FEPlotFluidSpecificEntropy::Save(FEDomain& dom, FEDataStream& a)
{
	FEThermoElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoElasticFluid>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->SpecificEntropy(mp_noconst);
		});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidHeatSupplyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEThermoFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEThermoFluidMaterial>();
    if (pfluid == 0) return false;
    FEThermoFluid* ptf = dom.GetMaterial()->ExtractProperty<FEThermoFluid>();

    // write solid element data
    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        FEFluidMaterialPoint* ppt = (mp_noconst.ExtractData<FEFluidMaterialPoint>());
        FEThermoFluidMaterialPoint* tfp = (mp_noconst.ExtractData<FEThermoFluidMaterialPoint>());
        double T = tfp->m_T + pfluid->GetViscous()->m_Tr;
        double dpdT = ptf->GetElastic()->Tangent_Temperature(mp_noconst);
        double r = (pfluid->GetViscous()->Stress(mp_noconst)*ppt->m_Lf).trace()
        -T*dpdT*ppt->m_Lf.trace();
        double rho = pfluid->Density(mp_noconst);
        return (tfp ? rho*r : 0.0);
    });
    
    return true;
}

