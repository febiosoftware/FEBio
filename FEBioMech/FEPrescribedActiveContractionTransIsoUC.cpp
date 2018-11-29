//
//  FEPrescribedActiveContractionTransIsoUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEPrescribedActiveContractionTransIsoUC.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPrescribedActiveContractionTransIsoUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_T0 , "T0"   );
	ADD_PARAMETER(m_thd, "theta");
	ADD_PARAMETER(m_phd, "phi"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionTransIsoUC::FEPrescribedActiveContractionTransIsoUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_thd = 0;
    m_phd = 90;
}

//-----------------------------------------------------------------------------
bool FEPrescribedActiveContractionTransIsoUC::Validate()
{
	if (FEUncoupledMaterial::Validate() == false) return false;
    // convert angles from degrees to radians
    double pi = 4*atan(1.0);
    double the = m_thd*pi/180.;
    double phi = m_phd*pi/180.;
    // fiber direction in local coordinate system (reference configuration)
    m_n0.x = cos(the)*sin(phi);
    m_n0.y = sin(the)*sin(phi);
    m_n0.z = cos(phi);
	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedActiveContractionTransIsoUC::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			ar << m_n0;
		}
		else
		{
			ar >> m_n0;
		}
	}
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedActiveContractionTransIsoUC::DevStress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F;
    mat3ds b = pt.LeftCauchyGreen();
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // evaluate fiber direction in global coordinate system
    vec3d n0 = Q*m_n0;
    
    // evaluate the deformed fiber direction
    vec3d nt = F*n0;
    mat3ds N = dyad(nt);
    
    // evaluate the active stress
    mat3ds s = (b-N)*(m_T0/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionTransIsoUC::DevTangent(FEMaterialPoint &mp)
{
    return tens4ds(0.0);
}
