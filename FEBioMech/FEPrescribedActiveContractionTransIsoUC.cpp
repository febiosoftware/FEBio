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
#include "FEPrescribedActiveContractionTransIsoUC.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPrescribedActiveContractionTransIsoUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_T0 , "T0"   );

    ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionTransIsoUC::FEPrescribedActiveContractionTransIsoUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_T0 = 0.0;
}

//-----------------------------------------------------------------------------
bool FEPrescribedActiveContractionTransIsoUC::Validate()
{
	if (FEUncoupledMaterial::Validate() == false) return false;

    // fiber direction in local coordinate system (reference configuration)
    m_n0.x = 1;
    m_n0.y = 0;
    m_n0.z = 0;
	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedActiveContractionTransIsoUC::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_n0;
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
	double T0 = m_T0(mp);
    mat3ds s = (b-N)*(T0/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionTransIsoUC::DevTangent(FEMaterialPoint &mp)
{
    return tens4ds(0.0);
}
