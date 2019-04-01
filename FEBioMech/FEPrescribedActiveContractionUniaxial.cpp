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
#include "FEPrescribedActiveContractionUniaxial.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPrescribedActiveContractionUniaxial, FEElasticMaterial)
	ADD_PARAMETER(m_T0 , "T0"   );
	ADD_PARAMETER(m_thd, "theta");
	ADD_PARAMETER(m_phd, "phi"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedActiveContractionUniaxial::FEPrescribedActiveContractionUniaxial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_thd = 0;
    m_phd = 90;
}

//-----------------------------------------------------------------------------
bool FEPrescribedActiveContractionUniaxial::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;

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
void FEPrescribedActiveContractionUniaxial::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_n0;
}

//-----------------------------------------------------------------------------
mat3ds FEPrescribedActiveContractionUniaxial::Stress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

    // evaluate fiber direction in global coordinate system
    vec3d n0 = Q*m_n0;

    // evaluate the deformed fiber direction
    vec3d nt = F*n0;
    mat3ds N = dyad(nt);
    
    // evaluate the active stress
    mat3ds s = N*(m_T0/J);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPrescribedActiveContractionUniaxial::Tangent(FEMaterialPoint &mp)
{
    tens4ds c;
    c.zero();

    return c;
}
