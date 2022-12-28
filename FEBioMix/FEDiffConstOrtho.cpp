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
#include "FEDiffConstOrtho.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEDiffConstOrtho, FESoluteDiffusivity)
	ADD_PARAMETER(m_free_diff,    FE_RANGE_GREATER(0.0)         , "free_diff")->setUnits(UNIT_DIFFUSIVITY)->setLongName("free diffusivity");
	ADD_PARAMETER(m_diff     , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "diff"     )->setUnits(UNIT_DIFFUSIVITY)->setLongName("diffusivity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstOrtho::FEDiffConstOrtho(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff[0] = m_diff[1] = m_diff[2] = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
bool FEDiffConstOrtho::Validate()
{
	if (FESoluteDiffusivity::Validate() == false) return false;
	if (m_free_diff < m_diff[0]) { feLogError("free_diff must be >= diff1"); return false; }
	if (m_free_diff < m_diff[1]) { feLogError("free_diff must be >= diff2"); return false; }
	if (m_free_diff < m_diff[2]) { feLogError("free_diff must be >= diff3"); return false; }
	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffConstOrtho::Free_Diffusivity(FEMaterialPoint& mp)
{
	return m_free_diff;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffConstOrtho::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor
mat3ds FEDiffConstOrtho::Diffusivity(FEMaterialPoint& mp)
{
	vec3d a0;				// texture direction in reference configuration
	mat3ds d(0,0,0,0,0,0);	// diffusion tensor

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// --- constant orthotropic diffusivity ---
	for (int i=0; i<3; i++) {	// Perform sum over all three texture directions
		
		// Copy the texture direction in the reference configuration to a0
		a0.x = Q[0][i]; a0.y = Q[1][i]; a0.z = Q[2][i];
		
		// Evaluate the texture tensor in the current configuration
		d += dyad(a0)*m_diff[i];
	}
	
	return d;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4dmm FEDiffConstOrtho::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
	tens4dmm D;
	D.zero();
	return D;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffConstOrtho::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
	mat3ds d;
	d.zero();
	return d;
}
