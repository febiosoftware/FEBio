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
#include "FEDiffConstIso.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEDiffConstIso, FESoluteDiffusivity)
	ADD_PARAMETER(m_free_diff, FE_RANGE_GREATER_OR_EQUAL(0.0), "free_diff")->setUnits(UNIT_DIFFUSIVITY)->setLongName("free diffusivity");
	ADD_PARAMETER(m_diff     , FE_RANGE_GREATER_OR_EQUAL(0.0), "diff"     )->setUnits(UNIT_DIFFUSIVITY)->setLongName("diffusivity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstIso::FEDiffConstIso(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff = 1;
}

//-----------------------------------------------------------------------------
//! Validation
bool FEDiffConstIso::Validate()
{
	if (FESoluteDiffusivity::Validate() == false) return false;
	if (m_free_diff < m_diff) {
		feLogError("free_diff must be >= diff");
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffConstIso::Free_Diffusivity(FEMaterialPoint& mp)
{
	return m_free_diff;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffConstIso::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor
mat3ds FEDiffConstIso::Diffusivity(FEMaterialPoint& mp)
{
	// --- constant isotropic diffusivity ---
	
	return mat3dd(m_diff);
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4dmm FEDiffConstIso::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
	tens4dmm D;
	D.zero();
	return D;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffConstIso::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
	mat3ds d;
	d.zero();
	return d;
}
