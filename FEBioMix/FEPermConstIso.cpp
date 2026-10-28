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
#include "FEPermConstIso.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEPermConstIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm")->setUnits(UNIT_PERMEABILITY);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermConstIso::FEPermConstIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm = 1;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermConstIso::Permeability(FEMaterialPoint& mp)
{
	// --- constant isotropic permeability ---
	
	return mat3dd(m_perm(mp));
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4dmm FEPermConstIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	// Even though the spatial permeability k is constant, the tangent is the push-forward
	// of 2*dK/dC, where K = J F^-1 k F^-T is the referential permeability. Therefore it
	// includes the geometric terms k (I x I - 2 I o I), consistent with the other
	// permeability materials (e.g. Holmes-Mow with k' = 0).
	mat3dd I(1);
	double k = m_perm(mp);
	return dyad1mm(I, mat3ds(I*k)) - dyad4s(I)*(2*k);
}
