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



#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is isotropic in the reference state, but exhibits
// strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEBIOMIX_API FEPermRefIso :	public FEHydraulicPermeability
{
public:
	//! constructor
	FEPermRefIso(FEModel* pfem);
		
	//! permeability
	mat3ds Permeability(FEMaterialPoint& pt) override;
		
	//! Tangent of permeability
	tens4dmm Tangent_Permeability_Strain(FEMaterialPoint& mp) override;
		
	//! data initialization and checking
	bool Validate() override;
		
public:
	FEParamDouble m_perm0;		//!< permeability for I term
	FEParamDouble m_perm1;		//!< permeability for b term
	FEParamDouble m_perm2;		//!< permeability for b^2 term
	double	m_M;			//!< nonlinear exponential coefficient
	double	m_alpha;		//!< nonlinear power exponent
		
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
