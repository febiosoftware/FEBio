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
#include "FEElasticMaterial.h"


class FEPRLig : public FEElasticMaterial
{
public:
	FEPRLig(FEModel* pfem);

public:
	
	double	m_c1;	 //!< fiber constant c1
	double	m_c2;	 //!< fiber constant c2 
	double  m_u;     //!< Lame Coefficient mu of the matrix
	double	m_m;	 //!< Poisson's ratio slope 
	double	m_v0;	 //!< initial Poisson's ratio
	double	m_k;	 //!< Penalty for the volumetric strain energy


public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
