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
#include "FEUncoupledMaterial.h"
#include "FEViscoElasticMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation uncoupled visco-elastic material
//
class FEUncoupledViscoElasticMaterial :	public FEUncoupledMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };
	
public:
	//! default constructor
	FEUncoupledViscoElasticMaterial(FEModel* pfem);

	// get the elastic base material
	FEUncoupledMaterial* GetBaseMaterial() { return m_pBase; }

	// set the elastic base material
	void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }

public:
	//! data initialization and checking
	bool Init() override;
	
	//! deviatoric stress function
	mat3ds DevStress(FEMaterialPoint& pt) override;
	
	//! deviatoric tangent function
	tens4ds DevTangent(FEMaterialPoint& pt) override;
	
	//! deviatoric strain energy density function
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! calculate exponent of right-stretch tensor in series spring
    bool SeriesStretchExponent(FEMaterialPoint& pt);

	//! returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
	
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	
private:
	FEUncoupledMaterial*	m_pBase;	//!< pointer to elastic solid material
	bool					m_binit;	//!< initialization flag
	
public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
