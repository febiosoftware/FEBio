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
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! This is a neo-Hookean material whose Young's modulus is evaluated from the density
//! according to the power-law relation proposed by Carter and Hayes for trabecular bone

class FEBIOMIX_API FECarterHayes : public FEElasticMaterial, public FERemodelingInterface
{
public:
	FECarterHayes(FEModel* pfem) : FEElasticMaterial(pfem) { m_E0 = 0; m_rho0 = 1; m_sbm = -1; m_lsbm = -1; m_g = 0; }
	
public: // parameters
	double	m_E0;	//!< Young's modulus at reference sbm density
	double	m_g;	//!< gamma exponent for calculation of Young's modulus
	double	m_v;	//!< prescribed Poisson's ratio
    double  m_rho0; //!< reference sbm density
	int		m_sbm;	//!< global id of solid-bound molecule

protected:
	int		m_lsbm;	//!< local id of solid-bound molecule

public:
	//! data initialization and checking
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! Create material point data
	FEMaterialPointData* CreateMaterialPointData() override;
	
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override { return StrainEnergy(pt); }
    
public: // --- remodeling interface ---

	//! calculate strain energy density at material point
	double StrainEnergy(FEMaterialPoint& pt) override;

	//! calculate tangent of strain energy density with mass density
	double Tangent_SE_Density(FEMaterialPoint& pt) override;
	
	//! calculate tangent of stress with mass density
	mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) override;

	//! return Young's modulus
	double YoungModulus(double rhor) { return m_E0*pow(rhor/m_rho0, m_g);}
	
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
