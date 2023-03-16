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
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
//! Base class for uncoupled hyperelastic material formulations.

//! In FEBio, for uncoupled materials it is assumed that the strain energy function 
//! is a sum of two terms, a deviatoric strain energy term, which only depends on the 
//! deviatoric right Cauchy-Green tensor C_tilde, and a volumetric term which only 
//! depends on J, the determinant of the deformation gradient. The total Cauchy stress 
//! is therefore also a sum of two contributions, namely the deviatoric Cauchy stress
//! and a pressure term where p = dU/dJ. 

//! One of the main motivations for the alternative material interface is that some
//! finite element implementations can take advantage of the uncoupling. For example,
//! the three-field formulation integrates the two different material terms differently
//! in order to avoid locking problems for (nearly-) incompressible materials. 

//! When implementing a new material derived from this base class, the developer needs
//! to provide the deviatoric stress function and the pressure function as well as their
//! derivatives. 

class FEBIOMECH_API FEUncoupledMaterial : public FEElasticMaterial
{
public:
	//! constructor
	FEUncoupledMaterial(FEModel* pfem);

	//! initialization
	bool Init() override;

public:

//----------------->
	// The following functions need to be overloaded
	// for each material derived from this class.

	//! Deviatoric Cauchy stress
	virtual mat3ds DevStress(FEMaterialPoint& mp) = 0;

	//! Deviatoric spatial Tangent
	virtual tens4ds DevTangent(FEMaterialPoint& mp) = 0;

	//! Deviatoric strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp) { return 0; }
    
public:
    virtual double StrongBondDevSED(FEMaterialPoint& pt) { return DevStrainEnergyDensity(pt); }
    virtual double WeakBondDevSED(FEMaterialPoint& pt) { return 0; }

public:
    // TODO: removing virtual from the following 3 functions causes changes
    // to the convergence criteria on macOS, despite these functions not
    // being overridden anywhere.
	//! strain energy density U(J)
    virtual double U(double J) {
        switch (m_npmodel) {
            case 0: return 0.5*m_K*pow(log(J),2); break;    // FEBio default
            case 1: return 0.25*m_K*(J*J - 2.0*log(J) - 1.0); break;    // NIKE3D's Ogden material
            case 2: return 0.5*m_K*(J-1)*(J-1); break;      // ABAQUS
            case 3: return 0.5*m_K*((J*J-1)/2-log(J)); break;      // ABAQUS - GOH
            default: { assert(false); return 0; }
        }
    }
	//! pressure, i.e. first derivative of U(J)
	virtual double UJ(double J) {
        switch (m_npmodel) {
            case 0: return m_K*log(J)/J; break;
            case 1: return 0.5*m_K*(J - 1.0/J); break;
            case 2: return m_K*(J-1); break;
            case 3: return 0.5*m_K*(J-1.0/J); break;
			default: { assert(false); return 0; }
		}
    }

	//! second derivative of U(J) 
	virtual double UJJ(double J) {
        switch (m_npmodel) {
            case 0: return m_K*(1-log(J))/(J*J); break;
            case 1: return 0.5*m_K*(1 + 1.0/(J*J)); break;
            case 2: return m_K; break;
            case 3: return 0.5*m_K*(1+1.0/(J*J)); break;
			default: { assert(false); return 0; }
		}
    }

public:
	// incompressibility constraint fnc and derivs
	double h  (double J) { return log(J); }
	double hp (double J) { return 1.0 / J; }
	double hpp(double J) { return -1.0 / (J*J); }

public:
	//! total Cauchy stress (do not overload!)
	mat3ds Stress(FEMaterialPoint& mp) final;

	//! total spatial tangent (do not overload!)
	tens4ds Tangent(FEMaterialPoint& mp) final;

	//! calculate strain energy (do not overload!)
	double StrainEnergyDensity(FEMaterialPoint& pt) final;
    double StrongBondSED(FEMaterialPoint& pt) final;
    double WeakBondSED(FEMaterialPoint& pt) final;

	// Create material point data
	FEMaterialPointData* CreateMaterialPointData() override;
    
public:
	double	m_K;			//!< bulk modulus
	int     m_npmodel;      //!< pressure model for U(J)

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEUncoupledMaterial)
};
