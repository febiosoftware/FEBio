/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEMultiphasic.h"
#include "FESoluteInterface.h"

//-----------------------------------------------------------------------------
//! Base class for triphasic materials.

class FEBIOMIX_API FETriphasic : public FEMaterial, public FESoluteInterface
{
public:
	FETriphasic(FEModel* pfem);

	// initialization
	bool Init() override;

	// serialization
	void Serialize(DumpStream& ar) override;
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override
	{ 
		return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
	}

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }

public:

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	
	//! calculate fluid (solvent) flux
	vec3d FluidFlux(FEMaterialPoint& pt);
	
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt, const int ion);
	
	//! actual fluid pressure (as opposed to effective pressure)
	double Pressure(FEMaterialPoint& pt);
	
	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int ion);
	
	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fixed charge density
	double FixedChargeDensity(FEMaterialPoint& pt);
	
	//! electric potential
	double ElectricPotential(FEMaterialPoint& pt, const bool eform=false);
	
	//! current density
	vec3d CurrentDensity(FEMaterialPoint& pt);
	
    //! partition coefficient
    double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
    
    //! partition coefficient derivatives
    void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                       vector<double>& dkdJ,
                                       vector< vector<double> >& dkdc);
    //! fluid density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity(const int ion) { return m_pSolute[ion]->Density(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int ion) { return m_pSolute[ion]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int ion) { return m_pSolute[ion]->ChargeNumber(); }
	
	//! Add a solute component
	void AddSolute(FESolute* ps);

// solute interface
public:
	int Solutes() override { return (int)m_pSolute.size(); }
	FESolute* GetSolute(int i) override { return m_pSolute[i]; }

public:
    FEElasticMaterial*			GetSolid()				{ return m_pSolid; }
    FEHydraulicPermeability*	GetPermeability()		{ return m_pPerm;  }
    FEOsmoticCoefficient*		GetOsmoticCoefficient() { return m_pOsmC;  }
    
public: // material parameters
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_rhoTw;		//!< true fluid density
	double						m_cFr;			//!< fixed charge density in reference configurations
	double						m_penalty;		//!< penalty for enforcing electroneutrality

public:
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	double						m_Fc;			//!< Faraday's constant

public: // material properties
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	std::vector<FESolute*>		m_pSolute;		//!< pointer to solute materials

	DECLARE_FECORE_CLASS();
};
