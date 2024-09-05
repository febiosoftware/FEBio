/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2024 University of Utah, The Trustees of Columbia University in
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
#include <FEBioMix/FEBiphasic.h>
#include <FEBioMix/FESolutesMaterialPoint.h>
#include <FEBioMix/FESolute.h>
#include "FEChemicalReactionERD.h"
#include "FEElasticReactionDiffusionInterface.h"
#include <FECore/FEModelParam.h>
#include "FEKinematicGrowthRateDependent.h"
#include "FEGrowthTensorERD.h"

//-----------------------------------------------------------------------------
//! Base class for elastic reaction diffusion materials.

class FEBIOERD_API FEElasticReactionDiffusion : public FEMaterial, public FEElasticReactionDiffusionInterface_T<FESolutesMaterialPoint>
{
public:
	//! constructor
	FEElasticReactionDiffusion(FEModel* pfem);

	//! initialization
	bool Init() override;
    
    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	// return elastic material component
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }

	FEKinematicGrowthRateDependent* GetKinematicGrowthMaterial() { return dynamic_cast<FEKinematicGrowthRateDependent*>(m_pSolid); }

public:
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
		
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt, const int sol);
	
    //! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int sol);

	//! porosity
	double Porosity(FEMaterialPoint& pt);

	//! solute density
	double SoluteDensity(const int sol) { return m_pSolute[sol]->Density(); }

	//! solute molar mass
	double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }

	//! Add a chemical reaction
	void AddChemicalReaction(FEChemicalReactionERD* pcr);

	// solute interface
public:
	int Solutes() override { return (int)m_pSolute.size(); }
	FESolute* GetSolute(int i) override { return m_pSolute[i]; }

public:
	FEElasticMaterial* GetSolid() { return m_pSolid; }
	FEChemicalReactionERD* GetReaction(int i) { return m_pReact[i]; }

public:
	int Reactions() { return (int)m_pReact.size(); }

public: // parameters
	FEParamDouble       m_phi0;     //!< solid volume fraction in reference configuration

public:
	double	m_Rgas;			//!< universal gas constant
	double	m_Tabs;			//!< absolute temperature

protected:
	// material properties
	FEElasticMaterial* m_pSolid;		//!< pointer to elastic solid material
	std::vector<FESolute*>				m_pSolute;		//!< pointer to solute materials
	std::vector<FEChemicalReactionERD*>	m_pReact;		//!< pointer to chemical reactions
	
	DECLARE_FECORE_CLASS();
};