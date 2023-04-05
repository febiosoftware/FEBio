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
#include <FECore/FEMaterial.h>
#include <FEBioMix/FESolute.h>
#include <FEBioMix/FESoluteInterface.h>
#include <FEBioMix/FEOsmoticCoefficient.h>
#include <FEBioMix/FEChemicalReaction.h>
#include "febiofluid_api.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
class FEBIOFLUID_API FESolutesMaterial : public FEMaterial, public FESoluteInterface
{
public:
	class Point : public FEMaterialPointData
	{
	public:
		//! constructor
		Point(FEMaterialPointData* pt);

		//! create a shallow copy
		FEMaterialPointData* Copy();

		//! data serialization
		void Serialize(DumpStream& ar);

		//! Data initialization
		void Init();

        //! Osmolarity        
        double Osmolarity() const;

	public:
		vec3d	m_vft;		// fluid velocity at integration point
		double	m_JfdotoJf;		// divergence of fluid velocity

		// solutes material data
		int                 m_nsol;     //!< number of solutes
		vector<double>      m_c;        //!< solute concentration
        vector<double>      m_ca;        //!< effective solute concentration
		vector<vec3d>       m_gradc;    //!< spatial gradient of solute concentration
		vector<vec3d>       m_j;        //!< solute molar flux
		vector<double>      m_cdot;     //!< material time derivative of solute concentration following fluid
        vector<double>    m_k;        //!< solute partition coefficient
        vector<double>    m_dkdJ;        //!< 1st deriv of m_k with strain (J)
        vector< vector<double> >    m_dkdc;            //!< 1st deriv of m_k with effective concentration
	};

public:
	FESolutesMaterial(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

	//! performs initialization
	bool Init() override;

    //! data serialization
    void Serialize(DumpStream& ar) override;

public:
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt, const int sol);

	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int sol);
    
    //! actual concentration (as opposed to effective concentration)
    double ConcentrationActual(FEMaterialPoint& pt, const int sol);
    
    //! actual fluid pressure (as opposed to effective pressure)
    double PressureActual(FEMaterialPoint& pt);
    
    //! partition coefficient
    double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
    
    //! partition coefficients and their derivatives
    void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                       vector<double>& dkdJ,
                                       vector< vector<double> >& dkdc);
    
    //! solute density
    double SoluteDensity(const int sol) { return m_pSolute[sol]->Density(); }
    
    //! solute molar mass
    double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }
    
    //! Add a chemical reaction
    void AddChemicalReaction(FEChemicalReaction* pcr);

	// solute interface
public:
    typedef FESolutesMaterial::Point SolutesMaterial_t;

	int Solutes() override { return (int)m_pSolute.size(); }
	FESolute* GetSolute(int i) override { return m_pSolute[i]; }
    double GetEffectiveSoluteConcentration(FEMaterialPoint& mp, int soluteIndex) override {
        SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_c[soluteIndex];
    };
    double GetActualSoluteConcentration(FEMaterialPoint& mp, int soluteIndex) override {
        SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_ca[soluteIndex];
    };
    double GetPartitionCoefficient(FEMaterialPoint& mp, int soluteIndex) override {
        SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_k[soluteIndex];
    };
    vec3d GetSoluteFlux(FEMaterialPoint& mp, int soluteIndex) override {
        SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_j[soluteIndex];
    };
    double GetOsmolarity(const FEMaterialPoint& mp) override {
        const SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->Osmolarity();
    }
    double dkdc(const FEMaterialPoint& mp, int i, int j) override {
        const SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_dkdc[i][j];
    }
    double dkdJ(const FEMaterialPoint& mp, int soluteIndex) override {
        const SolutesMaterial_t* spt = (mp.ExtractData<SolutesMaterial_t>());
        return spt->m_dkdJ[soluteIndex];
    }
    FEOsmoticCoefficient* GetOsmoticCoefficient() override { return m_pOsmC; }

public:
    FEChemicalReaction*            GetReaction            (int i) { return m_pReact[i];  }
    
    int Reactions         () { return (int) m_pReact.size();    }
public:
    double    m_Rgas;            //!< universal gas constant
    double    m_Tabs;            //!< absolute temperature
    double    m_Fc;              //!< Faraday's constant

private: // material properties
	std::vector<FESolute*>  m_pSolute;      //!< pointer to solute materials
    FEOsmoticCoefficient*        m_pOsmC;        //!< pointer to osmotic coefficient material
    std::vector<FEChemicalReaction*>    m_pReact;        //!< pointer to chemical reactions

	DECLARE_FECORE_CLASS();
};
