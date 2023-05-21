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
#include <FEBioMech/FEElasticMaterial.h>
#include "FEFluid.h"
#include <FEBioMix/FESolute.h>
#include <FEBioMix/FESoluteInterface.h>
#include <FEBioMix/FEOsmoticCoefficient.h>
#include <FEBioMix/FEChemicalReaction.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEFluidSolutesMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEFluidSolutesMaterialPoint(FEMaterialPointData* pt);
    
    //! create a shallow copy
	FEMaterialPointData* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();

public:
    double Osmolarity() const;
    
public:
    // solutes material data
    int                 m_nsol;     //!< number of solutes
    vector<double>      m_c;        //!< effective solute concentration
    vector<double>      m_ca;       //!< actual solute concentration
    vector<vec3d>       m_gradc;    //!< spatial gradient of solute concentration
    vector<vec3d>       m_j;        //!< solute molar flux
    vector<double>      m_cdot;     //!< material time derivative of solute concentration following fluid
    double            m_psi;        //!< electric potential
    vec3d            m_Ie;          //!< current density
    double           m_pe;          //!< effective fluid pressure
    vector<double>    m_k;          //!< solute partition coefficient
    vector<double>    m_dkdJ;       //!< 1st deriv of m_k with strain (J)
    vector<double>    m_dkdJJ;      //!< 2nd deriv of m_k with strain (J)
    vector< vector<double> >    m_dkdc;            //!< 1st deriv of m_k with effective concentration
    vector< vector<double> >    m_dkdJc;        //!< cross deriv of m_k with J and c
    vector< vector< vector<double> > > m_dkdcc;    // 2nd deriv of m_k with c
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEFluidSolutes : public FEMaterial, public FESoluteInterface_T<FEFluidSolutesMaterialPoint>
{
public:
    FEFluidSolutes(FEModel* pfem);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
    
    //! performs initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
public:
    FEFluid* Fluid() { return m_pFluid; }
    
    //! calculate solute molar flux
    vec3d SoluteFlux(const FEMaterialPoint& pt, const int sol);
    
    //! calculate diffusive solute molar flux
    vec3d SoluteDiffusiveFlux(const FEMaterialPoint& pt, const int sol);
    
    //! actual concentration (as opposed to effective concentration)
    double ConcentrationActual(const FEMaterialPoint& pt, const int sol);
    
    //! actual fluid pressure (as opposed to effective pressure)
    double PressureActual(const FEMaterialPoint& pt);
    
    //! partition coefficient
    double PartitionCoefficient(const FEMaterialPoint& pt, const int sol);
    
    //! partition coefficients and their derivatives
    void PartitionCoefficientFunctions(const FEMaterialPoint& mp, vector<double>& kappa,
                                       vector<double>& dkdJ,
                                       vector< vector<double> >& dkdc);
    
    //! electric potential
    double ElectricPotential(const FEMaterialPoint& pt, const bool eform=false);
    
    //! current density
    vec3d CurrentDensity(const FEMaterialPoint& pt);
    
    //! solute density
    double SoluteDensity(const int sol) { return m_pSolute[sol]->Density(); }
    
    //! solute molar mass
    double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }
    
    //! solute charge number
    int SoluteChargeNumber(const int sol) { return m_pSolute[sol]->ChargeNumber(); }
    
    //! Add a chemical reaction
    void AddChemicalReaction(FEChemicalReaction* pcr);
    
    // solute interface
public:

    typedef FEFluidSolutesMaterialPoint SoluteMaterialPoint_t;

    int Solutes() override { return (int)m_pSolute.size(); }
    FESolute* GetSolute(int i) override { return m_pSolute[i]; }
    FEOsmoticCoefficient* GetOsmoticCoefficient() override { return m_pOsmC;  }
    double GetEffectiveSoluteConcentration(FEMaterialPoint& mp, int soluteIndex) override;
    double GetActualSoluteConcentration(FEMaterialPoint& mp, int soluteIndex) override { return ConcentrationActual(mp, soluteIndex); }
    double GetFreeDiffusivity(FEMaterialPoint& mp, int soluteIndex) override;
    double GetPartitionCoefficient(FEMaterialPoint& mp, int soluteIndex) override;
    vec3d GetSoluteFlux(FEMaterialPoint& mp, int soluteIndex) override { return SoluteFlux(mp, soluteIndex); }
    double GetOsmolarity(const FEMaterialPoint& mp) override;
    double GetElectricPotential(const FEMaterialPoint& mp) override { return ElectricPotential(mp); }
    vec3d GetCurrentDensity(const FEMaterialPoint& mp) override { return CurrentDensity(mp); }
    double dkdc(const FEMaterialPoint& mp, int i, int j) override;

public:
    FEChemicalReaction*            GetReaction            (int i) { return m_pReact[i];  }
    
    int Reactions         () { return (int) m_pReact.size();    }
    
public:
    double    m_Rgas;            //!< universal gas constant
    double    m_Tabs;            //!< absolute temperature
    double    m_Fc;              //!< Faraday's constant
    bool      m_diffMtmSupp;     //!< Toggle on or off diffusive mtm supply for fluid
    int        m_zmin;            //!< minimum charge number in mixture
    int        m_ndeg;            //!< polynomial degree of zeta in electroneutrality
    double              m_penalty;  //!< penalty for enforcing electroneutrality
    
private: // material properties
    FEFluid*                m_pFluid;       //!< pointer to fluid material
    std::vector<FESolute*>  m_pSolute;      //!< pointer to solute materials
    FEOsmoticCoefficient*        m_pOsmC;        //!< pointer to osmotic coefficient material
    std::vector<FEChemicalReaction*>    m_pReact;        //!< pointer to chemical reactions
    

    DECLARE_FECORE_CLASS();
};
