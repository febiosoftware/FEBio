/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEFluidFSI.h"
#include "FEBiphasicFSI.h"
#include "FEMultiphasicFSI.h"
#include <FEBioMix/FEHydraulicPermeability.h>
#include <FEBioMech/FEBodyForce.h>
#include "FEFluid.h"
#include <FEBioMix/FESolute.h>
#include <FEBioMix/FESoluteInterface.h>
#include <FEBioMix/FEOsmoticCoefficient.h>
#include <FEBioMix/FEChemicalReaction.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! FSI material point class.
//
class FEBIOFLUID_API FEMultiphasicFSIMaterialPoint : public FEMaterialPointData
{
public:
    //! constructor
    FEMultiphasicFSIMaterialPoint(FEMaterialPointData* pt);
    
    //! create a shallow copy
	FEMaterialPointData* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();

public:
    double Osmolarity() const;
    
public:
    // Multiphasic FSI material data
    int                 m_nsol;     //!< number of solutes
    vector<double>      m_c;        //!< effective solute concentration
    vector<double>      m_ca;        //!< effective solute concentration
    vector<vec3d>       m_gradc;    //!< spatial gradient of solute concentration
    vector<vec3d>       m_j;        //!< solute molar flux
    vector<double>      m_cdot;     //!< material time derivative of solute concentration following fluid
    double            m_psi;        //!< electric potential
    vec3d            m_Ie;        //!< current density
    double            m_pe;        //!< effective fluid pressure
    double            m_cF;        //!< fixed charge density in current configuration
    vector<double>    m_k;        //!< solute partition coefficient
    vector<double>    m_dkdJ;        //!< 1st deriv of m_k with strain (J)
    vector<double>    m_dkdJJ;    //!< 2nd deriv of m_k with strain (J)
    vector< vector<double> >    m_dkdc;            //!< 1st deriv of m_k with effective concentration
    vector< vector<double> >    m_dkdJc;        //!< cross deriv of m_k with J and c
    vector< vector< vector<double> > > m_dkdcc;    // 2nd deriv of m_k with c
};

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEMultiphasicFSI : public FEBiphasicFSI, public FESoluteInterface_T<FEMultiphasicFSIMaterialPoint>
{
public:
    FEMultiphasicFSI(FEModel* pfem);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
    
    //! performs initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
public:
    //! calculate inner stress at material point
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! return the diffusivity tensor as a matrix
    void Diffusivity(double k[3][3], FEMaterialPoint& pt, int sol);
    
    //! return the diffusivity as a tensor
    mat3ds Diffusivity(FEMaterialPoint& pt, int sol);
    
    //! return the inverse diffusivity as a tensor
    mat3ds InvDiffusivity(FEMaterialPoint& pt, int sol);
    
    //! return the tangent diffusivity tensor
    tens4dmm Diffusivity_Tangent_Strain(FEMaterialPoint& pt, int isol);
    
    //! return the tangent diffusivity tensor
    mat3ds Diffusivity_Tangent_Concentration(FEMaterialPoint& pt, int isol, int jsol);
    
    //! return the diffusivity property
    FESoluteDiffusivity* GetDiffusivity(int sol) { return m_pSolute[sol]->m_pDiff; }
    
    //! calculate solute molar flux
    vec3d SoluteFlux(FEMaterialPoint& pt, const int sol);
    
    //! actual concentration (as opposed to effective concentration)
    double ConcentrationActual(FEMaterialPoint& pt, const int sol);
    
    //! actual fluid pressure (as opposed to effective pressure)
    double PressureActual(FEMaterialPoint& pt);
    
    //! fixed charge density
    virtual double FixedChargeDensity(FEMaterialPoint& pt);
    
    //! partition coefficient
    double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
    
    //! partition coefficients and their derivatives
    void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                       vector<double>& dkdJ,
                                       vector< vector<double> >& dkdc);
    
    //! electric potential
    double ElectricPotential(FEMaterialPoint& pt, const bool eform=false);
    
    //! current density
    vec3d CurrentDensity(FEMaterialPoint& pt);
    
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
    int Solutes() override { return (int)m_pSolute.size(); }
    FESolute* GetSolute(int i) override { return m_pSolute[i]; }
    double GetReferentialFixedChargeDensity(const FEMaterialPoint& mp) override;
    FEOsmoticCoefficient* GetOsmoticCoefficient() override { return m_pOsmC;  }
    double GetFixedChargeDensity(const FEMaterialPoint& mp) override {
        const FEMultiphasicFSIMaterialPoint* spt = (mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        return spt->m_cF;
    }

public: // from FEBiphasicInterface
    double GetActualFluidPressure(const FEMaterialPoint& mp) override {
        const FEMultiphasicFSIMaterialPoint* pt = (mp.ExtractData<FEMultiphasicFSIMaterialPoint>());
        return pt->m_pe;
    }

public:
    FEChemicalReaction*            GetReaction            (int i) { return m_pReact[i];  }
    
    int Reactions         () { return (int) m_pReact.size();    }
    
public: // material parameters
    FEParamDouble       m_cFr;      //!< fixed charge density in reference configurations TODO: gradphisr
    vector<FEBodyForce*>    m_mbf;       //!< body forces acting on this multiphasic material solutes
    double    m_Rgas;            //!< universal gas constant
    double    m_Tabs;            //!< absolute temperature
    double    m_Fc;              //!< Faraday's constant
    bool      m_diffMtmSupp;     //!< Toggle on or off diffusive mtm supply for fluid
    int        m_zmin;            //!< minimum charge number in mixture
    int        m_ndeg;            //!< polynomial degree of zeta in electroneutrality
    double              m_penalty;  //!< penalty for enforcing electroneutrality
    
protected: // material properties
    std::vector<FESolute*>  m_pSolute;      //!< pointer to solute materials
    FEOsmoticCoefficient*        m_pOsmC;        //!< pointer to osmotic coefficient material
    std::vector<FEChemicalReaction*>    m_pReact;        //!< pointer to chemical reactions
    
    DECLARE_FECORE_CLASS();
};
