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
#include "FEReaction.h"
#include "febiomix_api.h"
#include "FESolute.h"

//-----------------------------------------------------------------------------
//! Base class for membrane reaction rates.

class FEBIOMIX_API FEMembraneReactionRate : public FEMaterialProperty
{
public:
    //! constructor
    FEMembraneReactionRate(FEModel* pfem) : FEMaterialProperty(pfem), m_pReact(nullptr) {}
    
    //! reaction rate at material point
    virtual double ReactionRate(FEMaterialPoint& pt) = 0;
    
    //! tangent of reaction rate with area strain at material point
    virtual double Tangent_ReactionRate_Strain(FEMaterialPoint& pt) = 0;
    
    //! tangent of reaction rate with effective fluid pressure at material point
    virtual double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) = 0;
    virtual double Tangent_ReactionRate_Pe(FEMaterialPoint& pt) = 0;
    virtual double Tangent_ReactionRate_Pi(FEMaterialPoint& pt) = 0;

    //! tangent of reaction rate with effective solute concentration at material point
    virtual double Tangent_ReactionRate_Concentration(FEMaterialPoint& pt, const int isol) = 0;
    virtual double Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol) = 0;
    virtual double Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol) = 0;
    
    //! reset, initialize and update chemical reaction data in the FESolutesMaterialPoint
    virtual void ResetElementData(FEMaterialPoint& mp) {}
    virtual void InitializeElementData(FEMaterialPoint& mp) {}
    virtual void UpdateElementData(FEMaterialPoint& mp) {}
    
public:
    FEReaction*    m_pReact;    //!< pointer to parent reaction

    FECORE_BASE_CLASS(FEMembraneReactionRate)
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FEInternalReactantSpeciesRef : public FEReactionSpeciesRef
{
public: FEInternalReactantSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEInternalReactantSpeciesRef)
};

class FEBIOMIX_API FEInternalProductSpeciesRef : public FEReactionSpeciesRef {
public: FEInternalProductSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEInternalProductSpeciesRef)
};

class FEBIOMIX_API FEExternalReactantSpeciesRef : public FEReactionSpeciesRef
{
public: FEExternalReactantSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEExternalReactantSpeciesRef)
};

class FEBIOMIX_API FEExternalProductSpeciesRef : public FEReactionSpeciesRef {
public: FEExternalProductSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEExternalProductSpeciesRef)
};

//-----------------------------------------------------------------------------
//! Base class for membrane reactions.
class FEBIOMIX_API FEMembraneReaction : public FEReaction
{
public:
    //! constructor
    FEMembraneReaction(FEModel* pfem);
    
    //! get solute (use only during initialization)
    FESoluteData* GetSolute(int nsol);
    
    //! initialization
    bool Init() override;
    
public:
    //! set the forward reaction rate
    void SetForwardReactionRate(FEMembraneReactionRate* pfwd) { m_pFwd = pfwd; }
    
    //! set the reverse reaction rate
    void SetReverseReactionRate(FEMembraneReactionRate* prev) { m_pRev = prev; }
    
public:
    //! reset, initialize and update optional chemical reaction data in the FESolutesMaterialPoint
    void ResetElementData(FEMaterialPoint& mp)
    {
        if (m_pFwd) m_pFwd->ResetElementData(mp);
        if (m_pRev) m_pRev->ResetElementData(mp);
    }
    void InitializeElementData(FEMaterialPoint& mp)
    {
        if (m_pFwd) m_pFwd->InitializeElementData(mp);
        if (m_pRev) m_pRev->InitializeElementData(mp);
    }
    void UpdateElementData(FEMaterialPoint& mp)
    {
        if (m_pFwd) m_pFwd->UpdateElementData(mp);
        if (m_pRev) m_pRev->UpdateElementData(mp);
    }
    
public:
    //! molar supply at material point
    virtual double ReactionSupply(FEMaterialPoint& pt) = 0;
    
    //! tangent of molar supply with strain at material point
    virtual double Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) = 0;
    
    //! tangent of molar supply with effective pressure at material point
    virtual double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) = 0;
    virtual double Tangent_ReactionSupply_Pe(FEMaterialPoint& pt) = 0;
    virtual double Tangent_ReactionSupply_Pi(FEMaterialPoint& pt) = 0;

    //! tangent of molar supply with effective concentration at material point
    virtual double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) = 0;
    virtual double Tangent_ReactionSupply_Ce(FEMaterialPoint& pt, const int sol) = 0;
    virtual double Tangent_ReactionSupply_Ci(FEMaterialPoint& pt, const int sol) = 0;

public:
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
public:
    FEMembraneReactionRate*    m_pFwd;        //!< pointer to forward reaction rate
    FEMembraneReactionRate*    m_pRev;        //!< pointer to reverse reaction rate
    
    vector<FEReactantSpeciesRef*> m_vRtmp;	//!< helper variable for reading in stoichiometric coefficients for reactants
    vector<FEProductSpeciesRef*> m_vPtmp;	//!< helper variable for reading in stoichiometric coefficients for products
    vector<FEInternalReactantSpeciesRef*> m_vRitmp;	//!< helper variable for reading in stoichiometric coefficients for internal reactants
    vector<FEInternalProductSpeciesRef*> m_vPitmp;	//!< helper variable for reading in stoichiometric coefficients for internal products
    vector<FEExternalReactantSpeciesRef*> m_vRetmp;	//!< helper variable for reading in stoichiometric coefficients for external reactants
    vector<FEExternalProductSpeciesRef*> m_vPetmp;	//!< helper variable for reading in stoichiometric coefficients for external products

    intmap          m_solR;         //!< stoichiometric coefficients of solute reactants
    intmap          m_solP;         //!< stoichiometric coefficients of solute products
    intmap          m_sbmR;         //!< stoichiometric coefficients of solid-bound reactants
    intmap          m_sbmP;         //!< stoichiometric coefficients of solid-bound products
    intmap          m_solRi;        //!< stoichiometric coefficients of internal solute reactants
    intmap          m_solPi;        //!< stoichiometric coefficients of internal solute products
    intmap          m_solRe;        //!< stoichiometric coefficients of external solute reactants
    intmap          m_solPe;        //!< stoichiometric coefficients of external solute products

public:
    int             m_nsol;         //!< number of solutes in the mixture
    double          m_Vbar;         //!< weighted molar volume of reactants and products
    bool            m_Vovr;         //!< override flag for m_Vbar
    vector<int>     m_vR;           //!< stoichiometric coefficients of reactants
    vector<int>     m_vP;           //!< stoichiometric coefficients of products
    vector<int>     m_v;            //!< net stoichiometric coefficients of reactants and products
    int             m_NSOL;         //!< number of solutes in the model
    vector<int>     m_z;            //!< charge number of all solutes
    vector<int>     m_vRi;          //!< stoichiometric coefficients of reactants
    vector<int>     m_vPi;          //!< stoichiometric coefficients of products
    vector<int>     m_vi;           //!< net stoichiometric coefficients of reactants and products
    vector<int>     m_vRe;          //!< stoichiometric coefficients of reactants
    vector<int>     m_vPe;          //!< stoichiometric coefficients of products
    vector<int>     m_ve;           //!< net stoichiometric coefficients of reactants and products

    DECLARE_FECORE_CLASS();
    FECORE_BASE_CLASS(FEMembraneReaction)
};
