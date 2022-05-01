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

//-----------------------------------------------------------------------------
//! Base class for reaction rates.

class FEBIOMIX_API FEReactionRate : public FEMaterialProperty
{
public:
    //! constructor
    FEReactionRate(FEModel* pfem) : FEMaterialProperty(pfem), m_pReact(nullptr) {}
    
    //! reaction rate at material point
    virtual double ReactionRate(FEMaterialPoint& pt) = 0;
    
    //! tangent of reaction rate with strain at material point
    virtual mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) = 0;
    
    //! tangent of reaction rate with effective fluid pressure at material point
    virtual double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) = 0;
    
    //! reset, initialize and update chemical reaction data in the FESolutesMaterialPoint
    virtual void ResetElementData(FEMaterialPoint& mp) {}
    virtual void InitializeElementData(FEMaterialPoint& mp) {}
    virtual void UpdateElementData(FEMaterialPoint& mp) {}
    
public:
    FEReaction*    m_pReact;    //!< pointer to parent reaction

    FECORE_BASE_CLASS(FEReactionRate)
};

//-----------------------------------------------------------------------------
//! Base class for chemical reactions.

class FEBIOMIX_API FEChemicalReaction : public FEReaction
{
public:
	//! constructor
    FEChemicalReaction(FEModel* pfem);
    
	//! initialization
	bool Init() override;

public:
    //! set the forward reaction rate
    void SetForwardReactionRate(FEReactionRate* pfwd) { m_pFwd = pfwd; }
    
    //! set the reverse reaction rate
    void SetReverseReactionRate(FEReactionRate* prev) { m_pRev = prev; }
    
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
    virtual mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) = 0;
    
    //! tangent of molar supply with effective pressure at material point
    virtual double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) = 0;
    
    //! tangent of molar supply with effective concentration at material point
    virtual double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) = 0;
    
public:
	//! Serialization
	void Serialize(DumpStream& ar) override;

public:
    vector<FEReactantSpeciesRef*> m_vRtmp;	//!< helper variable for reading in stoichiometric coefficients for reactants
    vector<FEProductSpeciesRef*> m_vPtmp;	//!< helper variable for reading in stoichiometric coefficients for products

    FEReactionRate*    m_pFwd;        //!< pointer to forward reaction rate
    FEReactionRate*    m_pRev;        //!< pointer to reverse reaction rate
    
public:
	intmap			m_solR;		//!< stoichiometric coefficients of solute reactants (input)
	intmap			m_solP;		//!< stoichiometric coefficients of solute products (input)
	intmap			m_sbmR;		//!< stoichiometric coefficients of solid-bound reactants (input)
	intmap			m_sbmP;		//!< stoichiometric coefficients of solid-bound products (input)
	
public:
	int				m_nsol;		//!< number of solutes in the mixture
	vector<int>		m_vR;		//!< stoichiometric coefficients of reactants
	vector<int>		m_vP;		//!< stoichiometric coefficients of products
	vector<int>		m_v;		//!< net stoichiometric coefficients of reactants and products
    double          m_Vbar;     //!< weighted molar volume of reactants and products
    bool            m_Vovr;     //!< override flag for m_Vbar

	DECLARE_FECORE_CLASS();
    FECORE_BASE_CLASS(FEChemicalReaction)
};
