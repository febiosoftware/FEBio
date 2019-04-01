#pragma once
#include "FECore/FEMaterial.h"
#include "FEBioMix/FESolutesMaterialPoint.h"
#include <FEBioMech/FESSIShellDomain.h>
#include "FEReaction.h"
#include "FESolute.h"

//-----------------------------------------------------------------------------
//! Base class for membrane reaction rates.

class FEBIOMIX_API FEMembraneReactionRate : public FEMaterial
{
public:
    //! constructor
    FEMembraneReactionRate(FEModel* pfem) : FEMaterial(pfem) {}
    
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
};

//-----------------------------------------------------------------------------
//! Base class for membrane reactions.
class FEBIOMIX_API FEMembraneReaction : public FEReaction
{
public:
    //! constructor
    FEMembraneReaction(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
public:
    void SetParameter(FEParam& p) override;
    
    bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) override;
    
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
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
    //! get solute data
    FESoluteData* FindSoluteData(int nid);
    
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
    FEMembraneReactionRate*    m_pFwd;        //!< pointer to forward reaction rate
    FEMembraneReactionRate*    m_pRev;        //!< pointer to reverse reaction rate
    
public:
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
    int             m_vRtmp;        //!< helper variable for reading in stoichiometric coefficients for reactants
    int             m_vPtmp;        //!< helper variable for reading in stoichiometric coefficients for products
    int             m_NSOL;         //!< number of solutes in the model
    vector<int>     m_vRi;          //!< stoichiometric coefficients of reactants
    vector<int>     m_vPi;          //!< stoichiometric coefficients of products
    vector<int>     m_vi;           //!< net stoichiometric coefficients of reactants and products
    int             m_vRitmp;       //!< helper variable for reading in stoichiometric coefficients for reactants
    int             m_vPitmp;       //!< helper variable for reading in stoichiometric coefficients for products
    vector<int>     m_vRe;          //!< stoichiometric coefficients of reactants
    vector<int>     m_vPe;          //!< stoichiometric coefficients of products
    vector<int>     m_ve;           //!< net stoichiometric coefficients of reactants and products
    int             m_vRetmp;       //!< helper variable for reading in stoichiometric coefficients for reactants
    int             m_vPetmp;       //!< helper variable for reading in stoichiometric coefficients for products

    DECLARE_FECORE_CLASS();
};
