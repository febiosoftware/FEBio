//
//  FESolubManning.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 1/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FESolubManning_hpp
#define FESolubManning_hpp
#include <FECore/FEFunction1D.h>
#include "FESolute.h"

class FEMultiphasic;

//-----------------------------------------------------------------------------
// This class implements a material that has a solute solubility that follows
// the Wells-Manning theory.  The mobile ion-mobile ion correction from Wells
// is provided by a loadcurve.

//-----------------------------------------------------------------------------
class FESolubManning : public FESoluteSolubility
{
public:
    //! constructor
    FESolubManning(FEModel* pfem);
    
    //! Initialization
    bool Init() override;
    
    //! solubility
    double Solubility(FEMaterialPoint& pt) override;
    
    //! Tangent of solubility with respect to strain
    double Tangent_Solubility_Strain(FEMaterialPoint& mp) override;
    
    //! Tangent of solubility with respect to concentration
    double Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol) override;
    
    //! Cross derivative of solubility with respect to strain and concentration
    double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp, const int isol) override;
    
    //! Second derivative of solubility with respect to strain
    double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! Second derivative of solubility with respect to concentration
    double Tangent_Solubility_Concentration_Concentration(FEMaterialPoint& mp, const int isol, const int jsol) override;
    
    //! Manning response
    double Solubility_Manning(FEMaterialPoint& mp);
    double Tangent_Solubility_Strain_Manning(FEMaterialPoint& mp);
    double Tangent_Solubility_Concentration_Manning (FEMaterialPoint& mp, const int isol);
    
    //! Wells response
    double Solubility_Wells(FEMaterialPoint& mp);
    double Tangent_Solubility_Strain_Wells(FEMaterialPoint& mp);
    double Tangent_Solubility_Concentration_Wells(FEMaterialPoint& mp, const int isol);

public:
    double  m_ksi;              //!< Manning parameter
    int		m_sol;              //!< global id of co-ion
    int		m_lsol;             //!< local id of co-ion
    bool    m_bcoi;             //!< true if this solute is the co-ion
    FEFunction1D	m_solub;    //!< solubility from Wells correction

    FEMultiphasic*  m_pMP;      //!< pointer to ancestor multiphasic material
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* FESolubManning_hpp */
