//
//  FEOsmCoefManning.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 1/8/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEOsmCoefManning_hpp
#define FEOsmCoefManning_hpp

#include "FEMultiphasic.h"
#include "FECore/LoadCurve.h"

//-----------------------------------------------------------------------------
// This class implements a material that has an osmotic coefficient behaving
// according to the Wells-Manning theory.  The Wells correction is provided
// by a loadcurve.

class FEOsmCoefManning : public FEOsmoticCoefficient
{
public:
    //! constructor
    FEOsmCoefManning(FEModel* pfem);
    
    //! destructor
    ~FEOsmCoefManning() {}
    
    //! Initialization
    bool Init();
    
    //! osmotic coefficient
    double OsmoticCoefficient(FEMaterialPoint& pt);
    
    //! Tangent of osmotic coefficient with respect to strain (J=detF)
    double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp);
    
    //! Tangent of osmotic coefficient with respect to concentration
    double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp, const int isol);

    //! Manning response
    double OsmoticCoefficient_Manning(FEMaterialPoint& pt);
    double Tangent_OsmoticCoefficient_Strain_Manning(FEMaterialPoint& mp);
    double Tangent_OsmoticCoefficient_Concentration_Manning(FEMaterialPoint& mp, const int isol);

    //! Wells response
    double OsmoticCoefficient_Wells(FEMaterialPoint& pt);
    double Tangent_OsmoticCoefficient_Strain_Wells(FEMaterialPoint& mp);
    double Tangent_OsmoticCoefficient_Concentration_Wells(FEMaterialPoint& mp, const int isol);
    
    //! bookkeeping
    void Serialize(DumpStream& ar);
    bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval);
    
public:
    double  m_ksi;              //!< Manning parameter
    int		m_sol;              //!< global id of co-ion
    int		m_lsol;             //!< local id of co-ion
    double  m_osmc;             //!< osmotic coefficient for Wells correction (mobile ion - mobile interaction)
    int     m_lc;               //!< load curve number for Wells correction
    FEMultiphasic*  m_pMP;      //!< pointer to ancestor multiphasic material
    FELoadCurve*    m_pLC;      //!< pointer to load curve
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* FEOsmCoefManning_hpp */
