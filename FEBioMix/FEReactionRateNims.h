//
//  FEReactionRateNims.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEReactionRateNims__
#define __FEBioMix__FEReactionRateNims__

#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Concentration-history-dependent reaction rate.
//! Reaction rate depends on concentration of a solute (e.g., growth factor)
//! and whether solute has been released (removed) at some release time.
//! Before release, reaction rate varies linearly with current solute
//! concentration c, from k0 at c=0 to kc at c=cc, then holds constant at kc.
//! After release, reaction rate depends on past history of exposure to solute.
//! Increases from k0 at cmax=0 to kr at cmax=cr, then holds constant at kr.
//! Release time is trel.

class FEReactionRateNims : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateNims(FEModel* pfem) : FEReactionRate(pfem)
    { m_trel = 0; m_lid = m_cmax = -1; }
	
	//! data initialization and checking
	void Init();
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt);
	
    //! reset, initialize and update chemical reaction data in the FESolutesMaterialPoint
    void ResetElementData(FEMaterialPoint& mp);
    void InitializeElementData(FEMaterialPoint& mp);
    void UpdateElementData(FEMaterialPoint& mp);
    
public:
    int     m_sol;                  //!< solute id (1-based)
    int     m_lid;                  //!< local id of solute (zero-based)
    double  m_k0;                   //!< reaction rate at zero concentration
    double  m_cc;                   //!< concentration cc
    double  m_kc;                   //!< reaction rate at cc;
    double  m_cr;                   //!< concentration cr;
    double  m_kr;                   //!< reaction rate at cr;
    double  m_trel;                 //!< release time;
    int     m_cmax;                 //!< index of entry in m_crd
	
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMix__FEReactionRateNims__) */
