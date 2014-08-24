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
//! Material point data for Nims reaction rate material
class FENimsMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FENimsMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
	//! copy material point data
	FEMaterialPoint* Copy();
    
	//! Initialize material point data
	void Init(bool bflag);
    
	//! Serialize data to archive
	void Serialize(DumpFile& ar);
    
	//! data streaming
	void ShallowCopy(DumpStream& dmp, bool bsave);
    
public:
    int     m_lid;  //!< local id of solute
	double	m_cmax;	//!< maximum solute concentration over past history
};


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
	FEReactionRateNims(FEModel* pfem) : FEReactionRate(pfem) { m_trel = 0; m_lid = -1; }
	
	//! data initialization and checking
	void Init();
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt);
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
    
public:
    int     m_sol;                  //!< solute id (1-based)
    int     m_lid;                  //!< local id of solute (zero-based)
    double  m_k0;                   //!< reaction rate at zero concentration
    double  m_cc;                   //!< concentration cc
    double  m_kc;                   //!< reaction rate at cc;
    double  m_cr;                   //!< concentration cr;
    double  m_kr;                   //!< reaction rate at cr;
    double  m_trel;                 //!< release time;
	
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMix__FEReactionRateNims__) */
