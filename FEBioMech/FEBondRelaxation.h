//
//  FEBondRelaxation.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEBondRelaxation__
#define __FEBioMech__FEBondRelaxation__

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for bond relaxation of reactive viscoelastic materials.
//! These materials need to define a relaxation function.
//!
class FEBondRelaxation : public FEMaterial
{
public:
	FEBondRelaxation(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEBondRelaxation() {}
    
	//! relaxation
	virtual double Relaxation(FEMaterialPoint& pt, const double t) = 0;
    
	void Init();
};

//-----------------------------------------------------------------------------
// This class implements constant relaxation triggered by any strain

class FEBondRelaxationConstIso : public FEBondRelaxation
{
public:
	enum { MAX_TERMS = 6 };
    
public:
    //! constructor
    FEBondRelaxationConstIso(FEModel* pfem);
    
    //! relaxation
    double Relaxation(FEMaterialPoint& pt, const double t);
    
    //! data initialization and checking
    void Init();
    
public:
	double	m_t[MAX_TERMS];	//!< relaxation times
    int     m_nt;           //!< number of non-zero relaxation times
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
#endif /* defined(__FEBioMech__FEBondRelaxation__) */
