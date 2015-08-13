//
//  FEDamageCDF.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEDamageCDF__
#define __FEBioMech__FEDamageCDF__

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Virtual base class for damage cumulative distribution functions

class FEDamageCDF : public FEMaterial
{
public:
	FEDamageCDF(FEModel* pfem) : FEMaterial(pfem) {}
    
	//! damage
	virtual double Damage(FEMaterialPoint& pt) = 0;
};

//-----------------------------------------------------------------------------
// Simo damage cumulative distribution function
// Simo, CMAME 60 (1987), 153-173

class FEDamageCDFSimo : public FEDamageCDF
{
public:
	FEDamageCDFSimo(FEModel* pfem);
	~FEDamageCDFSimo() {}
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
public:
	double	m_alpha;			//!< parameter alpha
	double	m_beta;             //!< parameter beta
    
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Log-normal damage cumulative distribution function

class FEDamageCDFLogNormal : public FEDamageCDF
{
public:
	FEDamageCDFLogNormal(FEModel* pfem);
	~FEDamageCDFLogNormal() {}
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
public:
	double	m_mu;               //!< mean on log scale
	double	m_sigma;            //!< standard deviation on log scale
    double  m_Dmax;              //!< maximum allowable damage
    
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Weibull damage cumulative distribution function

class FEDamageCDFWeibull : public FEDamageCDF
{
public:
	FEDamageCDFWeibull(FEModel* pfem);
	~FEDamageCDFWeibull() {}
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
public:
	double	m_alpha;            //!< exponent alpha
	double	m_mu;               //!< mean mu
    double  m_Dmax;              //!< maximum allowable damage
    
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Step cumulative distribution function (sudden fracture)

class FEDamageCDFStep : public FEDamageCDF
{
public:
	FEDamageCDFStep(FEModel* pfem);
	~FEDamageCDFStep() {}
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
public:
	double	m_mu;               //!< threshold mu
    double  m_Dmax;              //!< maximum allowable damage
    
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Piecewise S-shaped quintic polynomial cumulative distribution function

class FEDamageCDFPQP : public FEDamageCDF
{
public:
	FEDamageCDFPQP(FEModel* pfem);
	~FEDamageCDFPQP() {}
    
	//! damage
	double Damage(FEMaterialPoint& pt);
    
	void Init();
    
public:
	double	m_mumin;            //!< mu threshold
	double	m_mumax;            //!< mu cap
    double  m_Dmax;              //!< maximum allowable damage
    
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEDamageCDF__) */
