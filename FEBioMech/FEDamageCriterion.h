//
//  FEDamageCriterion.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/18/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEDamageCriterion__
#define __FEBioMech__FEDamageCriterion__

#include "FECore/FEMaterial.h"
#include "FEDamageMaterial.h"

//-----------------------------------------------------------------------------
// Virtual base class for damage criterion

class FEDamageCriterion : public FEMaterial
{
public:
	FEDamageCriterion(FEModel* pfem) : FEMaterial(pfem) {}
    
	//! damage
	virtual double DamageCriterion(FEMaterialPoint& pt) = 0;
    
};

//-----------------------------------------------------------------------------
// Simo's damage criterion

class FEDamageCriterionSimo : public FEDamageCriterion
{
public:
	FEDamageCriterionSimo(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// Strain energy density as damage criterion

class FEDamageCriterionSED : public FEDamageCriterion
{
public:
	FEDamageCriterionSED(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// von Mises stress as damage criterion

class FEDamageCriterionVMS : public FEDamageCriterion
{
public:
	FEDamageCriterionVMS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max shear stress as damage criterion

class FEDamageCriterionMSS : public FEDamageCriterion
{
public:
	FEDamageCriterionMSS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal stress as damage criterion

class FEDamageCriterionMNS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal Lagrange strain as damage criterion

class FEDamageCriterionMNLS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNLS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

#endif /* defined(__FEBioMech__FEDamageCriterion__) */
