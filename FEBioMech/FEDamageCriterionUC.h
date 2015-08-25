//
//  FEDamageCriterionUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 9/19/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEDamageCriterionUC__
#define __FEBioMech__FEDamageCriterionUC__

#include "FECore/FEMaterial.h"
#include "FEDamageMaterialUC.h"

//-----------------------------------------------------------------------------
// Virtual base class for damage criterion

class FEDamageCriterionUC : public FEMaterial
{
public:
	FEDamageCriterionUC(FEModel* pfem) : FEMaterial(pfem) {}
    
	//! damage
	virtual double DamageCriterion(FEMaterialPoint& pt) = 0;
    
};

//-----------------------------------------------------------------------------
// Simo's damage criterion

class FEDamageCriterionUCSimo : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCSimo(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// Strain energy density as damage criterion

class FEDamageCriterionUCSED : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCSED(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// von Mises stress as damage criterion

class FEDamageCriterionUCVMS : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCVMS(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max shear stress as damage criterion

class FEDamageCriterionUCMSS : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCMSS(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal stress as damage criterion

class FEDamageCriterionUCMNS : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCMNS(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max normal Lagrange strain as damage criterion

class FEDamageCriterionUCMNLS : public FEDamageCriterionUC
{
public:
	FEDamageCriterionUCMNLS(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
	//! damage
	double DamageCriterion(FEMaterialPoint& pt);
};

//-----------------------------------------------------------------------------
// max principal stretch ratio as damage criterion

class FEDamageCriterionUCMPSR : public FEDamageCriterionUC
{
public:
    FEDamageCriterionUCMPSR(FEModel* pfem) : FEDamageCriterionUC(pfem) {}
    
    //! damage
    double DamageCriterion(FEMaterialPoint& pt);
};

#endif /* defined(__FEBioMech__FEDamageCriterionUC__) */
