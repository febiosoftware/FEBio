#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for solvent supply.
//! These materials need to define the supply and tangent supply functions.
//!
class FESolventSupply : public FEMaterial
{
public:
	FESolventSupply(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FESolventSupply(){}
	
	//! solvent supply
	virtual double Supply(FEMaterialPoint& pt) = 0;
	
	//! tangent of solvent supply with respect to strain
	virtual mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of solvent supply with respect to pressure
	virtual double Tangent_Supply_Pressure(FEMaterialPoint& mp) = 0;
	
	//! tangent of solvent supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp, const int isol);
	
	void Init();
};

