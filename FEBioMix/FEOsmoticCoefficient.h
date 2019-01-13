#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for osmotic coefficient.
//! These materials need to define the osmotic coefficient and tangent functions.
//!
class FECORE_API FEOsmoticCoefficient : public FEMaterial
{
public:
	//! constructor
	FEOsmoticCoefficient(FEModel* pfem) : FEMaterial(pfem) {}
    
	//! osmotic coefficient
	virtual double OsmoticCoefficient(FEMaterialPoint& pt) = 0;
	
	//! tangent of osmotic coefficient with respect to strain
	virtual double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of osmotic coefficient with respect to concentration
	virtual double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp, const int isol) = 0;
};

