#pragma once
#include "FESolidMaterial.h"
#include "FEElasticMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEBIOMECH_API FEElasticMaterial : public FESolidMaterial
{
public:
	//! constructor 
	FEElasticMaterial(FEModel* pfem);

	//! destructor
	~FEElasticMaterial();

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override { return new FEElasticMaterialPoint; }

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);

protected:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEElasticStress : public FEDomainParameter
{
public:
	FEElasticStress();
	FEParamValue value(FEMaterialPoint& mp) override;
};
