#pragma once
#include "FESolidMaterial.h"
#include "FEELasticMaterialPoint.h"

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEElasticMaterial : public FESolidMaterial
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
	FEMatAxis	m_Q;

	DECLARE_FECORE_CLASS();
};
