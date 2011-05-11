#pragma once

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
class FETrussMaterial : public FEMaterial
{
public:
	FETrussMaterial(){}
	~FETrussMaterial(){}

public:
	double	m_E;	// Elastic modulus

public:
	//! calculate Kirchhoff stress of truss
	virtual double Stress(FEMaterialPoint& pt);

	//! calculate elastic tangent
	virtual double Tangent(FEMaterialPoint& pt);

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() { return new FETrussMaterialPoint; }

	// declare as registered
	DECLARE_REGISTERED(FETrussMaterial);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
