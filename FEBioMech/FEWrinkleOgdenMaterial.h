#pragma once
#include "FEMembraneMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a wrinkle model for an Ogden material
//
class FEWrinkleOgdenMaterial : public FEMembraneMaterial
{
public:
	FEWrinkleOgdenMaterial(FEModel* pfem);

public: // material parameters
	double	m_u;
	double	m_a;
	double	m_l0[2];	// pre-stretch
	bool	m_bwrinkle;	// wrinkle flag

public: // material interface

	//! calculate stress at material point
	void Stress(FEMaterialPoint& pt, double s[3]);

	//! calculate tangent stiffness at material point
	void Tangent(FEMaterialPoint& pt, double D[3][3]);

protected:
	void principals(FEMaterialPoint& pt, double l[2], double v[4]);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
