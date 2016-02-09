#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
class FEMembraneMaterialPoint : public FEMaterialPoint
{
public:
	FEMembraneMaterialPoint()
	{
		s[0] = s[1] = s[2] = 0;
	}

	FEMaterialPoint* Copy()
	{
		return new FEMembraneMaterialPoint();
	}

	void Serialize(DumpStream& ar)
	{
		FEMaterialPoint::Serialize(ar);
	}

	void Init(bool bflag)
	{
		FEMaterialPoint::Init(bflag);
	}

public:
	// calculate membrane strain
	void strain(double e[3]);

public:
	double	g[6];	// deformation gradient
	double	s[3];	// in-plane PK2 stress
};

//-----------------------------------------------------------------------------
class FEMembraneMaterial : public FEMaterial
{
public:
	FEMembraneMaterial(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEMembraneMaterial() {}

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() { return new FEMembraneMaterialPoint; }

public:
	//! calculate in-plane membrane stress
	virtual void Stress(FEMaterialPoint& mp, double s[3]) = 0;

	//! calculate in-plane membrane tangent
	virtual void Tangent(FEMaterialPoint& mp, double D[3][3]) = 0;
};

//-----------------------------------------------------------------------------
// class for triangular membranes
class FEElasticMembrane : public FEMembraneMaterial
{
public:
	//! constructor
	FEElasticMembrane(FEModel* pfem) : FEMembraneMaterial(pfem) {}

public:
	//! calculate in-plane membrane stress
	void Stress(FEMaterialPoint& mp, double s[3]);

	//! calculate in-plane membrane tangent
	void Tangent(FEMaterialPoint& mp, double D[3][3]);

public:
	double	m_E;
	double	m_v;

	DECLARE_PARAMETER_LIST();
};
