#pragma once

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Material point class for truss materials
class FETrussMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FETrussMaterialPoint* pt = new FETrussMaterialPoint(*this);
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

	void Serialize(DumpStream& ar)
	{
		FEMaterialPoint::Serialize(ar);
		ar & m_l & m_tau;
	}

	void Init()
	{
		FEMaterialPoint::Init();
		m_l = 1;
		m_tau = 0;
	}

public:
	double	m_l;	// strech
	double	m_tau;	// Kirchoff stress
};

//-----------------------------------------------------------------------------
// Base class for truss element materials
class FETrussMaterial : public FEMaterial
{
public:
	FETrussMaterial(FEModel* pfem) : FEMaterial(pfem) {}
	~FETrussMaterial(){}

public:
	double	m_E;	// Elastic modulus

public:
	//! calculate Kirchhoff stress of truss
	virtual double Stress(FEMaterialPoint& pt);

	//! calculate elastic tangent
	virtual double Tangent(FEMaterialPoint& pt);

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() override { return new FETrussMaterialPoint; }

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
