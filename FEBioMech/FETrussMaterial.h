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
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (m_pt) m_pt->Init(bflag);
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

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
