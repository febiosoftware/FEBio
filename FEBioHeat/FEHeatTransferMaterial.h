#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Material point class for heat transfer materials.
class FEHeatMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FEHeatMaterialPoint* pt = new FEHeatMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void ShallowCopy(DumpStream& dmp, bool bsave)
	{
		if (bsave) { dmp << m_q; } else { dmp >> m_q; }
		if (m_pt) m_pt->ShallowCopy(dmp, bsave);
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (m_pt) m_pt->Init(bflag);
	}

public:
	vec3d	m_q;	//!< heat flux
};

//-----------------------------------------------------------------------------
// Base class for heat-transfer problems
class FEHeatTransferMaterial : public FEMaterial
{
public:
	//! constructor
	FEHeatTransferMaterial(FEModel* pfem) : FEMaterial(pfem) {}

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() { return new FEHeatMaterialPoint; }

	//! get the material's conductivity
	virtual void Conductivity(double D[3][3]) = 0;

	//! get the material's capacitance
	virtual double Capacitance() = 0;

	//! get the material's density
	virtual double Density() = 0;

	//! get the heat flux
	virtual vec3d HeatFlux(vec3d gradT) = 0;
};
