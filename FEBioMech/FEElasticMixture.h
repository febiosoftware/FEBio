#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for mixtures
//!
class FEElasticMixtureMaterialPoint : public FEMaterialPoint
{
public:
	FEElasticMixtureMaterialPoint() { m_pt = new FEElasticMaterialPoint; }
	FEMaterialPoint* Copy()
	{
		FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint;
		pt->m_w = m_w;
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		if (bflag)
		{
			for (int i=0; i<(int) m_w.size(); ++i) m_w[i] = 1.0;
		}
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_w;
		}
		else
		{
			ar >> m_w;
		}
	}

public:
	vector<double>	m_w;	//!< material weights
};

//-----------------------------------------------------------------------------
//! Elastic mixtures

//! This class describes a mixture of elastic solids.  The user must declare
//! elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.

class FEElasticMixture : public FEElasticMaterial
{
public:
	FEElasticMixture();

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// return number of materials
	int Materials() const { return m_pMat.size(); }

	// return a material component
	FEElasticMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEElasticMaterial* pm) { m_pMat.push_back(pm); }

	// get a material parameter
	FEParam* GetParameter(const ParamString& s);
	
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();

private:
	vector <FEElasticMaterial*>	m_pMat;	//!< pointers to elastic materials
};
