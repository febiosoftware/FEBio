#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include "FEHeatTransferMaterial.h"
#include "FEThermalElastic.h"
#include "FEThermalConductivity.h"

//-----------------------------------------------------------------------------
//! Class describing a thermo-elastic material
class FEThermoElasticMaterial : public FEMaterial
{
public:
	FEThermoElasticMaterial(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pElastic->GetElasticMaterial(); }

	// find a material parameter
	FEParam* GetParameter(const ParamString& s);
	
public:
	void Init();
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);

	//! return the thermal conductivity
	mat3ds Conductivity(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate the spatial thermal tangent (derivative of stress w.r.t. temperature)
	mat3ds ThermalTangent(FEMaterialPoint& pt);

	//! Serialization
	void Serialize(DumpFile& ar);

public:
	//! return number of material properties
	int Properties();

	//! return a material property
	FECoreBase* GetProperty(int n);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);

private: // material properties
	FEThermalElastic*			m_pElastic;	//!< elastic material property
	FEThermalConductivity*		m_pCond;	//!< thermal conductivity property
};
