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

	//! calculate the conductivity gradient
	tens4ds ConductivityGradient(FEMaterialPoint& pt);

private: // material properties
	FEPropertyT<FEThermalElastic>			m_pElastic;	//!< elastic material property
	FEPropertyT<FEThermalConductivity>		m_pCond;	//!< thermal conductivity property
};
