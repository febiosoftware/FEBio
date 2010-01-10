#pragma once
#include "FEMaterial.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial :	public FEElasticMaterial
{
public:
	FEMicroMaterial(void);
	~FEMicroMaterial(void);

public:
	char	m_szrve[256];	//!< filename for RVE file

protected:
	FEM		m_rve;	//!< the RVE (Representive Volume Element)
	double	m_V0;	//!< initial volume of RVE

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	double BulkModulus () { return 0; }

protected:
	void PrepRVE();
	mat3ds AveragedStress(FEMaterialPoint& pt);

public:
	// declare as registered
	DECLARE_REGISTERED(FEMicroMaterial);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
