#pragma once
#include "FEElasticMaterial.h"


class FEPRLig : public FEElasticMaterial
{
public:
	FEPRLig(FEModel* pfem);

public:
	
	double	m_c1;	 //!< fiber constant c1
	double	m_c2;	 //!< fiber constant c2 
	double  m_u;     //!< Lame Coefficient mu of the matrix
	double	m_m;	 //!< Poisson's ratio slope 
	double	m_v0;	 //!< initial Poisson's ratio
	double	m_k;	 //!< Penalty for the volumetric strain energy


public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
