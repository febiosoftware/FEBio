#pragma once
#include "FETransverselyIsotropic.h"

//-----------------------------------------------------------------------------
class FEPreStrainMaterialPoint : public FEMaterialPoint
{
public:
	FEPreStrainMaterialPoint(FEMaterialPoint* pt);

	void Init(bool bflag);

	FEMaterialPoint* Copy();

	void Serialize(DumpFile& ar);

	void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	double	m_ltrg;	//!< target fiber stretch
	double	m_lam;	//!< in-situ fiber stretch
	double	m_lamp;	//!< previous in-situ fiber stretch
};

//-----------------------------------------------------------------------------
class FEPreStrainTransIsoMR: public FETransverselyIsotropic
{
public:
	FEPreStrainTransIsoMR(FEModel* pfem) : FETransverselyIsotropic(pfem) {}

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2

	double	m_ltrg;	//!< target stretch

public:
	//! create material point data for this material
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEPreStrainMaterialPoint(new FEElasticMaterialPoint); }

	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	//! target fiber stretch
	double FiberStretch(FEMaterialPoint& pt);

protected:
	// calculate the prestrain deformation gradient
	mat3d PreStrainDeformationGradient(FEMaterialPoint& mp);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
