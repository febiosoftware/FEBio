#pragma once
#include "FEUncoupledMaterial.h"
#include "FEUncoupledFiberExpLinear.h"

//-----------------------------------------------------------------------------
class FEMRVonMisesMaterialPoint : public FEMaterialPoint
{
public:
	FEMRVonMisesMaterialPoint(FEMaterialPoint* mp = 0);

	FEMaterialPoint* Copy() override;

	void Serialize(DumpStream& ar) override;

public:
	double	m_kf;
	double	m_tp;
};

//-----------------------------------------------------------------------------
//! Transversely Isotropic Multiple material

//! This material has an isotopric Multiple basis and single preferred
//! fiber direction.

class FEMRVonMisesFibers: public FEUncoupledMaterial
{
public:
	FEMRVonMisesFibers (FEModel* pfem);

public:
	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2
	double	kf;	//!< Fiber Concentration Factor
	double	tp;	//!< Preferred Fiber Orientation IN RADIANS
	int	gipt;	//!< Gauss Integration Points
	int	vmc;	//!< Choice of von Mises distribution 
	double var_n;	//!< Exponent for the constrained von Mises distribution

public:
	//! calculate stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData() override;

	// declare parameter list
	DECLARE_FECORE_CLASS();

protected:
	FEUncoupledFiberExpLinear	m_fib;
};
