#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! material class for springs
class FEDiscreteMaterial : public FEMaterial
{
public:
	virtual double force    (double dl) = 0;
	virtual double stiffness(double dl) = 0;

	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEDiscreteMaterialPoint; }
};

//-----------------------------------------------------------------------------
//! linear spring
class FELinearSpring : public FEDiscreteMaterial
{
public:
	double force    (double dl);
	double stiffness(double dl);
	void Init();

public:
	double m_E;	//!< spring constant

	// declare as registered
	DECLARE_REGISTERED(FELinearSpring);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! tension-only linear spring
class FETensionOnlyLinearSpring : public FEDiscreteMaterial
{
public:
	double force    (double dl);
	double stiffness(double dl);
	void Init();

public:
	double m_E;	//!< spring constant

	// declare as registered
	DECLARE_REGISTERED(FETensionOnlyLinearSpring);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! general purpose nonlinear spring
class FENonLinearSpring : public FEDiscreteMaterial
{
public:
	FENonLinearSpring() { m_nlc = -1; m_F = 1; }

	double force    (double dl);
	double stiffness(double dl);
	void Init();

	void Serialize(DumpFile& ar);

public:
	double			m_F;	// force scale factor
	int				m_nlc; // load curve ID
	FELoadCurve*	m_plc; // force-displacement curve

	// declare as registered
	DECLARE_REGISTERED(FENonLinearSpring);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
