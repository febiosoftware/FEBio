#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! material class for springs
class FEDiscreteMaterial : public FEMaterial
{
public:
	virtual double force    (double dl) = 0;
	virtual double stiffness(double dl) = 0;
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

public:
	double			m_F;	// force scale factor
	int				m_nlc; // load curve ID
	FELoadCurve*	m_plc; // force-displacement curve
};
