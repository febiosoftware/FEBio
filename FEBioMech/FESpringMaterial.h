#pragma once
#include <FECore/FEDiscreteMaterial.h>
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
//! material class for discrete elements
class FESpringMaterial : public FEDiscreteMaterial
{
public:
	FESpringMaterial(FEModel* pfem) : FEDiscreteMaterial(pfem) {}

	virtual double force    (double dl) = 0;
	virtual double stiffness(double dl) = 0;
};

//-----------------------------------------------------------------------------
//! linear spring
class FELinearSpring : public FESpringMaterial
{
public:
	FELinearSpring(FEModel* pfem) : FESpringMaterial(pfem){}
	double force    (double dl) override;
	double stiffness(double dl) override;

public:
	double m_E;	//!< spring constant

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! tension-only linear spring
class FETensionOnlyLinearSpring : public FESpringMaterial
{
public:
	FETensionOnlyLinearSpring(FEModel* pfem) : FESpringMaterial(pfem){}
	double force    (double dl) override;
	double stiffness(double dl) override;

public:
	double m_E;	//!< spring constant

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! general purpose nonlinear spring
class FENonLinearSpring : public FESpringMaterial
{
public:
	FENonLinearSpring(FEModel* pfem);

	double force    (double dl) override;
	double stiffness(double dl) override;

public:
	FEFunction1D*	m_F;	//!< force-displacement function

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEExperimentalSpring : public FESpringMaterial
{
public:
	FEExperimentalSpring(FEModel* fem);

	double force(double dl) override;
	double stiffness(double dl) override;

public:
	double	m_E;
	double	m_sM, m_sm;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
