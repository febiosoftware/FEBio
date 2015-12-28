#pragma once
#include "FECore/FEDiscreteMaterial.h"
#include "FECore/LoadCurve.h"		

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
	double force    (double dl);
	double stiffness(double dl);

public:
	double m_E;	//!< spring constant

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! tension-only linear spring
class FETensionOnlyLinearSpring : public FESpringMaterial
{
public:
	FETensionOnlyLinearSpring(FEModel* pfem) : FESpringMaterial(pfem){}
	double force    (double dl);
	double stiffness(double dl);
	bool Init();

public:
	double m_E;	//!< spring constant

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! general purpose nonlinear spring
class FENonLinearSpring : public FESpringMaterial
{
public:
	FENonLinearSpring(FEModel* pfem);

	double force    (double dl);
	double stiffness(double dl);
	bool Init();

	void Serialize(DumpFile& ar);

	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval);

public:
	double			m_F;	// force scale factor
	int				m_nlc; // load curve ID
	FELoadCurve*	m_plc; // force-displacement curve

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
