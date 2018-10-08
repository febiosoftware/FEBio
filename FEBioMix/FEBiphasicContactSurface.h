#pragma once
#include "FEBioMech/FEContactSurface.h"

//-----------------------------------------------------------------------------
class FEBiphasicContactPoint : public FEContactMaterialPoint
{
public:
	double	m_Lmp;	//!< lagrange multipliers for fluid pressures
	double	m_pg;	//!< pressure "gap" for biphasic contact
};

//-----------------------------------------------------------------------------
//! This class describes a contact surface used in a biphasic/multiphasic analysis.
class FEBiphasicContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEBiphasicContactSurface(FEModel* pfem);

	//! destructor
	~FEBiphasicContactSurface();

	//! initialization
	bool Init();

	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! Get the total force exerted by the fluid
    virtual vec3d GetFluidForce();

    //! Get the total force exerted by the fluid
    virtual double GetFluidLoadSupport();
    
protected:
	int	m_dofP;
};
