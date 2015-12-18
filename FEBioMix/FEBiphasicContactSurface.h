#pragma once
#include "FEBioMech/FEContactSurface.h"

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
	//! Get the fluid pressure gap
	virtual void GetNodalPressureGap(int nface, double* pg);
    
	//! Get the total force exerted by the fluid
    virtual vec3d GetFluidForce();

protected:
	int	m_dofP;
};
