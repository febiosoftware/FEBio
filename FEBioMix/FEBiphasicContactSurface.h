#pragma once
#include "FECore/FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class describes a contact surface used in a biphasic/multiphasic analysis.
class FEBiphasicContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEBiphasicContactSurface(FEMesh* pm = 0);

	//! destructor
	~FEBiphasicContactSurface();

public:
	//! Get the total force exerted by the fluid
    virtual vec3d GetFluidForce();
};
