#pragma once

#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface

//!	this class is used in contact analyses to describe a contacting
//! surface in a contact interface.

class FEContactSurface : public FESurface
{
public:
	//! constructor
	FEContactSurface(FEMesh* pm=0) : FESurface(pm) {}

	//! data initialization
	void Init();
};
