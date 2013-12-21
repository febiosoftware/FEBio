#include "stdafx.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::FEBiphasicContactSurface(FEMesh* pm) : FEContactSurface(pm)
{
}

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::~FEBiphasicContactSurface()
{
}

//-----------------------------------------------------------------------------
vec3d FEBiphasicContactSurface::GetFluidForce()
{
	return vec3d(0,0,0);
}

