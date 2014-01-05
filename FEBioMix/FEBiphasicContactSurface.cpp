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
void FEBiphasicContactSurface::GetNodalPressureGap(int nface, double* pg) { assert(false); }

//-----------------------------------------------------------------------------
vec3d FEBiphasicContactSurface::GetFluidForce()
{
	assert(false);
    return vec3d(0,0,0);
}
