#include "stdafx.h"
#include "FEContactSurface.h"
#include <assert.h>

//-----------------------------------------------------------------------------
FEContactSurface::FEContactSurface(FEMesh* pm) : FESurface(pm) { m_pSibling = 0; }

//-----------------------------------------------------------------------------
FEContactSurface::~FEContactSurface() { m_pSibling = 0; }

//-----------------------------------------------------------------------------
void FEContactSurface::SetSibling(FEContactSurface* ps) { m_pSibling = ps; }

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactGap(int nface, double* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactPressure(int nface, double* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactTraction(int nface, vec3d* pt) {}

//-----------------------------------------------------------------------------
vec3d FEContactSurface::GetContactForce() { return vec3d(0,0,0); }

//-----------------------------------------------------------------------------
double FEContactSurface::GetContactArea() { return 0; }
