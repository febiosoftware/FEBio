#include "stdafx.h"
#include "FEDGFluxJump.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
//! constructor
FEDGFluxJump::FEDGFluxJump(FEModel* pfem) : FESurfaceLoad(pfem)
{

}

//-----------------------------------------------------------------------------
//! construct the surface representing all inter-element boundaries.
bool FEDGFluxJump::Init()
{
	if (FESurfaceLoad::Init() == false) return false;

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// create the internal element boundary surface
	FESurface* ps = mesh.ElementBoundarySurface(false, true);

	return true;
}

//-----------------------------------------------------------------------------
//! calculate residual
void FEDGFluxJump::Residual(FEGlobalVector& R)
{
}

//-----------------------------------------------------------------------------
//! calculate stiffness matrix
void FEDGFluxJump::StiffnessMatrix(FESolver* psolver)
{
	// TODO: implement this
}
