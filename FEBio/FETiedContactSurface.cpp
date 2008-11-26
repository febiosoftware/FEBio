#include "StdAfx.h"
#include "FETiedContactSurface.h"


//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

void FETiedContactSurface::Init()
{
	// always intialize base class first!
	FESurface::Init();

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	gap.create(nn);		// gap funtion
	pme.create(nn);		// penetrated master element
	rs.create(nn);		// natural coords of projected slave node on master element
	Lm.create(nn);		// Lagrangian multipliers

	// set initial values
	gap.zero();
	pme.set(0);
	Lm.zero();
}
