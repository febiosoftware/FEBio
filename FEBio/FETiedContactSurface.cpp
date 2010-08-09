#include "stdafx.h"
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
	gap.resize(nn);		// gap funtion
	pme.assign(nn, NULL);	// penetrated master element
	rs.resize(nn);		// natural coords of projected slave node on master element
	Lm.resize(nn);		// Lagrangian multipliers

	// set initial values
	zero(gap);
	zero(Lm);
}
