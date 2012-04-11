#include "stdafx.h"
#include "FETiedContactSurface.h"
#include "NumCore/vector.h"

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FETiedContactSurface::Init()
{
	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	gap.resize(nn);		// gap funtion
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));	// penetrated master element
	rs.resize(nn);		// natural coords of projected slave node on master element
	Lm.resize(nn);		// Lagrangian multipliers

	// set initial values
	zero(gap);
	zero(Lm);

	return true;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::Serialize(DumpFile &ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << gap;
		ar << rs;
		ar << Lm;
	}
	else
	{
		ar >> gap;
		ar >> rs;
		ar >> Lm;
	}
}
