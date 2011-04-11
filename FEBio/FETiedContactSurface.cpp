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
	FEContactSurface::Init();

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

		int ne = (int) m_pme.size();
		ar << ne;
		for (int i=0; i<ne; ++i)
		{
			FESurfaceElement* pe = m_pme[i];
			if (pe) ar << pe->m_lid; else ar << -1;
		}
	}
	else
	{
		ar >> gap;
		ar >> rs;
		ar >> Lm;

		assert(m_pSibling);

		int ne = 0, id;
		ar >> ne;
		assert(ne == m_pme.size());
		for (int i=0; i<ne; ++i)
		{
			ar >> id;
			if (id < 0) m_pme[i] = 0; 
			else 
			{
				m_pme[i] = &m_pSibling->Element(id);
				assert(m_pme[i]->m_lid == id);
			}
		}
	}
}
