#include "stdafx.h"
#include "FETiedContactSurface.h"
#include "FECore/vector.h"

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
	m_gap.resize(nn);		// gap funtion
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));	// penetrated master element
	m_rs.resize(nn);		// natural coords of projected slave node on master element
	m_Lm.resize(nn);		// Lagrangian multipliers

	// set initial values
	zero(m_gap);
	zero(m_Lm);

	return true;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Lm << m_gap;
	}
	else
	{
		dmp >> m_Lm >> m_gap;
	}
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::Serialize(DumpFile &ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_gap;
		ar << m_rs;
		ar << m_Lm;
	}
	else
	{
		ar >> m_gap;
		ar >> m_rs;
		ar >> m_Lm;
	}
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactGap(int nface, double* gn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.m_lnode.size();
	for (int j= 0; j< ne; ++j) gn[j] = m_gap[f.m_lnode[j]].norm();
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.m_lnode.size();
	for (int j= 0; j< ne; ++j) pn[j] = m_Lm[f.m_lnode[j]].norm();
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.m_lnode.size();
	for (int j= 0; j< ne; ++j) tn[j] = m_Lm[f.m_lnode[j]];
}
