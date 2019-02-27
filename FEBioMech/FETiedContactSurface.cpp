#include "stdafx.h"
#include "FETiedContactSurface.h"
#include "FECore/FEShellDomain.h"
#include "FECore/vector.h"
#include <FECore/FEMesh.h>

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
	m_Tc.resize(nn);		// contact forces
	m_off.resize(nn);		// surface offset values

	// set initial values
	zero(m_gap);
	zero(m_Lm);
	zero(m_off);
	zero(m_Tc);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	if (m_boffset)
	{
		FEMesh& m = *m_pMesh;
		vector<double> tag(m.Nodes());
		zero(tag);
		for (int nd = 0; nd<m.Domains(); ++nd)
		{
			FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
			if (psd)
			{
				for (int i = 0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					int n = el.Nodes();
					for (int j = 0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
				}
			}
		}
		for (int i = 0; i<nn; ++i) m_off[i] = tag[NodeIndex(i)];
	}

	return true;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::Serialize(DumpStream &ar)
{
	FEContactSurface::Serialize(ar);
	ar & m_Lm & m_gap & m_Tc;
	if (ar.IsShallow() == false)
	{
		ar & m_rs;
		ar & m_off;
	}
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ne = el.Nodes();
    pt = vec3d(0,0,0);
    for (int k=0; k<ne; ++k) pt += m_Tc[el.m_lnode[k]];
    pt /= ne;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	for (int j= 0; j< ne; ++j) pn[j] = m_Tc[f.m_lnode[j]].norm();
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	for (int j= 0; j< ne; ++j) 
	{
		tn[j] = m_Tc[f.m_lnode[j]];
	}
}
