#include "stdafx.h"
#include "FETiedContactSurface.h"
#include "FECore/FEShellDomain.h"
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
		for (int i = 0; i<nn; ++i) m_off[i] = tag[m_node[i]];
	}

	return true;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Lm << m_gap << m_Tc;
	}
	else
	{
		dmp >> m_Lm >> m_gap >> m_Tc;
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
		ar << m_off;
		ar << m_Tc;
	}
	else
	{
		ar >> m_gap;
		ar >> m_rs;
		ar >> m_Lm;
		ar >> m_off;
		ar >> m_Tc;
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
	for (int j= 0; j< ne; ++j) 
	{
		tn[j] = m_Tc[f.m_lnode[j]];
	}
}
