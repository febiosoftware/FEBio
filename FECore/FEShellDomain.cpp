#include "stdafx.h"
#include "FEShellDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! constructor
FEShellDomain::FEShellDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_SHELL, pm)
{
}

//-----------------------------------------------------------------------------
void FEShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEShellElement& el = Element(i);
		int n = el.GaussPoints();
		for (int j = 0; j<n; ++j) el.GetMaterialPoint(j)->Update(timeInfo);
	}
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEShellElement& el = Element(i);
		int n = el.GaussPoints();
		for (int j = 0; j<n; ++j) el.GetMaterialPoint(j)->Init();
	}
}

//=================================================================================================

FEShellDomainOld::FEShellDomainOld(FEMesh* pm) : FEShellDomain(pm) 
{
}

//-----------------------------------------------------------------------------
void FEShellDomainOld::Create(int nelems, int elemType)
{
	m_Elem.resize(nelems);
	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
void FEShellDomainOld::Serialize(DumpStream &ar)
{
	if (ar.IsShallow())
	{
		int NEL = (int) m_Elem.size();
		for (int i=0; i<NEL; ++i)
		{
			FEShellElement& el = m_Elem[i];
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
            el.Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_Node;
		
			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FEShellElementOld& el = m_Elem[i];
				ar << el.Type();

				ar << el.GetMatID();
				ar << el.GetID();
				ar << el.m_node;

				ar << el.m_h0;
				ar << el.m_D0;

				for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			int n, mat, nid;

			ar >> m_Node;

			FEMaterial* pmat = GetMaterial();
			assert(pmat);

			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FEShellElementOld& el = m_Elem[i];
				ar >> n;

				el.SetType(n);

				ar >> mat; el.SetMatID(mat);
				ar >> nid; el.SetID(nid);
				ar >> el.m_node;

				ar >> el.m_h0;
				ar >> el.m_D0;

				for (int j=0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate all shell normals (i.e. the shell directors).
//! And find shell nodes
void FEShellDomainOld::InitShells(FEMesh& mesh)
{
	// zero initial directors for shell nodes
	int NN = mesh.Nodes();
	vector<vec3d> D(NN, vec3d(0, 0, 0));

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (mesh.Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomainOld& sd = static_cast<FEShellDomainOld&>(mesh.Domain(nd));
			vec3d r0[FEElement::MAX_NODES];
			for (int i = 0; i<sd.Elements(); ++i)
			{
				FEShellElementOld& el = sd.ShellElement(i);

				int n = el.Nodes();
				int* en = &el.m_node[0];

				// get the nodes
				for (int j = 0; j<n; ++j) r0[j] = mesh.Node(en[j]).m_r0;

				for (int j = 0; j<n; ++j)
				{
					int m0 = j;
					int m1 = (j + 1) % n;
					int m2 = (j == 0 ? n - 1 : j - 1);

					vec3d a = r0[m0];
					vec3d b = r0[m1];
					vec3d c = r0[m2];

					D[en[m0]] += (b - a) ^ (c - a);
				}
			}
		}
	}

	// make sure we start with unit directors
	for (int i = 0; i<NN; ++i) D[i].unit();

	// assign directors to shells 
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (mesh.Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomainOld& sd = static_cast<FEShellDomainOld&>(mesh.Domain(nd));
			for (int i = 0; i<sd.Elements(); ++i)
			{
				FEShellElementOld& el = sd.ShellElement(i);
				int ne = el.Nodes();
				for (int j = 0; j<ne; ++j) el.m_D0[j] = D[el.m_node[j]] * el.m_h0[j];
			}
		}
	}
}

//=================================================================================================

FEShellDomainNew::FEShellDomainNew(FEMesh* pm) : FEShellDomain(pm) 
{
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::Create(int nelems, int elemType)
{
	m_Elem.resize(nelems);
	if (elemType != -1)
		for (int i = 0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::Serialize(DumpStream &ar)
{
	if (ar.IsShallow())
	{
		int NEL = (int)m_Elem.size();
		for (int i = 0; i<NEL; ++i)
		{
			FEShellElement& el = m_Elem[i];
			int nint = el.GaussPoints();
			for (int j = 0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
			el.Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_Node;

			for (size_t i = 0; i<m_Elem.size(); ++i)
			{
				FEShellElement& el = m_Elem[i];
				ar << el.Type();

				ar << el.GetMatID();
				ar << el.GetID();
				ar << el.m_node;

				ar << el.m_h0;

				for (int j = 0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			int n, mat, nid;

			ar >> m_Node;

			FEMaterial* pmat = GetMaterial();
			assert(pmat);

			for (size_t i = 0; i<m_Elem.size(); ++i)
			{
				FEShellElement& el = m_Elem[i];
				ar >> n;

				el.SetType(n);

				ar >> mat; el.SetMatID(mat);
				ar >> nid; el.SetID(nid);
				ar >> el.m_node;

				ar >> el.m_h0;

				for (int j = 0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEShellDomainNew::InitShells(FEMesh& mesh)
{
	// zero initial directors for shell nodes
	int NN = mesh.Nodes();
	vector<vec3d> D(NN, vec3d(0, 0, 0));
	vector<int> ND(NN, 0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// Calculate the shell directors as the local node normals
		if (mesh.Domain(nd).Class() == FE_DOMAIN_SHELL)
		{
			FEShellDomain& sd = static_cast<FEShellDomain&>(mesh.Domain(nd));
			vec3d r0[FEElement::MAX_NODES];
			for (int i = 0; i<sd.Elements(); ++i)
			{
				FEShellElement& el = sd.Element(i);

				int n = el.Nodes();
				int* en = &el.m_node[0];

				// get the nodes
				for (int j = 0; j<n; ++j) r0[j] = mesh.Node(en[j]).m_r0;

				for (int j = 0; j<n; ++j)
				{
					int m0 = j;
					int m1 = (j + 1) % n;
					int m2 = (j == 0 ? n - 1 : j - 1);

					vec3d a = r0[m0];
					vec3d b = r0[m1];
					vec3d c = r0[m2];
					vec3d d = (b - a) ^ (c - a); d.unit();

					D[en[m0]] += d*el.m_h0[j];
					++ND[en[m0]];
				}
			}
		}
	}

	// assign initial directors to shell nodes
	// make sure we average the directors
	for (int i = 0; i<NN; ++i)
		if (ND[i] > 0) mesh.Node(i).m_d0 = D[i] / ND[i];
}
