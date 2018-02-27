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

//=================================================================================================

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
