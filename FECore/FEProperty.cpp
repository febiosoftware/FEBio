#include "stdafx.h"
#include "FEProperty.h"
#include "FECoreBase.h"
#include "DumpStream.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
FEProperty::FEProperty(int classID) : m_szname(nullptr), m_flags(0), m_classID(classID) {}

//-----------------------------------------------------------------------------
FEProperty::~FEProperty(){}

//-----------------------------------------------------------------------------
//! Set the name of the property.
//! Note that the name is not copied so it must point to a static string.
FEProperty& FEProperty::SetName(const char* sz)
{
	m_szname = sz;
	return *this;
}

//-----------------------------------------------------------------------------
//! Return the name of this property
const char* FEProperty::GetName() const { return m_szname; }

//-----------------------------------------------------------------------------
void FEProperty::Write(DumpStream& ar, FECoreBase* pc)
{
	int nflag = (pc == 0 ? 0 : 1);
	ar << nflag;
	if (nflag)
	{
		int ntype = pc->GetSuperClassID();
		ar << pc->GetTypeStr();
		ar << ntype;
		pc->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
FECoreBase* FEProperty::Read(DumpStream& ar)
{
	int nflag = 0;
	FECoreBase* pm = 0;
	ar >> nflag;
	if (nflag)
	{
		char sz[256];
		int ntype = 0;
		ar >> sz;
		ar >> ntype;
		pm = fecore_new<FECoreBase>(ntype, sz, &ar.GetFEModel());
		pm->SetParent(GetParent());
		pm->Serialize(ar);

		// TODO: Do I really need to do this here?
		//pm->Init();
	}
	return pm;
}
