#pragma once
#include "FE_enum.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class FEDataExport
{
public:
	FEDataExport(Var_Type itype, Storage_Fmt ifmt, void* pd, const char* szname)
	{
		m_pd = pd;
		m_type = itype;
		m_fmt = ifmt;
		m_szname = szname;
	}

	virtual ~FEDataExport(){}

	virtual void Serialize(vector<float>& d);

public:
	Var_Type	m_type;
	Storage_Fmt	m_fmt;
	void*		m_pd;		//!< pointer to data field
	const char*	m_szname;
};
