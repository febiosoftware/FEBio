#pragma once
#include "FE_enum.h"
#include "FEDataStream.h"

//-----------------------------------------------------------------------------
// This class is used by domain classes to define their data exports.
// This is part of an experimental feature that allows domain classes to handle
// the data that they want to export. 
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

	virtual void Serialize(FEDataStream& d);

public:
	Var_Type	m_type;
	Storage_Fmt	m_fmt;
	void*		m_pd;		//!< pointer to data field
	const char*	m_szname;
};

#define EXPORT_DATA(iType, iFmt, pVar, Name) AddDataExport(new FEDataExport(iType, iFmt, pVar, Name));