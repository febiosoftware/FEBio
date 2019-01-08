// DataStore.h: interface for the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
#define AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DataRecord.h"
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class FECORE_API DataStore
{
public:
	DataStore();
	virtual ~DataStore();

	void Clear();

	void Write();

	void AddRecord(DataRecord* prec);

	int Size() { return (int) m_data.size(); }

	DataRecord* GetDataRecord(int i) { return m_data[i]; }

protected:
	std::vector<DataRecord*>	m_data;	//!< the data records
};


#endif // !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
