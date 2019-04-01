#pragma once
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
