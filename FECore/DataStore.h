// DataStore.h: interface for the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
#define AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DumpFile.h"
#include "FEMesh.h"
#include "FEModel.h"
#include <stdio.h>
#include <vector>
using namespace std;


#define FE_DATA_NODE	1
#define FE_DATA_ELEM	2
#define FE_DATA_RB		3

//-----------------------------------------------------------------------------
// Exception thrown when parsing fails
class UnknownDataField 
{
public:
	UnknownDataField(const char* sz);
	char	m_szdata[64];
};

//-----------------------------------------------------------------------------

class DataRecord
{
public:
	enum {MAX_DELIM=16, MAX_STRING=128};
public:
	DataRecord(FEModel* pfem, const char* szfile);
	virtual ~DataRecord();

	bool Write();

	void SetItemList(const char* szlist);

	void SetName(const char* sz);
	void SetDelim(const char* sz);
	void SetComments(bool b) { m_bcomm = b; }

public:
	virtual double Evaluate(int item, int ndata) = 0;
	virtual void SelectAllItems() = 0;
	virtual void Serialize(DumpFile& ar);
	virtual void Parse(const char* sz) = 0;
	virtual int Size() = 0;

public:
	int			m_nid;				//!< ID of data record
	vector<int>	m_item;				//!< item list

protected:
	bool	m_bcomm;				//!< export comments or not
	char	m_szname[MAX_STRING];	//!< name of expression
	char	m_szdelim[MAX_DELIM];	//!< data delimitor
	char	m_szdata[MAX_STRING];	//!< data expression

protected:
	char	m_szfile[MAX_STRING];	//!< file name of data record

	FEModel*	m_pfem;
	FILE*		m_fp;
};

//-----------------------------------------------------------------------------
class DataStore  
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
	vector<DataRecord*>	m_data;	//!< the data records
};


#endif // !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
