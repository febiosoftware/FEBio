// DataStore.h: interface for the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
#define AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include "FECore/vector.h"
#include "MathParser.h"

class FEM;
class FENodeSet;

//-----------------------------------------------------------------------------

class DataRecord
{
public:
	enum {MAX_DELIM=16, MAX_STRING=128};
public:
	DataRecord(FEM* pfem, const char* szfile);
	virtual ~DataRecord();

	bool Write();

	void SetItemList(const char* szlist);

	virtual double Evaluate(int item, const char* szexpr) = 0;

	virtual void SelectAllItems() = 0;

public:
	int		m_nid;					//!< ID of data record
	char	m_szdata[MAX_STRING];	//!< expression of data record
	char	m_szname[MAX_STRING];	//!< name of expression
	char	m_szdelim[MAX_DELIM];	//!< delimiter used
	vector<int>	m_item;				//!< item list
	MathParser	m_calc;
	bool	m_bcomm;			//!< export comments or not

protected:
	char	m_szfile[MAX_STRING];	//!< file name of data record

	FEM*	m_pfem;
	FILE*	m_fp;
};

//-----------------------------------------------------------------------------

class NodeDataRecord : public DataRecord
{
public:
	NodeDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, const char* szexpr);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
};

class ElementDataRecord : public DataRecord
{
public:
	ElementDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, const char* szexpr);
	void SelectAllItems();
};

class RigidBodyDataRecord : public DataRecord
{
public:
	RigidBodyDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, const char* szexpr);
	void SelectAllItems();
};

//-----------------------------------------------------------------------------

class DataStore  
{
public:
	DataStore();
	virtual ~DataStore();

	void Write();

	void AddRecord(DataRecord* prec);

protected:
	ptr_vector<DataRecord>	m_data;	//!< the data records
};

#endif // !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
