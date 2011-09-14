// DataStore.h: interface for the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
#define AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include "FECore/DumpFile.h"
#include <vector>
using namespace std;

class FEM;
class FENodeSet;

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
	DataRecord(FEM* pfem, const char* szfile);
	virtual ~DataRecord();

	bool Write();

	void SetItemList(const char* szlist);

	virtual double Evaluate(int item, int ndata) = 0;

	virtual void SelectAllItems() = 0;

	virtual void Serialize(DumpFile& ar);

	virtual void Parse(const char* sz) = 0;

	void SetName(const char* sz);
	void SetDelim(const char* sz);
	void SetComments(bool b) { m_bcomm = b; }

public:
	int		m_nid;					//!< ID of data record
	vector<int>	m_item;				//!< item list
	vector<int>	m_data;				//!< data list

protected:
	bool	m_bcomm;				//!< export comments or not
	char	m_szname[MAX_STRING];	//!< name of expression
	char	m_szdelim[MAX_DELIM];	//!< data delimitor

protected:
	char	m_szfile[MAX_STRING];	//!< file name of data record

	FEM*	m_pfem;
	FILE*	m_fp;
};

//-----------------------------------------------------------------------------

class NodeDataRecord : public DataRecord
{
	enum { X, Y, Z, UX, UY, UZ, VX, VY, VZ, RX, RY, RZ, P, C };

public:
	NodeDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
	void SetItemList(FENodeSet* pns);
};

class ElementDataRecord : public DataRecord
{
	enum {X, Y, Z, J, 
		EX, EY, EZ, EXY, EYZ, EXZ, 
		E1, E2, E3,
		SX, SY, SZ, SXY, SYZ, SXZ,
		S1, S2, S3,
		FX, FY, FZ, FYZ, FZX, FXY, FYX, FXZ, FZY, 
		P, WX, WY, WZ, C, JX, JY, JZ};

	struct ELEMREF
	{
		int	ndom;
		int	nid;
	};

public:
	ElementDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();

protected:
	void BuildELT();

protected:
	vector<ELEMREF>	m_ELT;
};

class RigidBodyDataRecord : public DataRecord
{
	enum {X, Y, Z, QX, QY, QZ, QW, FX, FY, FZ, MX, MY, MZ};

public:
	RigidBodyDataRecord(FEM* pfem, const char* szfile) :  DataRecord(pfem, szfile){}
	double Evaluate(int item, int ndata);
	void Parse(const char* sz);
	void SelectAllItems();
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

	void Serialize(DumpFile& ar);

protected:
	vector<DataRecord*>	m_data;	//!< the data records
};

#endif // !defined(AFX_DATASTORE_H__FC7861A3_2B1A_438C_AC7D_7ADD2F8DE6F4__INCLUDED_)
