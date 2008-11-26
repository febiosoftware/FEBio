// MyMap.h: interface for the CMyMap class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYMAP_H__D377BA90_B934_4751_BA46_EF8A15F55061__INCLUDED_)
#define AFX_MYMAP_H__D377BA90_B934_4751_BA46_EF8A15F55061__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string.h>

#define MAXNAME 32	// max length of variable name

template <class T>
class CMyMap  
{
protected:

	struct List
	{
		List(const char* szvar, const T& val) // constructor
		{
			strcpy(szname, szvar);
			value = val;

			pnext = 0;
		}

		char	szname[MAXNAME];	// name of variable
		T		value;				// value of variable
		List*	pnext;
	};

public:
	CMyMap(int nsize);
	virtual ~CMyMap();
	void Clear();

	T& operator [] (const char* szname); // retrieval of value (and insertion of element)

	int Find(const char* szname);	// returns the index of the key in the table, or -1 if not found
	
protected:
	unsigned long hash(const char* str);	// the hash function

protected:

	int		m_nsize;	// size of array
	List**	m_pTable;	// table of linked lists
};

template <class T>
CMyMap<T>::CMyMap(int nsize)
{	
	// create the new list
	m_pTable = new List* [nsize];
	m_nsize = nsize;

	// initialize with zeroes
	for (int i=0; i<nsize; i++) m_pTable[i] = 0;
}


template <class T>
CMyMap<T>::~CMyMap()
{
	Clear();
}


template <class T>
void CMyMap<T>::Clear()
{
	// delete the linked lists
	for (int i=0; i<m_nsize; i++) 
	{
		List* plist = m_pTable[i];
		while (plist)
		{
			List* ptemp = plist;
			plist = plist->pnext;
			delete ptemp;
		}
	}

	// delete the table
	if (m_nsize > 0) delete [] m_pTable;

	m_pTable = 0;
	m_nsize = 0;
}


template <class T>
T& CMyMap<T>::operator [](const char* szname)
{
	// hash the string
	unsigned long n = hash(szname) % m_nsize;

	// get the table list
	List* pl = m_pTable[n];

	if (pl == 0) 
	{
		// if there is no entry here, make one
		pl = m_pTable[n] = new List(szname, T(0));
	}
	else
	{
		// see if we can find the entry
		List* plast=0;
		while (pl) { if (strcmp(pl->szname, szname)==0) break; else { plast = pl; pl = pl->pnext; }}

		// if none is found, add one to the list
		if (pl == 0) pl = plast->pnext = new List(szname, T(0));
	}

	return pl->value;
}

template <class T> int CMyMap<T>::Find(const char* str)
{
	// hash the string
	unsigned long n = hash(str) % m_nsize;

	// get the table list
	List* pl = m_pTable[n];

	if (pl == 0) return -1;
	else
	{
		// see if we can find the entry
		List* plast;
		while (pl) { if (strcmp(pl->szname, str)==0) break; else { plast = pl; pl = pl->pnext; }}

		// if none is found, add one to the list
		if (pl == 0) return -1; else return n;
	}
}

template <class T>
unsigned long CMyMap<T>::hash(const char* str)
{
	unsigned long n = 5381;
	int c; 
	while (c = *str++) n = ((n << 5) + n) + c; // n*33 + c

	return n;
}


#endif // !defined(AFX_MYMAP_H__D377BA90_B934_4751_BA46_EF8A15F55061__INCLUDED_)
