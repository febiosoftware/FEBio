#include "stdafx.h"
#include "FEParameterList.h"
#include <cstring>
#include <assert.h>

//-----------------------------------------------------------------------------
// This function adds a parameter to the parameter list
void FEParameterList::AddParameter(void *pv, FEParamType itype, int ndim, const char *sz)
{
	// create a new parameter object
	FEParam p;

	// set the pointer to the value
	assert(pv);
	p.m_pv = pv;

	// set the type
	p.m_itype = itype;

	// set the dimension
	// make sure that the type is a vector type if the dimension > 1
	assert((ndim >= 1) && ( (ndim > 1 ? (itype >= 100) : true)));
	p.m_ndim = ndim;

	// set the name
	// note that we just copy the pointer, not the actual string
	// this is okay as long as the name strings are defined
	// as literal strings
	assert(sz);
	p.m_szname = sz;

	// add the parameter to the list
	m_pl.push_back(p);
}

//-----------------------------------------------------------------------------
// This function searches the parameters in the list for a parameter
// with the name given by the input argument
// \param sz name of parameter to find
FEParam* FEParameterList::Find(const char* sz)
{
	FEParam* pp = 0;
	if (m_pl.size() > 0)
	{
		list<FEParam>::iterator it;
		for (it = m_pl.begin(); it != m_pl.end(); ++it)
		{
			if (strcmp(it->m_szname, sz) == 0)
			{
				pp = &(*it);
				break;
			}
		}
	}

	return pp;
}

//-----------------------------------------------------------------------------
FEParameterList& FEParamContainer::GetParameterList()
{
	if (m_pParam == 0) 
	{
		m_pParam = new FEParameterList;
		BuildParamList();
	}
	return *m_pParam;
}
