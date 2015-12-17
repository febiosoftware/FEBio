//
//  DOFS.cpp
//  FECore
//
//  Created by Gerard Ateshian on 12/28/13.
//  Copyright (c) 2013 febio.org. All rights reserved.
//

#include "DOFS.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
using namespace std;

//-----------------------------------------------------------------------------
DOFS* DOFS::m_pdofs = 0;

//-----------------------------------------------------------------------------
DOFS* DOFS::GetInstance()
{
	if (m_pdofs == 0) m_pdofs = new DOFS();
	return m_pdofs;
}

//-----------------------------------------------------------------------------
// constructor for the DOFS class
DOFS::DOFS()
{
	Reset();
}

//-----------------------------------------------------------------------------
void DOFS::Reset()
{
	// clear the DOFS
	if (m_dof.empty() == false) m_dof.clear();
}

//-----------------------------------------------------------------------------
// destructor for the DOFS class
//
DOFS::~DOFS()
{
	m_pdofs = 0;
}

//-----------------------------------------------------------------------------
//! Add a degree of freedom.
//! Returns -1 if the degree of freedom exists or if the symbol is invalid
//! \sa DOFS::GetDOF
int DOFS::AddDOF(const char* sz, int nsize)
{
	// Make sure sz is a valid symbol
	if (sz    == 0) return -1;	// cannot be null
	if (sz[0] == 0) return -1;	// must have non-zero length

	// See if the dof is already defined
	int ndof = -1;
	for (int i=0; i<(int)m_dof.size(); ++i)
	{
		if (strcmp(sz, m_dof[i].sz) == 0)
		{
			// it's already defined so return -1
			return -1;
		}
	}
	
	// If we get here, the variable does not exist yet
	DOF_ITEM it = {sz, nsize, 0};
	m_dof.push_back(it);

	// update all dofs
	Update();

	// return a nonnegative number
	return (int) (m_dof.size() - 1);
}

//-----------------------------------------------------------------------------
//! Return the DOF index.
//! This index is used in the FENode::get(), FENode::set() functions to set 
//! the values of nodal values. This index is also used in the FENode::m_ID and FENode::m_BC arrays.
int DOFS::GetDOF(const char* sz, int n)
{
	const int NDOF = (int)m_dof.size();
	for (int i=0; i<NDOF; ++i)
	{
		DOF_ITEM& it = m_dof[i];
		if (strcmp(it.sz, sz) == 0) 
		{
			assert((n>=0)&&(n<it.nsize));
			return it.ndof + n;
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
//! Return the symbol assigned to a variable.
//! Returns 0 if the dof is undefined.
const char* DOFS::GetDOFSymbol(int nvar)
{
	if ((nvar >= 0) && (nvar < (int) m_dof.size())) return m_dof[nvar].sz;
	else return 0;
}

//-----------------------------------------------------------------------------
//! Change the size of the dof array of a variable
bool DOFS::ChangeDOFSize(const char* sz, int nsize)
{
	int n = (int) m_dof.size();
	for (int i=0; i<n; ++i)
	{
		if (strcmp(sz, m_dof[i].sz) == 0)
		{
			m_dof[i].nsize = nsize;
			Update();
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Change the size of the dof array of a variable
int DOFS::GetDOFSize(const char* sz)
{
	int n = (int) m_dof.size();
	for (int i=0; i<n; ++i)
	{
		if (strcmp(sz, m_dof[i].sz) == 0)
		{
			return m_dof[i].nsize;
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------
void DOFS::Update()
{
	m_maxdofs = 0;
	int n = (int) m_dof.size();
	for (int i=0; i<n; ++i)
	{
		DOF_ITEM& it = m_dof[i];
		it.ndof = m_maxdofs;
		m_maxdofs += it.nsize;
	}
}
