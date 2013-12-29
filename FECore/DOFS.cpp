//
//  DOFS.cpp
//  FECore
//
//  Created by Gerard Ateshian on 12/28/13.
//  Copyright (c) 2013 febio.org. All rights reserved.
//

#include "DOFS.h"

//=============================================================================
// The one-and-only DOFS
//DOFS& fedofs = *DOFS::GetInstance();

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
    MAX_NDOFS = 11;
    MAX_CDOFS = 0;
}

//-----------------------------------------------------------------------------
// destructor for the DOFS class
//
DOFS::~DOFS()
{
	m_pdofs = 0;
}

//-----------------------------------------------------------------------------
// set total number of DOFS
//
void DOFS::SetNDOFS(int ndofs)
{
    MAX_NDOFS = ndofs;
}

//-----------------------------------------------------------------------------
// return total number of DOFS
//
int DOFS::GetNDOFS()
{
    return MAX_NDOFS;
}

//-----------------------------------------------------------------------------
// set number of solute DOFS
//
void DOFS::SetCDOFS(int cdofs)
{
    MAX_CDOFS = cdofs;
}

//-----------------------------------------------------------------------------
// return number of solute DOFS
//
int DOFS::GetCDOFS()
{
    return MAX_CDOFS;
}
