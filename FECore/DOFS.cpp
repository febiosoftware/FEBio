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
    MAX_NDOFS = 15;
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

//-----------------------------------------------------------------------------
int DOFS::GetDOF(const char* sz)
{
	int bc = -1;
	if      (strcmp(sz, "x" ) == 0) bc = DOF_X;
	else if (strcmp(sz, "y" ) == 0) bc = DOF_Y;
	else if (strcmp(sz, "z" ) == 0) bc = DOF_Z;
	else if (strcmp(sz, "u" ) == 0) bc = DOF_U;
	else if (strcmp(sz, "v" ) == 0) bc = DOF_V;
	else if (strcmp(sz, "w" ) == 0) bc = DOF_W;
	else if (strcmp(sz, "p" ) == 0) bc = DOF_P;
	else if (strcmp(sz, "t" ) == 0) bc = DOF_T; 
    else if (strcmp(sz, "vx") == 0) bc = DOF_VX;
    else if (strcmp(sz, "vy") == 0) bc = DOF_VY;
    else if (strcmp(sz, "vz") == 0) bc = DOF_VZ;
    else if (strcmp(sz, "e" ) == 0) bc = DOF_E;
	else if (strcmp(sz, "c" ) == 0) bc = DOF_C;
	else if (strcmp(sz, "c1") == 0) bc = DOF_C;
	else if (strncmp(sz, "c", 1) == 0) bc = DOF_C + atoi(&sz[1]) - 1;
	return bc;
}
