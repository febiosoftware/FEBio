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

//-----------------------------------------------------------------------------
// Max nr of nodal degrees of freedom

// At this point the nodal dofs are used as follows:
//
#define DOF_X			0		// x-displacement
#define DOF_Y			1		// y-displacement
#define DOF_Z			2		// z-displacement
#define DOF_U			3		// x-rotation
#define DOF_V			4		// y-rotation
#define DOF_W			5		// z-rotation
#define DOF_P			6		// fluid pressure
#define DOF_RU			7		// rigid x-rotation
#define DOF_RV			8		// rigid y-rotation
#define DOF_RW			9		// rigid z-rotation
#define DOF_T			10		// temperature
#define DOF_VX			11		// x-fluid velocity
#define DOF_VY			12		// y-fluid velocity
#define DOF_VZ			13		// z-fluid velocity
#define DOF_E           14      // fluid dilatation
#define DOF_C			15		// solute concentration
//
// The rotational degrees of freedom are only used for rigid nodes and shells.
// The fluid pressure is only used for poroelastic problems.
// The rigid rotational degrees of freedom are only used for rigid nodes and only during the creation of the stiffenss matrix
// The temperature is only used during heat-conduction problems
// The solute concentration is only used in solute transport problems.

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
	else if (strcmp(sz, "Ru") == 0) bc = DOF_RU;
	else if (strcmp(sz, "Rv") == 0) bc = DOF_RV;
	else if (strcmp(sz, "Rw") == 0) bc = DOF_RW;
	else if (strcmp(sz, "t" ) == 0) bc = DOF_T; 
    else if (strcmp(sz, "vx") == 0) bc = DOF_VX;
    else if (strcmp(sz, "vy") == 0) bc = DOF_VY;
    else if (strcmp(sz, "vz") == 0) bc = DOF_VZ;
    else if (strcmp(sz, "e" ) == 0) bc = DOF_E;
	else if (strcmp(sz, "c" ) == 0) bc = DOF_C;
	else if (sz[0] == 'c') bc = DOF_C + atoi(&sz[1]) - 1;
	return bc;
}
