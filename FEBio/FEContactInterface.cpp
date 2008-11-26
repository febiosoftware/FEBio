// FEContactInterface.cpp: implementation of the FEContactInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEContactInterface.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEContactInterface::FEContactInterface(FEM* pfem)
{
	m_pfem = pfem;
	m_ntype = 0;
	m_nID = -1;
}

FEContactInterface::~FEContactInterface()
{

}
