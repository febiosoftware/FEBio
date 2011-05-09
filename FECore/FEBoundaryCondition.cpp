#include "stdafx.h"
#include "FEBoundaryCondition.h"

int FEBoundaryCondition::m_ncount = 0;

FEBoundaryCondition::FEBoundaryCondition()
{
	m_bactive = true; 
	m_nID = m_ncount++;
}

void FEBoundaryCondition::SetID(int nid)
{
	m_nID = nid;
	if (nid > m_ncount) m_ncount = nid+1;
}
