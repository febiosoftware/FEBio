#include "stdafx.h"
#include "FENLConstraint.h"


int FENLConstraint::m_ncount = 0;

//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint(FEModel* pfem) : FEModelComponent(FENLCONSTRAINT_ID, pfem)
{
	m_nID = m_ncount++;
}

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint(){}
