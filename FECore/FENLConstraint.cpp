#include "stdafx.h"
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint(FEModel* pfem) : FEModelComponent(FENLCONSTRAINT_ID, pfem)
{
	static int ncount = 1;
	SetID(ncount++);
}

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint()
{
}
