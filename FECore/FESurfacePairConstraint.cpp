#include "stdafx.h"
#include "FESurfacePairConstraint.h"

REGISTER_SUPER_CLASS(FESurfacePairConstraint, FESURFACEPAIRINTERACTION_ID);

//-----------------------------------------------------------------------------
FESurfacePairConstraint::FESurfacePairConstraint(FEModel* pfem) : FEModelComponent(pfem)
{
}
