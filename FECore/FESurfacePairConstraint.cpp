#include "stdafx.h"
#include "FESurfacePairConstraint.h"

REGISTER_SUPER_CLASS(FESurfacePairConstraint, FESURFACEPAIRINTERACTION_ID);

//-----------------------------------------------------------------------------
FESurfacePairConstraint::FESurfacePairConstraint(FEModel* pfem) : FEModelComponent(FESURFACEPAIRINTERACTION_ID, pfem)
{
}
