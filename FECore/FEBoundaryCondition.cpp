#include "stdafx.h"
#include "FEBoundaryCondition.h"

REGISTER_SUPER_CLASS(FEBoundaryCondition, FEBC_ID);

FEBoundaryCondition::FEBoundaryCondition(SUPER_CLASS_ID sid, FEModel* pfem) : FEModelComponent(sid, pfem)
{
}
