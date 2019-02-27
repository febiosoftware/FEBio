#include "stdafx.h"
#include "FEBoundaryCondition.h"

REGISTER_SUPER_CLASS(FEBoundaryCondition, FEBC_ID);

FEBoundaryCondition::FEBoundaryCondition(FEModel* pfem) : FEModelComponent(pfem)
{
}
