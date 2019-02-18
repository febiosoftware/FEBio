#include "stdafx.h"
#include "FEPrescribedBC.h"
#include "FEFacetSet.h"

//-----------------------------------------------------------------------------
FEPrescribedBC::FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
}

//-----------------------------------------------------------------------------
// assign a surface to the BC
void FEPrescribedBC::AddNodes(const FEFacetSet& surf)
{
	FENodeSet nset = surf.GetNodeSet();
	AddNodes(nset);
}
