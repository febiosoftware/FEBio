#include "stdafx.h"
#include "FEEdgeLoad.h"
#include "FEEdge.h"

//-----------------------------------------------------------------------------
FEEdgeLoad::FEEdgeLoad(FEModel* pfem) : FEBoundaryCondition(FEEDGELOAD_ID, pfem)
{
	m_pedge = 0;
}

FEEdgeLoad::~FEEdgeLoad(void)
{
}
