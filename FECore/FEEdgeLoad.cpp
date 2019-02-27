#include "stdafx.h"
#include "FEEdgeLoad.h"
#include "FEEdge.h"

REGISTER_SUPER_CLASS(FEEdgeLoad, FEEDGELOAD_ID);

//-----------------------------------------------------------------------------
FEEdgeLoad::FEEdgeLoad(FEModel* pfem) : FEBoundaryCondition(pfem)
{
	m_pedge = 0;
}

FEEdgeLoad::~FEEdgeLoad(void)
{
}


void FEEdgeLoad::Serialize(DumpStream& ar)
{
	FEEdge& e = Edge();
	e.Serialize(ar);
}
