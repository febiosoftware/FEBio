#include "stdafx.h"
#include "FESurfaceLoad.h"

FESurfaceLoad::FESurfaceLoad(FEModel* pfem) : FEBoundaryCondition(FESURFACELOAD_ID), m_pfem(pfem)
{
	m_psurf = 0;
}

FESurfaceLoad::~FESurfaceLoad(void)
{

}
