#include "stdafx.h"
#include "FESurfaceLoad.h"

FESurfaceLoad::FESurfaceLoad(FEModel* pfem) : FEBoundaryCondition(FESURFACELOAD_ID, pfem)
{
	m_psurf = 0;
}

FESurfaceLoad::~FESurfaceLoad(void)
{

}
