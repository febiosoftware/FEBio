#include "stdafx.h"
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
FEBodyForce::FEBodyForce(FEModel* pfem) : m_pfem(pfem)
{
}

//-----------------------------------------------------------------------------
//! Serialize body force
//! \todo serialize parameters
void FEBodyForce::Serialize(DumpFile& ar)
{
	
}
