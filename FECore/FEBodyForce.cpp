#include "stdafx.h"
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
FEBodyForce::FEBodyForce(FEModel* pfem) : m_pfem(pfem)
{
}

//-----------------------------------------------------------------------------
void FEBodyForce::Serialize(DumpFile& ar)
{
	// TODO: Serialize parameters
}
