#include "stdafx.h"
#include "FEBodyForce.h"


//-----------------------------------------------------------------------------
FEBodyForce::FEBodyForce(FEModel* pfem) : FEBodyLoad(pfem)
{
}

//-----------------------------------------------------------------------------
//! Serialize body force
//! \todo serialize parameters
void FEBodyForce::Serialize(DumpFile& ar)
{
	return;
}
