#include "stdafx.h"
#include "FENLConstraint.h"


//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint(FEModel* pfem) : m_pfem(pfem) { m_bactive = true; }

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint(){}

//-----------------------------------------------------------------------------
bool FENLConstraint::IsActive()
{ 
	return m_bactive; 
}

//-----------------------------------------------------------------------------
void FENLConstraint::Activate()
{ 
	m_bactive = true; 
}

//-----------------------------------------------------------------------------
void FENLConstraint::Deactivate()
{ 
	m_bactive = false; 
}
