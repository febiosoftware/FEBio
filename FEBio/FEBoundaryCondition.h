#pragma once

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

class FEBoundaryCondition
{
public:
	FEBoundaryCondition() { m_bactive = true; }

	bool IsActive() { return m_bactive; }

	void Activate() { m_bactive = true; }
	void Deactivate() { m_bactive = false; }

protected:
	bool	m_bactive;	//!< flag indicating whether the BC is active during the current step
};
