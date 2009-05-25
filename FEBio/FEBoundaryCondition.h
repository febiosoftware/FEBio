#pragma once

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

class FEBoundaryCondition
{
public:
	FEBoundaryCondition();
	virtual ~FEBoundaryCondition(){}

	bool IsActive() { return m_bactive; }

	void Activate() { m_bactive = true; }
	void Deactivate() { m_bactive = false; }

	int GetID() { return m_nID; }
	void SetID(int nid);

protected:
	bool	m_bactive;	//!< flag indicating whether the BC is active during the current step

	int		m_nID;

	static int	m_ncount;
};
