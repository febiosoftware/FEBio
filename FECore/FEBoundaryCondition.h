#pragma once
#include "FECoreBase.h"

namespace FECore {

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

//! Specific boundary conditions can be defined be inheriting from this class.
class FEBoundaryCondition : public FECoreBase
{
public:
	//! constructor
	FEBoundaryCondition(SUPER_CLASS_ID sid);

	//! desctructor
	virtual ~FEBoundaryCondition(){}

	//! return the active status of this BC
	bool IsActive() { return m_bactive; }

	//! activate BC
	void Activate() { m_bactive = true; }

	//! Deactive BC
	void Deactivate() { m_bactive = false; }

	//! Get the BC ID
	int GetID() { return m_nID; }

	//! Set the BC ID
	void SetID(int nid);

protected:
	bool	m_bactive;	//!< flag indicating whether the BC is active during the current step
	int		m_nID;		//!< unique ID for this BC.

	static int	m_ncount;
};

} // namespace FECore
