#pragma once
#include "FEModelComponent.h"

namespace FECore {

//-----------------------------------------------------------------------------
//! This class is the base class of all boundary conditions

//! Specific boundary conditions can be defined be inheriting from this class.
class FEBoundaryCondition : public FEModelComponent
{
public:
	//! constructor
	FEBoundaryCondition(SUPER_CLASS_ID sid, FEModel* pfem);

	//! desctructor
	virtual ~FEBoundaryCondition(){}

	//! Get the BC ID
	int GetID() { return m_nID; }

	//! Set the BC ID
	void SetID(int nid);

protected:
	int		m_nID;		//!< unique ID for this BC.

	static int	m_ncount;
};

} // namespace FECore
