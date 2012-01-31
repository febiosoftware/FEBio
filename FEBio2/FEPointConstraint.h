#pragma once
#include "FECore/FEConstraint.h"
#include "FECore/FEElement.h"

//-----------------------------------------------------------------------------
//! This class implements a point constraint. That is, it forces a node of a 
//! mesh in the same relative position with respect to the element in which the 
//! node is located.
class FEPointConstraint : public FEConstraint
{
public:
	//! constructor
	FEPointConstraint(FEModel* pfem);

	//! initialize data
	void Init();

	//! Calculate the constraint force
	void Residual(FENLSolver* psolver, vector<double>& R);

	//! calculate the constraint stiffness
	void StiffnessMatrix(FENLSolver* psolver);

	//! augmentations (TODO: implement this)
	bool Augment(int naug) { return true; }

public:
	FEModel*			m_pfem;
	int					m_node;		//!< node to which the constraint is applied
	FESolidElement*		m_pel;		//!< element in which the node is located.
	double				m_rs[3];	//!< natural coordinates in element m_pel
	double				m_eps;		//!< penalty parameter
};
