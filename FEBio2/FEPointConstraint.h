#pragma once
#include "FECore/FEMesh.h"

class FEM;

//-----------------------------------------------------------------------------
//! This class implements a point constraint. That is, it forces a node of a 
//! mesh in the same relative position with respect to the element in which the 
//! node is located.
class FEPointConstraint
{
public:
	FEPointConstraint(FEM* pfem);
	~FEPointConstraint(){}

	//! initialize data
	void Init();

	//! Calculate the constraint force
	void Residual(vector<double>& R);

	//! calculate the constraint stiffness
	void Stiffness();

public:
	FEM*				m_pfem;
	int					m_node;		//!< node to which the constraint is applied
	FESolidElement*		m_pel;		//!< element in which the node is located.
	double				m_rs[3];	//!< natural coordinates in element m_pel
	double				m_eps;		//!< penalty parameter
};
