#pragma once
#include "FEModelComponent.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! linear constraint
class FELinearConstraint : public FEModelComponent
{
public:
	class DOF
	{
	public:
		DOF()
		{
			node = dof = -1;
			val = 0.0;
		}

	public:
		int		node;	// node number
		int		dof;	// degree of freedom
		double	val;	// coefficient value (ignored for master)
	};

public:
	// constructors
	FELinearConstraint(FEModel* pfem);
	FELinearConstraint(const FELinearConstraint& LC);

	// copy data
	void CopyFrom(const FELinearConstraint& LC);

	// serialization
	void Serialize(DumpStream& ar);

	// initialize the linear constraint
	bool Init();

	// make the constraint active
	void Activate();
	void Deactivate();

public:
	DOF			master;	// master degree of freedom
	vector<DOF>	slave;	// list of slave nodes
};
