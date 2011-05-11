#pragma once
#include "FECore/DumpFile.h"
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
//! A degree of freedom structure
class DOF
{
public:
	DOF() { node = bc = neq = -1; }
public:
	int	node;	// the node to which this dof belongs to
	int	bc;		// the degree of freedom
	int	neq;	// the equation number (or -1 if none)
};

//-----------------------------------------------------------------------------
//! linear constraint
class FELinearConstraint
{
public:
	class SlaveDOF : public DOF
	{
	public:
		SlaveDOF() : val(0){}
		double	val;	// coefficient value
	};

public:
	FELinearConstraint(){}
	FELinearConstraint(const FELinearConstraint& LC);

	double FindDOF(int n);

	void Serialize(DumpFile& ar);

public:
	DOF				master;	// master degree of freedom
	list<SlaveDOF>	slave;	// list of slave nodes
};
