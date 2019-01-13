#pragma once
#include <FECore/FEMesh.h>
#include "febiomech_api.h"
#include <vector>

class FEModel;

class FEBIOMECH_API FEPeriodicLinearConstraint2O
{
	class NodeSetSet
	{
	public:
		NodeSetSet();
		NodeSetSet(const NodeSetSet& nss);
		void operator = (const NodeSetSet& nns);

	public:
		FENodeSet	master;
		FENodeSet	slave;
	};

public:
	FEPeriodicLinearConstraint2O();
	~FEPeriodicLinearConstraint2O();

	void AddNodeSetPair(const FENodeSet& ms, const FENodeSet& ss, bool push_back = true);

	bool GenerateConstraints(FEModel* fem);

private:
	int closestNode(FEMesh& mesh, const FENodeSet& set, const vec3d& r);
	void addLinearConstraint(FEModel& fem, int master, int slave);

private:
	std::vector<NodeSetSet>	m_set;	// list of node set pairs
};
