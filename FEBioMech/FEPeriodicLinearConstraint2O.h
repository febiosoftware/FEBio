#pragma once
#include <FECore/FEMesh.h>
#include <vector>
using namespace std;

class FEModel;

class FEPeriodicLinearConstraint2O
{
	class NodeSetSet
	{
	public:
		NodeSetSet();
		NodeSetSet(const NodeSetSet& nss);
		void operator = (const NodeSetSet& nns);

	public:
		FENodeSet			master;
		vector<FENodeSet>	slave;
	};

public:
	FEPeriodicLinearConstraint2O();
	~FEPeriodicLinearConstraint2O();

	void AddNodeSetPair(const FENodeSet& ms, const FENodeSet& ss, bool push_back = true);
	void AddNodeSetSet(FENodeSet** set, int ncount, bool push_back = true);

	bool GenerateConstraints(FEModel* fem);

	void ExcludeNodes(const FENodeSet& ps) { m_exclude = ps; }

private:
	vector<NodeSetSet>	m_set;	// list of node set pairs
	FENodeSet	m_exclude;		// nodes to exclude
};
