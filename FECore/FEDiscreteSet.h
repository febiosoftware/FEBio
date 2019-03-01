#pragma once
#include "fecore_api.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
// Forward declarations
class FEMesh;
class DumpStream;

//-----------------------------------------------------------------------------
//! Defines a discrete element set (i.e. node-pairs)
class FECORE_API FEDiscreteSet
{
public:
	struct NodePair
	{
		int	n0, n1;

		void Serialize(DumpStream& ar);
	};

public:
	FEDiscreteSet(FEMesh* pm);
	void create(int n);
	int size() const { return (int)m_pair.size(); }

	void add(int n0, int n1);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	const NodePair& Element(int i) const { return m_pair[i]; }

	void Serialize(DumpStream& ar);

private:
	FEMesh*					m_pmesh;
	std::vector<NodePair>	m_pair;		//!< list of discrete elements
	std::string				m_name;
};
