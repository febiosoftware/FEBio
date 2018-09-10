#pragma once
#include "fecore_api.h"
#include "DumpStream.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FECORE_API FEElementSet
{
public:
	//! constructor
	FEElementSet(FEMesh* pm);

	void create(int n);

	int size() { return (int)m_Elem.size(); }

	int& operator [] (int i) { return m_Elem[i]; }

	void SetName(const std::string& name);
	const std::string& GetName() const;

	void Serialize(DumpStream& ar);

protected:
	std::string			m_name;		//!< name of element set
	FEMesh*				m_mesh;		//!< pointer to parent mesh
	std::vector<int>	m_Elem;			//!< list of elements
};
