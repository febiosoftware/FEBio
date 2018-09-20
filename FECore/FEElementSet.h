#pragma once
#include "fecore_api.h"
#include "FEItemList.h"
#include "DumpStream.h"
#include "FEElement.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FECORE_API FEElementSet : public FEItemList
{
public:
	//! constructor
	FEElementSet(FEMesh* pm);

	// Create the element set
	void Create(const std::vector<int>& elemList);

	// Return number of elements in the set
	int Elements() const { return (int)m_Elem.size(); }

	int operator [] (int i) const { return m_Elem[i]; }

	void Serialize(DumpStream& ar);

	FEMesh* GetMesh() { return m_mesh; }
	const FEMesh* GetMesh() const { return m_mesh; }

	// return the local index of an element into the element set
	// returns -1 if the element is not part of element set
	int GetLocalIndex(const FEElement& el) const;

protected:
	FEMesh*				m_mesh;		//!< pointer to parent mesh
	std::vector<int>	m_Elem;		//!< list of elements' global ID

	// used for fast lookup in GetLocalIndex
	std::vector<int>	m_LUT;
	int					m_minID, m_maxID;
};

inline int FEElementSet::GetLocalIndex(const FEElement& el) const
{
	return m_LUT[el.GetID() - m_minID];
}
