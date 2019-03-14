#pragma once
#include "fecore_api.h"
#include "FEItemList.h"
#include "FEElement.h"
#include "FEDomainList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;
class DumpStream;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FECORE_API FEElementSet : public FEItemList
{
public:
	//! constructor
	FEElementSet(FEMesh* pm);

	// Create the element set
	void Create(const std::vector<int>& elemList);

	// Create the element set from a domain
	void Create(FEDomain* dom);

	// Create the element set from a domain
	void Create(FEDomainList& dom);

	// Return number of elements in the set
	int Elements() const { return (int)m_Elem.size(); }

	int operator [] (int i) const { return m_Elem[i]; }

	FEMesh* GetMesh() { return m_mesh; }
	const FEMesh* GetMesh() const { return m_mesh; }

	// return the local index of an element into the element set
	// returns -1 if the element is not part of element set
	int GetLocalIndex(const FEElement& el) const;

	// Get the element ID list
	const std::vector<int>& GetElementIDList() const { return m_Elem; }

	// get the domain list that generated the element set
	FEDomainList& GetDomainList() { return m_dom; }

	// Get an element
	FEElement& Element(int i);

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FEElementSet* p);
	static FEElementSet* LoadClass(DumpStream& ar, FEElementSet* p);

private:
	// Build the lookup table
	void BuildLUT();

protected:
	FEMesh*				m_mesh;		//!< pointer to parent mesh
	std::vector<int>	m_Elem;		//!< list of elements' global ID

	FEDomainList		m_dom;	//!< domain list that generated the element set

	// used for fast lookup in GetLocalIndex
	std::vector<int>	m_LUT;
	int					m_minID, m_maxID;
};

inline int FEElementSet::GetLocalIndex(const FEElement& el) const
{
	return m_LUT[el.GetID() - m_minID];
}
