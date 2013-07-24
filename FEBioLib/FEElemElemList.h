#pragma once

#include "FECore/FEMesh.h"
#include <vector>

//-----------------------------------------------------------------------------
//! This class finds for each element the neighbouring elements
//!
class FEElemElemList
{
public:
	FEElemElemList(void);
	~FEElemElemList(void);

	void Create(FEMesh* pmesh);

	FEElement* Neighbor(int n, int j) { return m_pel[ m_ref[n] + j]; }

protected:
	void Init();

protected:
	std::vector<int>	m_ref;		// start index into pel array
	std::vector<FEElement*>	m_pel;	// list of all neighbouring elements (or 0 if no neighbor)
	FEMesh*	m_pmesh;
};
