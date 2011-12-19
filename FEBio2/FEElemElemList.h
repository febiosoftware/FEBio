#pragma once

#include "NumCore/vector.h"
#include "FECore/FEMesh.h"

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
	vector<int>	m_ref;		// start index into pel array
	vector<FEElement*>	m_pel;	// list of all neighbouring elements (or 0 if no neighbor)
	FEMesh*	m_pmesh;
};
