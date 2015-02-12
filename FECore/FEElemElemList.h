#pragma once

#include "FEMesh.h"
#include <vector>

//-----------------------------------------------------------------------------
//! This class finds for each element the neighbouring elements
//!
class FEElemElemList
{
public:
	//! constructor
	FEElemElemList(void);

	//! destructor
	~FEElemElemList(void);

	//! create the element-element list
	void Create(FEMesh* pmesh);

	//! create the element-element list for a surface
	void Create(FESurface* psurf);

	//! Find the j-th neighbor element of element n
	FEElement* Neighbor(int n, int j) { return m_pel[ m_ref[n] + j]; }

protected:
	//! Initialization
	void Init();

protected:
	std::vector<int>	m_ref;		//!< start index into pel array
	std::vector<FEElement*>	m_pel;	//!< list of all neighbouring elements (or 0 if no neighbor)
	FEMesh*	m_pmesh;				//!< pointer to mesh that created this list
};
