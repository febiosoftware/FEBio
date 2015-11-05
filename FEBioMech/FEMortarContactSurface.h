#pragma once
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class represents a surface used by the mortar contact interface.
class FEMortarContactSurface : public FEContactSurface
{
public:
	FEMortarContactSurface(FEMesh* pm = 0);

	//! Initializes data structures
	bool Init();

	//! update nodal areas
	void UpdateNodalAreas();

public:
	vector<double>	m_A;		//!< nodal areas
	vector<vec3d>	m_gap;		//!< nodal gaps
};
