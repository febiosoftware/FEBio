#pragma once

#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface used for 
//! tied contact

//!	this class is used in contact analyses to describe a contacting
//! surface in a tied contact interface.

class FETiedContactSurface : public FESurface
{
public:
	//! constructor
	FETiedContactSurface(FEMesh* pm=0) : FESurface(pm) {}

	//! Initializes data structures
	void Init();

	//! shallow copy
	void ShallowCopy(FETiedContactSurface& s)
	{
		Lm  = s.Lm;
		gap = s.gap;
	}

	//! Update the surface data
	void Update();

public:
	vector<vec3d>		gap;	//!< gap function at nodes
	vector<FEElement*>	pme;	//!< master element a slave node penetrates
	vector<vec2d>		rs;		//!< natural coordinates of slave projection on master element
	vector<vec3d>		Lm;		//!< Lagrange multipliers
};
