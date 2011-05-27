#pragma once
#include "FEMaterialPoint.h"
#include "FEMesh.h"
#include "DumpFile.h"
#include "FEParameterList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//! Derived classes need to implement the force and stiffness functions.
//
class FEBodyForce : public FEParamContainer
{
public:
	//! constructor
	FEBodyForce(FEModel* pfem);

	//! calculate the body force at a material point
	virtual vec3d force(FEMaterialPoint& pt) = 0;

	//! calculate constribution to stiffness matrix
	virtual mat3ds stiffness(FEMaterialPoint& pt) = 0;

	//! serialize data to archive
	virtual void Serialize(DumpFile& ar);

	//! initialization
	virtual bool Init() { return true; }

	//! update
	virtual void Update(){}

public:
	FEModel*	m_pfem;	//!< point to model class
	double	s[3];		// force scale factor
	int		lc[3];		// load curve number
};
