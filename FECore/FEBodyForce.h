#pragma once
#include "FEMaterialPoint.h"
#include "FEMesh.h"
#include "DumpFile.h"
#include "FEParameterList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
// Base class for body-loads (TODO: work in progress)
class FEBodyLoad : public FEParamContainer
{
public:
	FEBodyLoad(){}
	virtual ~FEBodyLoad(){}

	//! initialization
	virtual bool Init() { return true; }
};

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//! Derived classes need to implement the force and stiffness functions.
//
class FEBodyForce : public FEBodyLoad
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

	//! update
	virtual void Update(){}

public:
	FEModel*	m_pfem;	//!< point to model class
};
