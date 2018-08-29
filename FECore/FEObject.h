#pragma once
#include "FEParameterList.h"
#include "DumpStream.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declaration of the model class.
class FEModel;

//-----------------------------------------------------------------------------
//! The FEObject class defines a physical object. An object can be for instance
//! a rigid body or a deformable object. Objects can be connected. For example,
//! a rigid body can be tied to a deformable object.

// NOTE: This is currently only used as a method to abstract the rigid body concept.

class FECORE_API FEObject : public FEParamContainer
{
public:
	//! constructor
	FEObject(FEModel* pfem) : m_fem(*pfem) {}

	//! destructor
	virtual ~FEObject(){}

	// object serialization
	virtual void Serialize(DumpStream& ar) = 0;

	//! initialize object
	virtual void Init() = 0;

	//! reset object data
	virtual void Reset() = 0;

	//! get the material ID
	virtual int GetMaterialID() { assert(false); return -1; }

protected:
	FEModel&	m_fem;	//!< Pointer to FE model
};
