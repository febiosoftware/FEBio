#pragma once
#include "DumpFile.h"
#include "FEParameterList.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declaration of the model class.
class FEModel;

//-----------------------------------------------------------------------------
//! The FEObject class defines a physical object. An object can be for instance
//! a rigid body or a deformable object. Objects can be connected. For example,
//! a rigid body can be tied to a deformable object.

// NOTE: This is currently only used as a method to abstract the rigid body concept.

class FEObject : public FEParamContainer
{
public:
	//! constructor
	FEObject(FEModel* pfem) : m_fem(*pfem) {}

	//! destructor
	virtual ~FEObject(){}

	// object serialization
	virtual void Serialize(DumpFile& ar) = 0;

	//! shallow copy
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	//! reset object data
	virtual void Reset() = 0;

	//! update solution
	virtual void Update(std::vector<double>& Ui, std::vector<double>& ui) = 0;
	
	//! get the material ID
	virtual int GetMaterialID() { assert(false); return -1; }

protected:
	FEModel&	m_fem;	//!< Pointer to FE model
};
