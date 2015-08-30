#pragma once

#include "FEElement.h"
#include "DumpFile.h"
#include "DumpStream.h"
#include "FE_enum.h"
#include "FESolver.h"
#include "FEGlobalVector.h"

//-----------------------------------------------------------------------------
// forward declaration of classes
class FEModel;
class FENode;
class FEMesh;
class FEMaterial;

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//! and material. 
class FEDomain : public FECoreBase
{
public:
	enum { MAX_DOMAIN_NAME = 64 };

public:
	//! constructor
	FEDomain(int nclass, FEMesh* pm);

	//! virtual destructor
	virtual ~FEDomain();

	//! return domain class
	int Class() { return m_nclass; }

	//! set the mesh of this domain
	void SetMesh(FEMesh* pm) { m_pMesh = pm; }

	//! get the mesh of this domain
	FEMesh* GetMesh() { return m_pMesh; }

	//! find the element with a specific ID
	FEElement* FindElementFromID(int nid);

public:
	// Set the domain name
	void SetName(const char* szname);

	// get the domain name
	const char* GetName();

public:
	//! get the material of this domain
	//! \todo Delete this.
	virtual FEMaterial* GetMaterial() { return 0; }

	// assign a material to this domain
	virtual void SetMaterial(FEMaterial* pm) {}

	//! set the material ID of all elements
	void SetMatID(int mid);

public: // interface for derived classes
	
	//! create a domain of n elements
	virtual void create(int n) = 0;

	//! return number of nodes
	virtual int Nodes() = 0;

	//! return a specific node
	virtual FENode& Node(int i) = 0;

	//! return number of elements
	virtual int Elements() = 0;

	//! return a reference to an element \todo this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;

	//! Unpack the LM data for an element of this domain
	virtual void UnpackLM(FEElement& el, vector<int>& lm) = 0;

public: // optional functions to overload

	//! reset the domain
	virtual void Reset() {}

	//! create a copy of this domain
	virtual FEDomain* Copy() { return 0; }

	//! stream domain data
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) {}

	//! initialize domain
	virtual bool Initialize(FEModel& fem) { return true; }

	//! Initialize elements of domain
	virtual void InitElements() {}

	//! Initialize material point data for the elements
	void InitMaterialPointData();

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of

protected:
	int	m_nclass;			//!< domain class

private:
	char	m_szname[MAX_DOMAIN_NAME];	//!< domain name
};
