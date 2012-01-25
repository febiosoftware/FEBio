#pragma once

#include "FEElement.h"
#include "DumpFile.h"
#include "FE_enum.h"
#include "FENLSolver.h"

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
class FEDomain
{
public:
	//! constructor
	FEDomain(int ntype, FEMesh* pm, FEMaterial* pmat) { m_pMesh = pm; m_ntype = ntype; m_pMat = pmat; }

	//! virtual destructor
	virtual ~FEDomain() {}

	//! return domain type
	int Type() { return m_ntype; }

	//! set the mesh of this domain
	void SetMesh(FEMesh* pm) { m_pMesh = pm; }

	//! get the mesh of this domain
	FEMesh* GetMesh() { return m_pMesh; }

	//! find the element with a specific ID
	FEElement* FindElementFromID(int nid);

	//! set the material of this domain
	void SetMaterial(FEMaterial* pmat) { m_pMat = pmat; }

	//! get the material of this domain
	FEMaterial* GetMaterial() { return m_pMat; }

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

	// return a reference to an element (TODO: this is not the preferred interface but I've added it for now)
	virtual FEElement& ElementRef(int i) = 0;

	//! Unpack the LM data for an element of this domain
	virtual void UnpackLM(FEElement& el, vector<int>& lm) = 0;

public: // optional functions to overload

	//! create a clone of this domain (used in running restarts)
	virtual FEDomain* Clone() { assert(false); return 0; }

	//! reset the domain
	virtual void Reset() {}

	//! serialize domain to archive
	virtual void Serialize(DumpFile& ar) {}

	//! initialize domain
	virtual bool Initialize(FEModel& fem) { return true; }

	//! Initialize elements of domain
	virtual void InitElements() {}

	//! Initialize material point data for the elements
	void InitMaterialPointData();

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of
	FEMaterial*	m_pMat;		//!< the material for this domain

protected:
	int	m_ntype;			//! type of domain
};
