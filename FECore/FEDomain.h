#pragma once

#include "FEElement.h"
#include "FE_enum.h"
#include "FESolver.h"
#include "FEGlobalVector.h"
#include "FETypes.h"

//-----------------------------------------------------------------------------
// forward declaration of classes
class FEModel;
class FENode;
class FEMesh;
class FEMaterial;
class FEDataExport;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//! and material. 
//! Creating a domain requires the following steps:
//! 1. Create an instance of domain class 
//! 2. Call its Create(int nsize, int elemType) member with the number of elements and element type
//! 3. For each element of the domain, call the SetMatID, and set the node numbers
//! 4. Attach the domain to the mesh FEMesh::AddDomain
//! 5. Call the CreateMaterialPointData
//! \todo I'd like to simplify these steps:
//! A. Creating a domain should automatically add it to the mesh (eliminate step 4)
//! B. Also, I want to eliminate the call to SetMatID. 
//! C. perhaps the call to CreateMaterialPointData can be done in Create? 
//! 
class FECOREDLL_EXPORT FEDomain : public FECoreBase
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

	//! serialization
	void Serialize(DumpStream& ar);

public:
	//! get the material of this domain
	//! \todo Delete this.
	virtual FEMaterial* GetMaterial() { return 0; }

	// assign a material to this domain
	virtual void SetMaterial(FEMaterial* pm) {}

	//! set the material ID of all elements
	void SetMatID(int mid);

	//! return number of nodes
	int Nodes() const { return (int) m_Node.size(); }

	//! return a specific node
	FENode& Node(int i);
	const FENode& Node(int i) const;

	//! return the global node index from a local index
	int NodeIndex(int i) const { return m_Node[i]; }

public: // interface for derived classes
	
	//! create a domain of n elements of type
	virtual void Create(int nelems, int elemType = -1) = 0;

	//! return number of elements
	virtual int Elements() const = 0;

	//! return a reference to an element \todo this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;

	//! Unpack the LM data for an element of this domain
	virtual void UnpackLM(FEElement& el, vector<int>& lm);

public: // optional functions to overload

	//! reset the domain
	virtual void Reset() {}

	//! create a copy of this domain
	virtual void CopyFrom(FEDomain* pd);

	//! initialize domain
	//! one-time initialization, called during model initialization
	bool Init() override;

	//! This function is called at the start of a solution step.
	//! Domain classes can use this to update time dependant quantities
	//! \todo replace this by a version of Update that takes a flag that indicates
	//! whether the update is final or not
	virtual void PreSolveUpdate(const FETimeInfo& timeInfo) {}

	//! Activate the domain
	virtual void Activate();

	//! Update domain data.
	//! This is called when the model state needs to be updated (i.e. at the end of each Newton iteration)
	virtual void Update(const FETimeInfo& tp) {}

	//! build the matrix profile
	virtual void BuildMatrixProfile(FEGlobalMatrix& M);

public:
	void SetDOFList(vector<int>& dof);

	int GetDOFCount() const { return (int) m_dof.size(); }

	const vector<int>& GetDOFList() const { return m_dof; }

	//! Allocate material point data for the elements
	//! This is called after elements get read in from the input file.
	//! And must be called before material point data can be accessed.
	//! \todo Perhaps I can make this part of the "creation" routine
	void CreateMaterialPointData();

public:
	// This is an experimental feature.
	// The idea is to let the class define what data it wants to export
	// The hope is to eliminate the need for special plot and log classes
	// and to automate the I/O completely.
	void AddDataExport(FEDataExport* pd);
	int DataExports() const { return (int) m_Data.size(); }
	FEDataExport* GetDataExport(int i) { return m_Data[i]; }

public:
	bool IsActive() const { return m_bactive; }
	void SetActive(bool b) { m_bactive = b; }

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of
	vector<int>	m_Node;		//!< list of nodes in this domain
	vector<int>	m_dof;		//!< list of active degrees of freedom for this domain

protected:
	int	m_nclass;			//!< domain class

	bool	m_bactive;

private:
	vector<FEDataExport*>	m_Data;	//!< list of data export classes
};
