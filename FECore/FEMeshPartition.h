/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEElement.h"
#include "fecore_enum.h"
#include "FESolver.h"
#include "FEGlobalVector.h"
#include "FETimeInfo.h"
#include <functional>

//-----------------------------------------------------------------------------
// forward declaration of classes
class FEModel;
class FENode;
class FEMesh;
class FEDataExport;
class FEGlobalMatrix;
class FEElementSet;

//-----------------------------------------------------------------------------
//! This class describes a mesh partition, that is, a group of elements that represent
//! a part of the mesh.
class FECORE_API FEMeshPartition : public FECoreBase
{
public:
	//! constructor
	FEMeshPartition(int nclass, FEModel* fem);

	//! virtual destructor
	virtual ~FEMeshPartition();

	//! return domain class
	int Class() { return m_nclass; }

	//! set the mesh of this domain
	void SetMesh(FEMesh* pm) { m_pMesh = pm; }

	//! get the mesh of this domain
	FEMesh* GetMesh() { return m_pMesh; }
	const FEMesh* GetMesh() const { return m_pMesh; }

	//! find the element with a specific ID
	FEElement* FindElementFromID(int nid);

	//! serialization
	void Serialize(DumpStream& ar) override;

public:
	//! return number of nodes
	int Nodes() const { return (int)m_Node.size(); }

	//! return a specific node
	FENode& Node(int i);
	const FENode& Node(int i) const;

	//! return the global node index from a local index
	int NodeIndex(int i) const { return m_Node[i]; }

public: // interface for derived classes

	//! return number of elements
	virtual int Elements() const = 0;

	//! return a reference to an element \todo this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;
	virtual const FEElement& ElementRef(int i) const = 0;

public: // optional functions to overload

	//! reset the domain
	virtual void Reset() {}

	//! create a copy of this domain
	virtual void CopyFrom(FEMeshPartition* pd);

	//! initialize domain
	//! one-time initialization, called during model initialization
	bool Init() override;

	//! This function is called at the start of a solution step.
	//! Domain classes can use this to update time dependant quantities
	//! \todo replace this by a version of Update that takes a flag that indicates
	//! whether the update is final or not
	virtual void PreSolveUpdate(const FETimeInfo& timeInfo) {}

	//! Update domain data.
	//! This is called when the model state needs to be updated (i.e. at the end of each Newton iteration)
	virtual void Update(const FETimeInfo& tp) {}

public:
	//! Initialize material points in the domain (optional)
	virtual void InitMaterialPoints() {}

	// Loop over all material points
	void ForEachMaterialPoint(std::function<void(FEMaterialPoint& mp)> f);

	// Loop over all elements
	void ForEachElement(std::function<void(FEElement& el)> f);

public:
	// This is an experimental feature.
	// The idea is to let the class define what data it wants to export
	// The hope is to eliminate the need for special plot and log classes
	// and to automate the I/O completely.
	void AddDataExport(FEDataExport* pd);
	int DataExports() const { return (int)m_Data.size(); }
	FEDataExport* GetDataExport(int i) { return m_Data[i]; }

public:
	bool IsActive() const { return m_bactive; }
	void SetActive(bool b) { m_bactive = b; }

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of
	vector<int>	m_Node;		//!< list of nodes in this domain

protected:
	int	m_nclass;			//!< domain class

	bool	m_bactive;

private:
	vector<FEDataExport*>	m_Data;	//!< list of data export classes
};
