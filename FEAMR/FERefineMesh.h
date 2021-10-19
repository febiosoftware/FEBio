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
#include <FECore/FEMeshAdaptor.h>
class FEModel;

class FEMeshTopo;
class FEBoundaryCondition;
class FESurfaceLoad;
class FESurfacePairConstraint;
class FESurface;
class FENodeSet;
class FEDomainMap;

//-----------------------------------------------------------------------------
// Base class for mesh refinement algorithms
class FERefineMesh : public FEMeshAdaptor
{
protected:
	// Supported transfer methods for mapping data between meshes
	enum TransferMethod {
		TRANSFER_SHAPE,
		TRANSFER_MLQ
	};
	
public:
	FERefineMesh(FEModel* fem);
	~FERefineMesh();

	// Apply mesh refinement
	bool Apply(int iteration) override;

protected:
	// Derived classes need to override this function
	// The return value should be true if the mesh was refined
	// or false otherwise. 
	virtual bool RefineMesh() = 0;

protected:
	bool BuildMeshTopo();
	void CopyMesh();

	bool BuildMapData();
	void TransferMapData();
	void ClearMapData();

	bool BuildDomainMapData();
	bool BuildDomainMapData(FEDomain& dom, int domIndex);
	void TransferDomainMapData();

	bool BuildUserMapData();
	void TransferUserMapData();

protected:
	FEMeshTopo*	m_topo;		//!< mesh topo structure

	int		m_maxiter;		// max nr of iterations per time step
	int		m_maxelem;		// max nr of elements

	int		m_transferMethod;		//!< method for transferring data between meshes
	bool	m_bmap_data;			//!< map data flag
	int		m_nnc;					//!< nearest-neighbor-count for MLQ transfer method
	int		m_nsdim;				//!< nearest-neighbor search dimension (2 or 3)

	FEMesh*	m_meshCopy;		//!< copy of "old" mesh, before refinement

	std::vector< std::vector<FEDomainMap*> >	m_domainMapList;	// list of nodal data for each domain
	std::vector< FEDomainMap* >	m_userDataList;						// list of nodal data for user-defined mesh data

	DECLARE_FECORE_CLASS();
};
