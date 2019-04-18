/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

#include "stdafx.h"
#include "FEMesh.h"
#include "FEException.h"
#include "FEDiscreteDomain.h"
#include "FETrussDomain.h"
#include "FEShellDomain.h"
#include "FESolidDomain.h"
#include "FEDomain2D.h"
#include "DOFS.h"
#include "FEElemElemList.h"
#include "FEElementList.h"
#include "FESurface.h"
#include "FEDataArray.h"
#include "FEDomainMap.h"
#include "FESurfaceMap.h"
#include "DumpStream.h"
#include <algorithm>

//-----------------------------------------------------------------------------
FEDataArray* CreateDataMap(int mapType)
{
	FEDataArray* map = nullptr;
	switch (mapType)
	{
	case FE_DOMAIN_MAP : map = new FEDomainMap; break;
	case FE_SURFACE_MAP: map = new FESurfaceMap; break;
	default:
		assert(false);
	}

	return map;
}

//=============================================================================
// FEMesh
//-----------------------------------------------------------------------------
FEMesh::FEMesh(FEModel* fem) : m_fem(fem)
{
	m_LUT = 0;
}

//-----------------------------------------------------------------------------
FEMesh::~FEMesh()
{
	Clear();
}

//-----------------------------------------------------------------------------
//! return number of nodes
int FEMesh::Nodes() const 
{ 
	return (int)m_Node.size(); 
}

//-----------------------------------------------------------------------------
void FEMesh::Serialize(DumpStream& ar)
{
	// clear the mesh if we are loading from an archive
	if ((ar.IsShallow() == false) && (ar.IsLoading())) Clear();

	// we don't want to store pointers to all the nodes
	// mostly for efficiency, so we tell the archive not to store the pointers
	ar.LockPointerTable();
	{
		// store the node list
		ar & m_Node;
	}
	ar.UnlockPointerTable();

	// stream domain data
	ar & m_Domain;

	// if this is a shallow archive, we're done
	if (ar.IsShallow()) return;

	// serialize node sets
	ar & m_NodeSet;

	if (ar.IsSaving())
	{
		// write segment sets
		int ssets = SegmentSets();
		ar << ssets;
		for (int i=0; i<ssets; ++i)
		{
			FESegmentSet& sset = SegmentSet(i);
			sset.Serialize(ar);
		}

		// write element sets
		ar << m_ElemSet;

		// write discrete sets
		int dsets = DiscreteSets();
		ar << dsets;
		for (int i=0; i<dsets; ++i)
		{
			FEDiscreteSet& dset = DiscreteSet(i);
			dset.Serialize(ar);
		}

		// write data maps
		int maps = DataArrays();
		ar << maps;
		for (int i = 0; i < maps; ++i)
		{
			FEDataArray* map = GetDataArray(i);
			ar << (int)map->DataMapType();
			ar << DataArrayName(i);
			map->Serialize(ar);
		}
	}
	else
	{
		FEModel* fem = &ar.GetFEModel();

		// read segment sets
		int ssets = 0;
		ar >> ssets;
		for (int i=0; i<ssets; ++i)
		{
			FESegmentSet* sset = new FESegmentSet(fem);
			AddSegmentSet(sset);
			sset->Serialize(ar);
		}

		// read element sets
		ar >> m_ElemSet;

		// read discrete sets
		int dsets = 0;
		ar >> dsets;
		for (int i=0; i<dsets; ++i)
		{
			FEDiscreteSet* dset = new FEDiscreteSet(this);
			AddDiscreteSet(dset);
			dset->Serialize(ar);
		}

		// write data maps
		ClearDataArrays();
		int maps = 0, mapType;
		string mapName;
		ar >> maps;
		for (int i = 0; i < maps; ++i)
		{
			ar >> mapType;

			FEDataArray* map = CreateDataMap(mapType);
			assert(map);

			ar >> mapName;
			map->Serialize(ar);
			AddDataArray(mapName, map);
		}

		UpdateBox();
	}
}

//-----------------------------------------------------------------------------
void FEMesh::SaveClass(DumpStream& ar, FEMesh* p)
{
	// we should never get here
	assert(false);
}

//-----------------------------------------------------------------------------
FEMesh* FEMesh::LoadClass(DumpStream& ar, FEMesh* p)
{
	// we should never get here
	assert(false);
	return nullptr;
}

//-----------------------------------------------------------------------------
//  Allocates storage for mesh data.
//
void FEMesh::CreateNodes(int nodes)
{
	assert(nodes);
	m_Node.resize(nodes);

	// set the default node IDs
	for (int i=0; i<nodes; ++i) Node(i).SetID(i+1);
}

//-----------------------------------------------------------------------------
// Make more room for nodes
void FEMesh::AddNodes(int nodes)
{
	assert(nodes);
	int N0 = (int) m_Node.size();

	// get the ID of the last node
	// (It is assumed that nodes are sorted according their ID
	//  so the last node should have the highest ID)
	int n0 = 1;
	if (N0 > 0) n0 = m_Node[N0-1].GetID() + 1;

	m_Node.resize(N0 + nodes);
	for (int i=0; i<nodes; ++i) m_Node[i+N0].SetID(n0+i);
}

//-----------------------------------------------------------------------------
void FEMesh::SetDOFS(int n)
{
	int NN = Nodes();
	for (int i=0; i<NN; ++i) m_Node[i].SetDOFS(n);
}

//-----------------------------------------------------------------------------
//! Return the total number elements
int FEMesh::Elements() const
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) 
	{
		N += m_Domain[i]->Elements();
	}
	return N;
}

//-----------------------------------------------------------------------------
//! Return the total number of elements of a specific domain type
int FEMesh::Elements(int ndom_type) const
{
	int N = 0;
	for (int i=0; i<(int) m_Domain.size(); ++i) 
	{
		FEDomain& dom = *m_Domain[i];
		if (dom.Class() == ndom_type) N += m_Domain[i]->Elements();
	}
	return N;
}

//-----------------------------------------------------------------------------
//! return reference to a node
FENode& FEMesh::Node(int i) { return m_Node[i]; }
const FENode& FEMesh::Node(int i) const { return m_Node[i]; }

//-----------------------------------------------------------------------------
//  Updates the bounding box of the mesh (using current coordinates)
//
void FEMesh::UpdateBox()
{
	if (Nodes() > 0)
	{
		m_box = FEBoundingBox(Node(0).m_rt);
		for (int i=1; i<Nodes(); ++i)
		{
			m_box.add(Node(i).m_rt);
		}
	}
	else m_box = FEBoundingBox(vec3d(0,0,0));
}

//-----------------------------------------------------------------------------
//  Counts the number of shell elements in the mesh
//
int FEMesh::RemoveIsolatedVertices()
{
	int i, j, k, N = Nodes(), n;

	// create a valence array
	vector<int> val; val.assign(N, 0);

	// count the nodal valences
	for (i=0; i<(int) m_Domain.size(); ++i)
	{
		FEDomain& d = Domain(i);
		for (j=0; j<d.Elements(); ++j)
		{
			FEElement& el = d.ElementRef(j);
			n = el.Nodes();
			for (k=0; k<n; ++k) ++val[el.m_node[k]];
		}
	}

	// See if there are any isolated nodes
	// Exclude them from the analysis
	int ni = 0;
	for (i=0; i<N; ++i)
		if (val[i] == 0)
		{
			++ni;
			FENode& node = Node(i);
			node.SetFlags(FENode::EXCLUDE);
		}

	return ni;
}

//-----------------------------------------------------------------------------
//! Does one-time initialization of the Mesh material point data.
void FEMesh::InitMaterialPoints()
{
    for (int i = 0; i<Domains(); ++i)
    {
        FEDomain& dom = Domain(i);
        dom.InitMaterialPoints();
    }
}

//-----------------------------------------------------------------------------
void FEMesh::Clear()
{
	m_Node.clear();
	for (size_t i=0; i<m_Domain.size (); ++i) delete m_Domain [i];

	// TODO: Surfaces are currently managed by the classes that use them so don't delete them
//	for (size_t i=0; i<m_Surf.size   (); ++i) delete m_Surf   [i];

	for (size_t i=0; i<m_NodeSet.size (); ++i) delete m_NodeSet [i];
	for (size_t i=0; i<m_LineSet.size (); ++i) delete m_LineSet [i];
	for (size_t i=0; i<m_ElemSet.size (); ++i) delete m_ElemSet [i];
	for (size_t i=0; i<m_DiscSet.size (); ++i) delete m_DiscSet [i];

	m_Domain.clear();
	m_Surf.clear();
	m_NodeSet.clear();
	m_LineSet.clear();
	m_ElemSet.clear();
	m_DiscSet.clear();

	m_NEL.Clear();
	if (m_LUT) delete m_LUT; m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Reset the mesh data. Return nodes to their intial position, reset their 
//! attributes and zero all element stresses.

void FEMesh::Reset()
{
	// reset nodal data
	for (int i=0; i<Nodes(); ++i) 
	{
		FENode& node = Node(i);

		node.m_rp = node.m_rt = node.m_r0;
		node.m_vp = vec3d(0,0,0);
		node.m_ap = node.m_at = vec3d(0,0,0);

        node.m_Fr = vec3d(0,0,0);

		// reset ID arrays
		int ndof = (int)node.dofs();
		for (int i=0; i<ndof; ++i) 
		{
			node.set_inactive(i);
			node.set_bc(i, DOF_OPEN);
			node.set(i, 0.0);
		}
	}

	// update the mesh
	UpdateBox();

	// reset domain data
	for (int n=0; n<(int) m_Domain.size(); ++n) m_Domain[n]->Reset();
}

//-----------------------------------------------------------------------------
//! This function calculates the (initial) volume of an element. In some case, the volume
//! may only be approximate.
double FEMesh::ElementVolume(FEElement &el)
{
	double V = 0;
	switch (el.Class())
	{
	case FE_ELEM_SOLID: V = SolidElementVolume(static_cast<FESolidElement&>(el)); break;
	case FE_ELEM_SHELL: V = ShellElementVolume(static_cast<FEShellElement&>(el)); break;
	}
	return V;
}

//-----------------------------------------------------------------------------
double FEMesh::SolidElementVolume(FESolidElement& el)
{
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(el.GetMeshPartition()); assert(dom);
	if (dom)
		return dom->Volume(el);
	else
		return 0.0;
}

//-----------------------------------------------------------------------------
double FEMesh::ShellElementVolume(FEShellElement& el)
{
	FEShellDomain* dom = dynamic_cast<FEShellDomain*>(el.GetMeshPartition()); assert(dom);
	if (dom)
		return dom->Volume(el);
	else
		return 0.0;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by ID

FENodeSet* FEMesh::FindNodeSet(int nid)
{
	for (size_t i=0; i<m_NodeSet.size(); ++i) if (m_NodeSet[i]->GetID() == nid) return m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a nodeset by name

FENodeSet* FEMesh::FindNodeSet(const std::string& name)
{
	for (size_t i=0; i<m_NodeSet.size(); ++i) if (m_NodeSet[i]->GetName() == name) return m_NodeSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a segment set set by name

FESegmentSet* FEMesh::FindSegmentSet(const std::string& name)
{
	for (size_t i=0; i<m_LineSet.size(); ++i) if (m_LineSet[i]->GetName() ==  name) return m_LineSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a surface set set by name

FESurface* FEMesh::FindSurface(const std::string& name)
{
	for (size_t i=0; i<m_Surf.size(); ++i) if (m_Surf[i]->GetName() == name) return m_Surf[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a discrete element set set by name

FEDiscreteSet* FEMesh::FindDiscreteSet(const std::string& name)
{
	for (size_t i=0; i<m_DiscSet.size(); ++i) if (m_DiscSet[i]->GetName() == name) return m_DiscSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
//! Find a element set by name

FEElementSet* FEMesh::FindElementSet(const std::string& name)
{
	for (size_t i=0; i<m_ElemSet.size(); ++i) if (m_ElemSet[i]->GetName() == name) return m_ElemSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::Domains() { return (int)m_Domain.size(); }

//-----------------------------------------------------------------------------
FEDomain& FEMesh::Domain(int n) { return *m_Domain[n]; }

//-----------------------------------------------------------------------------
void FEMesh::AddDomain(FEDomain* pd)
{ 
	int N = (int)m_Domain.size();
	pd->SetID(N);
	m_Domain.push_back(pd); 
	if (m_LUT) delete m_LUT; m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Find a domain

FEDomain* FEMesh::FindDomain(const std::string& name)
{
	for (size_t i = 0; i<m_Domain.size(); ++i) if (m_Domain[i]->GetName() == name) return m_Domain[i];
	return 0;
}

FEDomain* FEMesh::FindDomain(int domId)
{
	for (size_t i = 0; i<m_Domain.size(); ++i) if (m_Domain[i]->GetID() == domId) return m_Domain[i];
	return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::Faces(FEElement& el)
{
	switch (el.Type())
	{
	case FE_HEX8G8:
	case FE_HEX8RI:
	case FE_HEX8G1:
    case FE_HEX20G8:
	case FE_HEX20G27:
	case FE_HEX27G27: return 6;
	case FE_PENTA6G6:
    case FE_PENTA15G8:
    case FE_PENTA15G21:
	case FE_PYRA5G8: return 5;
	case FE_TET4G4:
	case FE_TET10G4:
	case FE_TET10G8:
	case FE_TET10GL11:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET20G15:
	case FE_TET4G1: return 4;
    case FE_SHELL_QUAD4G8:
    case FE_SHELL_QUAD4G12:
    case FE_SHELL_QUAD8G18:
    case FE_SHELL_QUAD8G27:
    case FE_SHELL_TRI3G6:
    case FE_SHELL_TRI3G9:
    case FE_SHELL_TRI6G14:
    case FE_SHELL_TRI6G21: return 1;
	default:
		assert(false);
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! This function returns the face connectivity from a certain element

int FEMesh::GetFace(FEElement& el, int n, int* nf)
{
	int nn = -1;
	int* en = &el.m_node[0];
	switch (el.Type())
	{
	case FE_HEX8G8:
	case FE_HEX8RI:
	case FE_HEX8G1:
		nn = 4;
		switch (n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; break;
		}
		break;
	case FE_PENTA6G6:
		switch(n)
		{
		case 0: nn = 4; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; break;
		case 1: nn = 4; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; break;
		case 2: nn = 4; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; break;
		case 3: nn = 3; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[1]; break;
		case 4: nn = 3; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[5]; break;
		}
		break;
    case FE_PENTA15G8:
    case FE_PENTA15G21:
        switch(n)
        {
            case 0: nn = 8; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; nf[3] = en[3]; nf[4] = en[ 6]; nf[5] = en[13]; nf[6] = en[ 9]; nf[7] = en[12]; break;
            case 1: nn = 8; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[ 7]; nf[5] = en[14]; nf[6] = en[10]; nf[7] = en[13]; break;
            case 2: nn = 8; nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[5]; nf[3] = en[2]; nf[4] = en[12]; nf[5] = en[11]; nf[6] = en[14]; nf[7] = en[ 8]; break;
            case 3: nn = 6; nf[0] = en[0]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[8]; nf[4] = en[ 7]; nf[5] = en[ 6]; break;
            case 4: nn = 6; nf[0] = en[3]; nf[1] = en[4]; nf[2] = en[5]; nf[3] = en[9]; nf[4] = en[10]; nf[5] = en[11]; break;
        }
        break;
	case FE_PYRA5G8:
		switch (n)
		{
			case 0: nn = 3; nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[4]; break;
			case 1: nn = 3; nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[4]; break;
			case 2: nn = 3; nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[4]; break;
			case 3: nn = 3; nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; break;
			case 4: nn = 4; nf[0] = en[3]; nf[1] = en[2]; nf[2] = en[1]; nf[3] = en[0]; break;
		}
		break;
    case FE_TET4G4:
	case FE_TET4G1:
		nn = 3;
		switch (n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = nf[3] = en[3]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = nf[3] = en[3]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = nf[3] = en[3]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = nf[3] = en[0]; break;
		}
		break;
	case FE_TET10G4:
	case FE_TET10G8:
	case FE_TET10GL11:
		nn = 6;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; break;
		}
		break;
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
		nn = 7;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[8]; nf[5] = en[7]; nf[6] = en[11]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[5]; nf[4] = en[9]; nf[5] = en[8]; nf[6] = en[12]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[9]; nf[6] = en[13]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[5]; nf[4] = en[4]; nf[5] = en[6]; nf[6] = en[10]; break;
		}
		break;
	case FE_TET20G15:
		nn = 10;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[3]; nf[3] = en[4]; nf[4] = en[5]; nf[5] = en[12]; nf[6] = en[13]; nf[7] = en[10]; nf[8] = en[11]; nf[9] = en[16]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[3]; nf[3] = en[6]; nf[4] = en[7]; nf[5] = en[14]; nf[6] = en[15]; nf[7] = en[13]; nf[8] = en[14]; nf[9] = en[17]; break;
		case 2: nf[0] = en[2]; nf[1] = en[0]; nf[2] = en[3]; nf[3] = en[9]; nf[4] = en[8]; nf[5] = en[10]; nf[6] = en[11]; nf[7] = en[14]; nf[8] = en[15]; nf[9] = en[18]; break;
		case 3: nf[0] = en[2]; nf[1] = en[1]; nf[2] = en[0]; nf[3] = en[7]; nf[4] = en[6]; nf[5] = en[ 5]; nf[6] = en[ 4]; nf[7] = en[10]; nf[8] = en[ 8]; nf[9] = en[19]; break;
		}
		break;
    case FE_HEX20G8:
	case FE_HEX20G27:
		nn = 8;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[ 8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[ 9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[ 9]; nf[7] = en[ 8]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; break;
		}
		break;
	case FE_HEX27G27:
		nn = 9;
		switch(n)
		{
		case 0: nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[5]; nf[3] = en[4]; nf[4] = en[ 8]; nf[5] = en[17]; nf[6] = en[12]; nf[7] = en[16]; nf[8] = en[20]; break;
		case 1: nf[0] = en[1]; nf[1] = en[2]; nf[2] = en[6]; nf[3] = en[5]; nf[4] = en[ 9]; nf[5] = en[18]; nf[6] = en[13]; nf[7] = en[17]; nf[8] = en[21]; break;
		case 2: nf[0] = en[2]; nf[1] = en[3]; nf[2] = en[7]; nf[3] = en[6]; nf[4] = en[10]; nf[5] = en[19]; nf[6] = en[14]; nf[7] = en[18]; nf[8] = en[22]; break;
		case 3: nf[0] = en[3]; nf[1] = en[0]; nf[2] = en[4]; nf[3] = en[7]; nf[4] = en[11]; nf[5] = en[16]; nf[6] = en[15]; nf[7] = en[19]; nf[8] = en[23]; break;
		case 4: nf[0] = en[0]; nf[1] = en[3]; nf[2] = en[2]; nf[3] = en[1]; nf[4] = en[11]; nf[5] = en[10]; nf[6] = en[ 9]; nf[7] = en[ 8]; nf[8] = en[24]; break;
		case 5: nf[0] = en[4]; nf[1] = en[5]; nf[2] = en[6]; nf[3] = en[7]; nf[4] = en[12]; nf[5] = en[13]; nf[6] = en[14]; nf[7] = en[15]; nf[8] = en[25]; break;
		}
		break;
    case FE_SHELL_QUAD4G8:
    case FE_SHELL_QUAD4G12:
        nn = 4;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3];
		break;
    case FE_SHELL_QUAD8G18:
    case FE_SHELL_QUAD8G27:
        nn = 8;
        nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5]; nf[6] = en[6]; nf[7] = en[7];
        break;
    case FE_SHELL_TRI3G6:
    case FE_SHELL_TRI3G9:
        nn = 3;
		nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2];
		break;
    case FE_SHELL_TRI6G14:
    case FE_SHELL_TRI6G21:
        nn = 6;
        nf[0] = en[0]; nf[1] = en[1]; nf[2] = en[2]; nf[3] = en[3]; nf[4] = en[4]; nf[5] = en[5];
        break;
	}

	return nn;
}

//-----------------------------------------------------------------------------
//! return an element
FEElement* FEMesh::Element(int n)
{
	if (n < 0) return nullptr;
	for (int i = 0; i < Domains(); ++i)
	{
		FEDomain& dom = Domain(i);
		int NEL = dom.Elements();
		if (n < NEL) return &dom.ElementRef(n); 
		else n -= NEL;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Find a node from a given ID. return 0 if the node cannot be found.

FENode* FEMesh::FindNodeFromID(int nid)
{
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.GetID() == nid) return &node;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Find an element from a given ID. return 0 if the element cannot be found.

FEElement* FEMesh::FindElementFromID(int nid)
{
	if (m_LUT == 0) m_LUT = new FEElementLUT(*this);
	return m_LUT->Find(nid);
}

/*
FEElement* FEMesh::FindElementFromID(int nid)
{
	FEElement* pe = 0;

	for (int i=0; i<Domains(); ++i)
	{
		FEDomain& d = Domain(i);
		pe = d.FindElementFromID(nid);
		if (pe) return pe;
	}

	return pe;
}
*/


//-----------------------------------------------------------------------------
//! See if all elements are of a particular shape
bool FEMesh::IsType(FE_Element_Shape eshape)
{
	FEElementList elemList(*this);
	for (FEElementList::iterator it = elemList.begin(); it != elemList.end(); ++it)
	{
		FEElement& el = *it;
		if (el.Shape() != eshape) return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// Find the element in which point y lies
FESolidElement* FEMesh::FindSolidElement(vec3d y, double r[3])
{
	int ND = (int) m_Domain.size();
	for (int i=0; i<ND; ++i)
	{
		if (m_Domain[i]->Class() == FE_DOMAIN_SOLID)
		{
			FESolidDomain& bd = static_cast<FESolidDomain&>(*m_Domain[i]);
			FESolidElement* pe = bd.FindElement(y, r);
			if (pe) return pe;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FEMesh::ClearDomains()
{
	int N = Domains();
	for (int i = 0; i < N; ++i) delete m_Domain[i];
	m_Domain.clear();
	if (m_LUT) delete m_LUT; m_LUT = 0;
}

//-----------------------------------------------------------------------------
//! Rebuild the LUT
void FEMesh::RebuildLUT()
{
	if (m_LUT) delete m_LUT;
	m_LUT = new FEElementLUT(*this);
}

//-----------------------------------------------------------------------------
//! Calculate the surface representing the element boundaries
//! boutside : include all exterior facets
//! binside  : include all interior facets
FESurface* FEMesh::ElementBoundarySurface(bool boutside, bool binside)
{
	if ((boutside == false) && (binside == false)) return 0;

	// create the element neighbor list
	FEElemElemList EEL;
	EEL.Create(this);

	// get the number of elements in this mesh
	int NE = Elements();

	// count the number of facets we have to create
	int NF = 0;
	FEElementList EL(*this);
	FEElementList::iterator it = EL.begin();
	for (int i=0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = Faces(el);
		for (int j=0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if ((pen == 0) && boutside) ++NF;
			if ((pen != 0) && (el.GetID() < pen->GetID()) && binside ) ++NF;
		}
	}
	// create the surface
	FESurface* ps = new FESurface(GetFEModel());
	if (NF == 0) return ps;
	ps->Create(NF);

	// build the surface elements
	int face[FEElement::MAX_NODES];
	NF = 0;
	it = EL.begin();
	for (int i=0; i<NE; ++i, ++it)
	{
		FEElement& el = *it;
		int nf = Faces(el);
		for (int j=0; j<nf; ++j)
		{
			FEElement* pen = EEL.Neighbor(i, j);
			if (((pen == 0) && boutside)||
				((pen != 0) && (el.GetID() < pen->GetID()) && binside ))
			{
				FESurfaceElement& se = ps->Element(NF++);
				GetFace(el, j, face);

				switch (el.Shape())
				{
				case ET_HEX8:
					se.SetType(FE_QUAD4G4); 
					break;
				case ET_HEX20:
                    se.SetType(FE_QUAD8G9);
                    break;
				case ET_HEX27:
					se.SetType(FE_QUAD9G9);
					break;
				case ET_TET4:
					se.SetType(FE_TRI3G1); 
					break;
				case ET_TET10:
				case ET_TET15:
                    se.SetType(FE_TRI6G7);
                    break;
				default:
					assert(false);
				}
				
				se.m_elem[0] = &el;
				if (pen) se.m_elem[1] = pen;
				
				int nn = se.Nodes();
				for (int k=0; k<nn; ++k)
				{
					se.m_node[k] = face[k];
				}
			}
		}
	}

	// initialize the surface. 
	// This will set the local surface element ID's and also set the m_nelem IDs.
	ps->Init();

	// all done
	return ps;
}
FESurface* FEMesh::ElementBoundarySurface(std::vector<FEDomain*> domains, bool boutside, bool binside)
{
	if ((boutside == false) && (binside == false)) return nullptr;

	// create the element neighbor list
	FEElemElemList EEL;
	EEL.Create(this);

	// get the number of elements in this mesh
	int NE = Elements();

	// count the number of facets we have to create
	int NF = 0;

	for (int i = 0; i < domains.size(); i++)
	{
		for (int j = 0; j < domains[i]->Elements(); j++)
		{
			FEElement& el = domains[i]->ElementRef(j);
			int nf = Faces(el);
			for (int k = 0; k<nf; ++k)
			{
				FEElement* pen = EEL.Neighbor(el.GetID()-1, k);
				if ((pen == nullptr) && boutside) ++NF;
				else if (pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) && boutside) ++NF;
				if ((pen != nullptr) && (el.GetID() < pen->GetID()) && binside && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) != domains.end())) ++NF;
			}
		}
	}

	// create the surface
	FESurface* ps = new FESurface(GetFEModel());
	if (NF == 0) return ps;
	ps->Create(NF);

	// build the surface elements
	int face[FEElement::MAX_NODES];
	NF = 0;
	for (int i = 0; i < domains.size(); i++)
	{
		for (int j = 0; j < domains[i]->Elements(); j++)
		{
			FEElement& el = domains[i]->ElementRef(j);
			int nf = Faces(el);
			for (int k = 0; k < nf; ++k)
			{
				FEElement* pen = EEL.Neighbor(el.GetID()-1, k);
				if (((pen == nullptr) && boutside) ||
					(pen && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) == domains.end()) && boutside) ||
					((pen != nullptr) && (el.GetID() < pen->GetID()) && binside && (std::find(domains.begin(), domains.end(), pen->GetMeshPartition()) != domains.end())))
				{
					FESurfaceElement& se = ps->Element(NF++);
					GetFace(el, k, face);

					switch (el.Shape())
					{
					case ET_HEX8:
						se.SetType(FE_QUAD4G4);
						break;
					case ET_HEX20:
						se.SetType(FE_QUAD8G9);
						break;
					case ET_HEX27:
						se.SetType(FE_QUAD9G9);
						break;
					case ET_TET4:
						se.SetType(FE_TRI3G1);
						break;
					case ET_TET10:
					case ET_TET15:
						se.SetType(FE_TRI6G7);
						break;
					default:
						assert(false);
					}

					se.m_elem[0] = &el;
					if (pen) se.m_elem[1] = pen;

					int nn = se.Nodes();
					for (int p = 0; p < nn; ++p)
					{
						se.m_node[p] = face[p];
					}
				}
			}
		}
	}

	// initialize the surface. 
	// This will set the local surface element ID's and also set the m_nelem IDs.
	ps->Init();

	// all done
	return ps;
}

//-----------------------------------------------------------------------------
//! Retrieve the nodal coordinates of an element in the reference configuration.
void FEMesh::GetInitialNodalCoordinates(const FEElement& el, vec3d* node)
{
	const int neln = el.Nodes();
	for (int i=0; i<neln; ++i) node[i] = Node(el.m_node[i]).m_r0;
}

//-----------------------------------------------------------------------------
//! Retrieve the nodal coordinates of an element in the current configuration.
void FEMesh::GetNodalCoordinates(const FEElement& el, vec3d* node)
{
	const int neln = el.Nodes();
	for (int i=0; i<neln; ++i) node[i] = Node(el.m_node[i]).m_rt;
}

//=============================================================================
FEElementLUT::FEElementLUT(FEMesh& mesh)
{
	// get the ID ranges
	m_minID = -1;
	m_maxID = -1;
	int NDOM = mesh.Domains();
	for (int i=0; i<NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j=0; j<NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int eid = el.GetID();
			if ((eid < m_minID) || (m_minID == -1)) m_minID = eid;
			if ((eid > m_maxID) || (m_maxID == -1)) m_maxID = eid;
		}
	}

	// allocate size
	int nsize = m_maxID - m_minID + 1;
	m_elem.resize(nsize, (FEElement*) 0);

	// fill the table
	for (int i = 0; i<NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j<NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int eid = el.GetID();
			m_elem[eid - m_minID] = &el;
		}
	}
}

// Find an element from its ID
FEElement* FEElementLUT::Find(int nid)
{
	if ((nid < m_minID) || (nid > m_maxID)) return 0;
	return m_elem[nid - m_minID];
}

// update the domains of the mesh
void FEMesh::Update(const FETimeInfo& tp)
{
	for (int i = 0; i<Domains(); ++i)
	{
		FEDomain& dom = Domain(i);
		if (dom.IsActive()) dom.Update(tp);
	}
}


//-----------------------------------------------------------------------------
void FEMesh::ClearDataArrays()
{
	// clear the surface maps
	for (int i = 0; i<(int)m_DataArray.size(); ++i) delete m_DataArray[i].second;
	m_DataArray.clear();
}

//-----------------------------------------------------------------------------
void FEMesh::AddDataArray(const std::string& name, FEDataArray* map)
{
	m_DataArray.push_back(pair<string, FEDataArray*>(name, map));
}

//-----------------------------------------------------------------------------
FEDataArray* FEMesh::FindDataArray(const std::string& map)
{
	for (int i = 0; i<(int)m_DataArray.size(); ++i)
	{
		if (m_DataArray[i].first == map) return m_DataArray[i].second;
	}
	return 0;
}

//-----------------------------------------------------------------------------
int FEMesh::DataArrays() const
{
	return (int)m_DataArray.size();
}

//-----------------------------------------------------------------------------
FEDataArray* FEMesh::GetDataArray(int i)
{
	return m_DataArray[i].second;
}

//-----------------------------------------------------------------------------
std::string FEMesh::DataArrayName(int i)
{
	return m_DataArray[i].first;
}
