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
#include "stdafx.h"
#include "FERefineMesh.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEEdgeList.h>
#include <FECore/FEElementList.h>
#include <FECore/FEFaceList.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FESurfacePairConstraint.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEDomainMap.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/DumpMemStream.h>
#include <FECore/log.h>
#include "FELeastSquaresInterpolator.h"
#include "FEMeshShapeInterpolator.h"
#include "FEDomainShapeInterpolator.h"

BEGIN_FECORE_CLASS(FERefineMesh, FEMeshAdaptor)
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elements");
	ADD_PARAMETER(m_bmap_data, "map_data");
	ADD_PARAMETER(m_nnc      , "nnc");
	ADD_PARAMETER(m_nsdim  , "nsdim");
	ADD_PARAMETER(m_transferMethod, "transfer_method");
END_FECORE_CLASS();

FERefineMesh::FERefineMesh(FEModel* fem) : FEMeshAdaptor(fem), m_topo(nullptr)
{
	m_meshCopy = nullptr;
	m_bmap_data = false;
	m_transferMethod = TRANSFER_SHAPE;
	m_nnc = 8;
	m_nsdim = 3;

	m_maxiter = -1;
	m_maxelem = -1;
}

FERefineMesh::~FERefineMesh()
{
	if (m_meshCopy) delete m_meshCopy;
	m_meshCopy = nullptr;

	ClearMapData();
}

// Apply mesh refinement
bool FERefineMesh::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// check for max iterations
	if ((m_maxiter > 0) && (iteration >= m_maxiter))
	{
		feLog("Skipping refinement: Max iterations reached.");
		return false;
	}

	// see if we reached max elements
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("Skipping refinement: Element limit reached.\n");
		return false;
	}

	// build the mesh-topo
	if (BuildMeshTopo() == false)
	{
		throw std::runtime_error("Error building topo structure.");
	}

	// build data maps
	feLog("-- Building map data:\n");
	if (BuildMapData() == false)
	{
		throw std::runtime_error("Failed mapping data.");
	}

	// refine the mesh (This is done by sub-classes)
	feLog("-- Starting Mesh refinement.\n");
	if (RefineMesh() == false)
	{
		feLog("Nothing to refine.");
		return false;
	}
	feLog("-- Mesh refinement completed.\n");

	// map data to new mesh
	feLog("-- Transferring map data to new mesh:\n");
	TransferMapData();

	// update the model
	UpdateModel();

	// print some mesh statistics
	int NN = mesh.Nodes();
	int NE = mesh.Elements();
	feLog("\n Mesh Statistics:\n");
	feLog(" \tNumber of nodes    : %d\n", NN);
	feLog(" \tNumber of elements : %d\n", NE);
	feLog("\n");

	// all done!
	return true;
}

void FERefineMesh::ClearMapData()
{
	// clear domain maps
	for (size_t i = 0; i < m_domainMapList.size(); ++i)
	{
		std::vector<FEDomainMap*>& map_i = m_domainMapList[i];
		for (size_t j = 0; j < map_i.size(); ++j) delete map_i[j];
	}
	m_domainMapList.clear();

	// clear user maps
	for (int i = 0; i < m_userDataList.size(); ++i) delete m_userDataList[i];
	m_userDataList.clear();
}

bool FERefineMesh::BuildMeshTopo()
{
	FEModel& fem = *GetFEModel();
	if (m_topo) { delete m_topo; m_topo = nullptr; }
	m_topo = new FEMeshTopo;
	return m_topo->Create(&fem.GetMesh());
}

void FERefineMesh::CopyMesh()
{
	if (m_meshCopy) delete m_meshCopy;

	m_meshCopy = new FEMesh(nullptr);
	FEMesh& mesh = GetFEModel()->GetMesh();
	m_meshCopy->CopyFrom(mesh);
}

FEDomainMap* createElemDataMap(FEModel& fem, FEDomain& dom, vector<vec3d>& nodePos, FEDomainMap* map, FEMeshDataInterpolator* dataMapper)
{
	assert(map->StorageFormat() == Storage_Fmt::FMT_NODE);

	FEDataType dataType = map->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		return nullptr;
	}

	// count nr of integration points
	int NMP = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		NMP += el.GaussPoints();
	}

	int N0 = nodePos.size();

	// create new domain map
	FEDomainMap* elemData = new FEDomainMap(map->DataType(), Storage_Fmt::FMT_MATPOINTS);
	FEElementSet* eset = new FEElementSet(&fem);
	eset->Create(&dom);
	elemData->Create(eset);

	vector<double> srcData(N0);
	vector<double> trgData(NMP);

	vector< vector<double> > mappedData(NMP, vector<double>(9, 0.0));

	// loop over all the new nodes
	for (int l = 0; l < dataSize; ++l)
	{
		for (int i = 0; i < N0; ++i)
		{
			double vm = 0.0;
			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: vm = map->value<double>(0, i); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v = map->value<vec3d>(0, i);
				if (l == 0) vm = v.x;
				if (l == 1) vm = v.y;
				if (l == 2) vm = v.z;
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
				mat3d v = map->value<mat3d>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
				mat3ds v = map->value<mat3ds>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			default:
				assert(false);
			}
			srcData[i] = vm;
		}

		dataMapper->Map(trgData, [&srcData](int sourcePoint) {
			return srcData[sourcePoint];
		});

		for (int i = 0; i < NMP; ++i)
		{
			mappedData[i][l] = trgData[i];
		}
	}

	// write mapped data to domain map
	int n = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		int nint = el.GaussPoints();
		for (int j = 0; j < nint; ++j)
		{
			vector<double>& vj = mappedData[n++];

			for (int l = 0; l < dataSize; ++l)
			{
				switch (dataType)
				{
				case FEDataType::FE_DOUBLE: elemData->setValue(i, j, vj[0]); break;
				case FEDataType::FE_VEC3D:
				{
					vec3d v;
					v.x = vj[0];
					v.y = vj[1];
					v.z = vj[2];
					elemData->setValue(i, j, v);
				}
				break;
				case FEDataType::FE_MAT3D:
				{
					mat3d v;
					v(0, 0) = vj[0]; v(0, 1) = vj[1]; v(0, 2) = vj[2];
					v(1, 0) = vj[3]; v(1, 1) = vj[4]; v(1, 2) = vj[5];
					v(2, 0) = vj[6]; v(2, 1) = vj[7]; v(2, 2) = vj[8];
					elemData->setValue(i, j, v);
				}
				break;
				case FEDataType::FE_MAT3DS:
				{
					mat3ds v;
					v(0, 0) = vj[0];
					v(0, 1) = vj[1];
					v(0, 2) = vj[2];
					v(1, 1) = vj[3];
					v(1, 2) = vj[4];
					v(2, 2) = vj[5];
					elemData->setValue(i, j, v);
				}
				break;
				default:
					assert(false);
				}
			}
		}
	}

	return elemData;
}

bool createNodeDataMap(FEDomain& dom, FEDomainMap* map, FEDomainMap* nodeMap)
{
	FEDataType dataType = map->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		return false;
	}

	// temp storage 
	double si[FEElement::MAX_INTPOINTS * 9];
	double sn[FEElement::MAX_NODES * 9];

	// allocate node data
	int NN = dom.Nodes();
	vector<double> nodeData(NN*dataSize);

	// build tag list
	vector<int> tag(NN, 0);
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		for (int k = 0; k < ne; ++k)
		{
			tag[e.m_lnode[k]]++;
		}
	}

	// get the data format
	int dataFormat = map->StorageFormat();
	if ((dataFormat != FMT_MATPOINTS) && (dataFormat != FMT_MULT)) return false;

	// loop over all elements
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();

		int ni = (dataFormat == FMT_MATPOINTS ? e.GaussPoints() : ne);

		for (int j = 0; j < dataSize; ++j)
		{
			// get the integration point values
			for (int k = 0; k < ni; ++k)
			{
				switch (dataType)
				{
				case FEDataType::FE_DOUBLE:
					si[k] = map->value<double>(i, k);
					break;
				case FEDataType::FE_VEC3D:
				{
					vec3d v = map->value<vec3d>(i, k);
					if (j == 0) si[k] = v.x;
					if (j == 1) si[k] = v.y;
					if (j == 2) si[k] = v.z;
				}
				break;
				case FEDataType::FE_MAT3D:
				{
					mat3d v = map->value<mat3d>(i, k);
					int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
					si[k] = v(LUT[j][0], LUT[j][1]);
				}
				break;
				case FEDataType::FE_MAT3DS:
				{
					mat3ds v = map->value<mat3ds>(i, k);
					int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
					si[k] = v(LUT[j][0], LUT[j][1]);
				}
				break;
				}
			}

			// project to nodes
			if (dataFormat == FMT_MATPOINTS)
			{
				e.project_to_nodes(si, sn);
			}
			else
			{
				for (int k = 0; k < ne; ++k) sn[k] = si[k];
			}

			for (int k = 0; k < ne; ++k)
			{
				nodeData[e.m_lnode[k] * dataSize + j] += sn[k];
			}
		}
	}

	// normalize data
	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0)
		{
			for (int j = 0; j < dataSize; ++j)
				nodeData[i*dataSize + j] /= (double)tag[i];
		}
	}

	// check node data map
	if (nodeMap->StorageFormat() != Storage_Fmt::FMT_NODE) return false;
	if (nodeMap->DataType() != dataType) return false;
	if (nodeMap->DataCount() != NN) return false;

	// assign data
	for (int i = 0; i < NN; ++i)
	{
		switch (dataType)
		{
		case FEDataType::FE_DOUBLE: nodeMap->setValue(i, nodeData[i]); break;
		case FEDataType::FE_VEC3D:
		{
			vec3d v;
			v.x = nodeData[i*dataSize];
			v.y = nodeData[i*dataSize + 1];
			v.z = nodeData[i*dataSize + 2];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3D:
		{
			mat3d v;
			v(0, 0) = nodeData[i*dataSize]; v(0, 1) = nodeData[i*dataSize + 1]; v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 0) = nodeData[i*dataSize + 3]; v(1, 1) = nodeData[i*dataSize + 4]; v(1, 2) = nodeData[i*dataSize + 5];
			v(2, 0) = nodeData[i*dataSize + 6]; v(2, 1) = nodeData[i*dataSize + 7]; v(2, 2) = nodeData[i*dataSize + 8];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3DS:
		{
			mat3ds v;
			v(0, 0) = nodeData[i*dataSize];
			v(0, 1) = nodeData[i*dataSize + 1];
			v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 1) = nodeData[i*dataSize + 3];
			v(1, 2) = nodeData[i*dataSize + 4];
			v(2, 2) = nodeData[i*dataSize + 5];
			nodeMap->setValue(i, v);
		}
		break;
		}
	}

	return true;
}

void NodeToElemData(FEModel& fem, FEDomain& dom, FEDomainMap* nodeMap, FEDomainMap* elemMap, FEMeshDataInterpolator* dataMapper)
{
	assert(nodeMap->StorageFormat() == Storage_Fmt::FMT_NODE);

	FEDataType dataType = nodeMap->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		throw std::runtime_error("Error in FEMMGRemesh::MMG::NodeToElemData");
		return;
	}

	// count nr of points
	int NN = nodeMap->DataCount();

	// count nr of target points
	int NE = dom.Elements();
	int NP = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		NP += el.Nodes();
	}

	vector<double> srcData(NN);
	vector<double> trgData(NP);

	vector< vector<double> > mappedData(NP, vector<double>(9, 0.0));

	// loop over all the new nodes
	for (int l = 0; l < dataSize; ++l)
	{
		for (int i = 0; i < NN; ++i)
		{
			double vm = 0.0;
			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: vm = nodeMap->value<double>(0, i); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v = nodeMap->value<vec3d>(0, i);
				if (l == 0) vm = v.x;
				if (l == 1) vm = v.y;
				if (l == 2) vm = v.z;
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
				mat3d v = nodeMap->value<mat3d>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
				mat3ds v = nodeMap->value<mat3ds>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			default:
				assert(false);
			}
			srcData[i] = vm;
		}

		dataMapper->Map(trgData, [&srcData](int sourcePoint) {
			return srcData[sourcePoint];
		});

		for (int i = 0; i < NP; ++i)
		{
			mappedData[i][l] = trgData[i];
		}
	}

	// create new element set
	FEElementSet* eset = new FEElementSet(&fem);
	eset->Create(&dom);
	elemMap->Create(eset);

	// write mapped data to domain map
	int n = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k = 0; k < el.Nodes(); ++k)
		{
			vector<double>& vj = mappedData[n++];

			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: elemMap->setValue(i, k, vj[0]); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v;
				v.x = vj[0];
				v.y = vj[1];
				v.z = vj[2];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				mat3d v;
				v(0, 0) = vj[0]; v(0, 1) = vj[1]; v(0, 2) = vj[2];
				v(1, 0) = vj[3]; v(1, 1) = vj[4]; v(1, 2) = vj[5];
				v(2, 0) = vj[6]; v(2, 1) = vj[7]; v(2, 2) = vj[8];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				mat3ds v;
				v(0, 0) = vj[0];
				v(0, 1) = vj[1];
				v(0, 2) = vj[2];
				v(1, 1) = vj[3];
				v(1, 2) = vj[4];
				v(2, 2) = vj[5];
				elemMap->setValue(i, k, v);
			}
			break;
			default:
				assert(false);
			}
		}
	}
}

bool FERefineMesh::BuildMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// make a copy of the old mesh
	// we need it for mapping data
	CopyMesh();

	// clear all map data
	ClearMapData();
	m_domainMapList.clear();
	m_domainMapList.resize(mesh.Domains());

	// only map domain data if requested
	if (m_bmap_data)
	{
		if (BuildDomainMapData() == false)
		{
			return false;
		}
	}

	// do the same thing for the user-defined mesh data
	return BuildUserMapData();
}

bool FERefineMesh::BuildDomainMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		feLog(" Processing domain: %s\n", dom.GetName().c_str());

		if (BuildDomainMapData(dom, i) == false)
		{
			return false;
		}
	}

	return true;
}

bool FERefineMesh::BuildDomainMapData(FEDomain& dom, int domIndex)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// write all material point data to a data stream
	DumpMemStream ar(fem);
	ar.Open(true, true);
	ar.WriteTypeInfo(true);

	// loop over all integration points
	int totalPoints = 0;
	for (int j = 0; j < dom.Elements(); ++j)
	{
		FEElement& el = dom.ElementRef(j);
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)
		{
			FEMaterialPoint* mp = el.GetMaterialPoint(k);
			mp->Serialize(ar);
		}
		totalPoints += nint;
	}

	// figure out how much data was written for each material point
	size_t bytes = ar.bytesSerialized();

	size_t bytesPerPoint = bytes / totalPoints;
	assert((bytes%totalPoints) == 0);

	// re-open for reading
	ar.Open(false, true);

	// create an element set (we need this for the domain map below)
	FEDomain& oldDomain = m_meshCopy->Domain(domIndex);
	FEElementSet* elemSet = new FEElementSet(&fem);
	elemSet->Create(&oldDomain);

	// next, we need to figure out the datamaps for each data item
	vector<FEDomainMap*> mapList;
	DumpStream::DataBlock d;
	while (ar.bytesSerialized() < bytesPerPoint)
	{
		ar.readBlock(d);

		const char* typeStr = nullptr;
		FEDomainMap* map = nullptr;
		switch (d.dataType())
		{
		case TypeID::TYPE_DOUBLE: map = new FEDomainMap(FEDataType::FE_DOUBLE, Storage_Fmt::FMT_MATPOINTS); typeStr = "double"; break;
		case TypeID::TYPE_VEC3D: map = new FEDomainMap(FEDataType::FE_VEC3D, Storage_Fmt::FMT_MATPOINTS); typeStr = "vec3d"; break;
		case TypeID::TYPE_MAT3D: map = new FEDomainMap(FEDataType::FE_MAT3D, Storage_Fmt::FMT_MATPOINTS); typeStr = "mat3d"; break;
		case TypeID::TYPE_MAT3DS: map = new FEDomainMap(FEDataType::FE_MAT3DS, Storage_Fmt::FMT_MATPOINTS); typeStr = "mat3ds"; break;
		default:
			assert(false);
			throw std::runtime_error("Error in mapping data.");
		}

		map->Create(elemSet, 0.0);

		feLog("\tData map %d: %s\n", mapList.size(), typeStr);

		mapList.push_back(map);
	}

	feLog(" %d data maps identified.\n", mapList.size());

	// rewind for processing
	ar.Open(false, true);

	for (int j = 0; j < dom.Elements(); ++j)
	{
		FEElement& el = dom.ElementRef(j);
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)
		{
			int m = 0;
			size_t bytesRead = 0;
			while (bytesRead < bytesPerPoint)
			{
				size_t size0 = ar.bytesSerialized();
				bool b = ar.readBlock(d); assert(b);
				size_t size1 = ar.bytesSerialized();
				bytesRead += size1 - size0;

				FEDomainMap* map = mapList[m];

				switch (d.dataType())
				{
				case TypeID::TYPE_DOUBLE: { double v = d.value<double>(); map->setValue(j, k, v); } break;
				case TypeID::TYPE_MAT3D: { mat3d  v = d.value<mat3d >(); map->setValue(j, k, v); } break;
				case TypeID::TYPE_MAT3DS: { mat3ds v = d.value<mat3ds>(); map->setValue(j, k, v); } break;
				}

				m++;
				assert(m <= mapList.size());
			}
		}
	}

	// Now, we need to project all the data onto the nodes
	for (int j = 0; j < mapList.size(); ++j)
	{
		feLog("\tProcessing data map %d ...", j);
		FEDomainMap* elemMap = mapList[j];

		FEDataType dataType = elemMap->DataType();
		FEDomainMap* nodeMap = new FEDomainMap(dataType, Storage_Fmt::FMT_NODE);

		FEElementSet* elset = const_cast<FEElementSet*>(elemMap->GetElementSet());
		nodeMap->Create(elset);

		bool bret = createNodeDataMap(dom, elemMap, nodeMap); assert(bret);
		m_domainMapList[domIndex].push_back(nodeMap);
		feLog("done.\n");
	}
    
    return true;
}

bool FERefineMesh::BuildUserMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	m_userDataList.clear();
	int dataMaps = mesh.DataMaps();
	if (dataMaps > 0) feLog(" Processing user data maps:\n");

	for (int i = 0; i < dataMaps; ++i)
	{
		FEDataMap* dataMap = mesh.GetDataMap(i);

		// process domain map
		FEDomainMap* dmap = dynamic_cast<FEDomainMap*>(dataMap);
		if (dmap)
		{
			FEDataType dataType = dmap->DataType();
			if ((dataType != FEDataType::FE_DOUBLE) && 
				(dataType != FEDataType::FE_VEC3D)) return false;

			feLog("\tProcessing user data map \"%s\" ...", dmap->GetName().c_str());

			const FEElementSet* elset = dmap->GetElementSet();
			const FEDomainList& domainList = elset->GetDomainList();
			if (domainList.Domains() != 1) return false;
			FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));

			FEElementSet* oldElemSet = m_meshCopy->FindElementSet(dom.GetName()); assert(oldElemSet);

			FEDomainMap* nodeMap = new FEDomainMap(dataType, Storage_Fmt::FMT_NODE);
			nodeMap->Create(oldElemSet);

			// create a node data map of this domain map
			createNodeDataMap(dom, dmap, nodeMap);
			m_userDataList.push_back(nodeMap);

			feLog("done.\n");
		}

		FESurfaceMap* smap = dynamic_cast<FESurfaceMap*>(dataMap);
		if (smap)
		{
			assert(false);
		}
	}

	return true;
}

// Transfer data to new mesh
void FERefineMesh::TransferMapData()
{
	// transfer domain data
	if (m_bmap_data) TransferDomainMapData();

	// transfer user-defined maps
	TransferUserMapData();
}

// Transfer domain data back to the new mesh
void FERefineMesh::TransferDomainMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	if (m_bmap_data && m_domainMapList.size())
	{
		for (int i = 0; i < mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);

			std::vector<FEDomainMap*>& nodeMap_i = m_domainMapList[i];
			int mapCount = nodeMap_i.size();

			if (mapCount > 0)
			{
				feLog(" Mapping data for domain \"%s\":\n", dom.GetName().c_str());

				// we need an element set for the domain maps below
				FEElementSet* elemSet = new FEElementSet(&fem);
				elemSet->Create(&dom);

				// build source point list
				FEDomain& oldDomain = m_meshCopy->Domain(i);
				vector<vec3d> srcPoints; srcPoints.reserve(oldDomain.Nodes());
				for (int i = 0; i < oldDomain.Nodes(); ++i)
				{
					vec3d r = oldDomain.Node(i).m_r0;
					srcPoints.push_back(r);
				}

				// build target node list
				vector<vec3d> trgPoints; trgPoints.reserve(dom.Elements());
				for (int i = 0; i < dom.Elements(); ++i)
				{
					FEElement& el = dom.ElementRef(i);
					int nint = el.GaussPoints();
					for (int j = 0; j < nint; ++j)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(j);
						vec3d r = mp.m_r0;
						trgPoints.push_back(r);
					}
				}

				// set up mapper
				FEMeshDataInterpolator* mapper = nullptr;
				switch (m_transferMethod)
				{
				case TRANSFER_SHAPE:
				{
					FEDomain* oldDomain = &m_meshCopy->Domain(i);
					FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain);
					dsm->SetTargetPoints(trgPoints);
					mapper = dsm;
				}
				break;
				case TRANSFER_MLQ:
				{
					FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
					MLQ->SetNearestNeighborCount(m_nnc);
					MLQ->SetDimension(m_nsdim);
					MLQ->SetSourcePoints(srcPoints);
					MLQ->SetTargetPoints(trgPoints);
					mapper = MLQ;
				}
				break;
				default:
					assert(false);
					return;
				}
				if (mapper->Init() == false)
				{
					assert(false);
					throw std::runtime_error("Failed to initialize LLQ");
				}

				// loop over all the domain maps
				vector<FEDomainMap*> elemMapList(mapCount);
				for (int j = 0; j < mapCount; ++j)
				{
					feLog("\tMapping map %d ...", j);
					FEDomainMap* nodeMap = nodeMap_i[j];

					// map node data to integration points
					FEDomainMap* elemMap = createElemDataMap(fem, dom, srcPoints, nodeMap, mapper);

					elemMapList[j] = elemMap;
					feLog("done.\n");
				}

				// now we need to reconstruct the data stream
				DumpMemStream ar(fem);
				ar.Open(true, true);

				for (int j = 0; j < dom.Elements(); ++j)
				{
					FEElement& el = dom.ElementRef(j);
					int nint = el.GaussPoints();

					for (int k = 0; k < nint; ++k)
					{
						for (int l = 0; l < mapCount; ++l)
						{
							FEDomainMap* map = elemMapList[l];
							switch (map->DataType())
							{
							case FEDataType::FE_DOUBLE: { double v = map->value<double>(j, k); ar << v; } break;
							case FEDataType::FE_VEC3D: { vec3d  v = map->value<vec3d >(j, k); ar << v; } break;
							case FEDataType::FE_MAT3D: { mat3d  v = map->value<mat3d >(j, k); ar << v; } break;
							case FEDataType::FE_MAT3DS: { mat3ds v = map->value<mat3ds>(j, k); ar << v; } break;
							default:
								assert(false);
							}
						}
					}
				}

				// time to serialize everything back to the new integration points
				ar.Open(false, true);
				for (int j = 0; j < dom.Elements(); ++j)
				{
					FEElement& el = dom.ElementRef(j);
					int nint = el.GaussPoints();

					for (int k = 0; k < nint; ++k)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(k);
						mp.Serialize(ar);
					}
				}
			}
		}
	}
}

// Transfer user data maps to new mesh
void FERefineMesh::TransferUserMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// map mesh data
	for (int i = 0; i < m_userDataList.size(); ++i)
	{
		FEDomainMap* elemMap = dynamic_cast<FEDomainMap*>(mesh.GetDataMap(i));
		FEDomainMap* nodeMap = m_userDataList[i];

		feLog("\tMapping user map \"%s\" ...", elemMap->GetName().c_str());

		// build the source point list
		vector<vec3d> srcPoints;
		const FEElementSet* oldSet = nodeMap->GetElementSet();
		FEMesh* oldMesh = oldSet->GetMesh(); assert(oldMesh == m_meshCopy);
		FENodeList oldNodeList = oldSet->GetNodeList();
		srcPoints.resize(oldNodeList.Size());
		for (int i = 0; i < oldNodeList.Size(); ++i)
		{
			FENode& oldNode_i = oldMesh->Node(oldNodeList[i]);
			vec3d r = oldNode_i.m_r0;
			srcPoints[i] = r;
		}

		// build target points.
		const FEDomainList& domainList = elemMap->GetElementSet()->GetDomainList();
		FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));
		int NE = dom.Elements();
		vector<vec3d> trgPoints; trgPoints.reserve(NE);
		for (int n = 0; n < NE; ++n)
		{
			FEElement& el = dom.ElementRef(n);
			int ne = el.Nodes();
			for (int l = 0; l < ne; ++l)
			{
				vec3d r = mesh.Node(el.m_node[l]).m_r0;
				trgPoints.push_back(r);
			}
		}

		// set up mapper
		FEMeshDataInterpolator* mapper = nullptr;
		switch (m_transferMethod)
		{
		case TRANSFER_SHAPE:
		{
			FEDomain* oldDomain = m_meshCopy->FindDomain(dom.GetName());
			FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain);
			dsm->SetTargetPoints(trgPoints);
			mapper = dsm;
		}
		break;
		case TRANSFER_MLQ:
		{
			FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
			MLQ->SetNearestNeighborCount(m_nnc);
			MLQ->SetDimension(m_nsdim);
			MLQ->SetSourcePoints(srcPoints);
			MLQ->SetTargetPoints(trgPoints);
			mapper = MLQ;
		}
		break;
		default:
			assert(false);
			return;
		}
		if (mapper->Init() == false)
		{
			assert(false);
			throw std::runtime_error("Failed to initialize LLQ");
		}

		// map node data to integration points
		NodeToElemData(fem, dom, nodeMap, elemMap, mapper);

		feLog("done.\n");
	}
}
