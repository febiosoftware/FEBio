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
#include "VTKPlotFile.h"
#include <FECore/FEModel.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/FEDomain.h>
#include <sstream>

enum VTK_CELLTYPE {
	VTK_VERTEX = 1,
	VTK_POLY_VERTEX = 2,
	VTK_LINE = 3,
	VTK_POLY_LINE = 4,
	VTK_TRIANGLE = 5,
	VTK_TRIANGLE_STRIP = 6,
	VTK_POLYGON = 7,
	VTK_PIXEL = 8,
	VTK_QUAD = 9,
	VTK_TETRA = 10,
	VTK_VOXEL = 11,
	VTK_HEXAHEDRON = 12,
	VTK_WEDGE = 13,
	VTK_PYRAMID = 14,
	VTK_QUADRATIC_EDGE = 21,
	VTK_QUADRATIC_TRIANGLE = 22,
	VTK_QUADRATIC_QUAD = 23,
	VTK_QUADRATIC_TETRA = 24,
	VTK_QUADRATIC_HEXAHEDRON = 25,
	VTK_QUADRATIC_WEDGE = 26,
	VTK_QUADRATIC_PYRAMID = 27
};


VTKPlotFile::VTKPlotFile(FEModel* fem) : PlotFile(fem)
{
	m_fp = nullptr;
	m_count = 0;
	m_valid = false;
}

//! Open the plot database
bool VTKPlotFile::Open(const char* szfile)
{
	m_filename = szfile;
	size_t n = m_filename.rfind('.');
	if (n != std::string::npos) m_filename.erase(n, std::string::npos);

	BuildDictionary();
	m_valid = true;
	return true;
}

//! Open for appending
bool VTKPlotFile::Append(const char* szfile)
{
	m_filename = szfile;
	size_t n = m_filename.rfind('.');
	if (n != std::string::npos) m_filename.erase(n, std::string::npos);

	BuildDictionary();
	m_valid = true;
	return true;
}

//! see if the plot file is valid
bool VTKPlotFile::IsValid() const
{
	return m_valid;
}

void VTKPlotFile::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar& m_count;
}

//! Write current FE state to plot database
bool VTKPlotFile::Write(float ftime, int flag)
{
	FEModel& fem = *GetFEModel();

	std::stringstream ss;
	ss << m_filename << "." << m_count++ << ".vtk";
	string fileName = ss.str();
	
	m_fp = fopen(fileName.c_str(), "wt");
	if (m_fp == nullptr) return false;

	WriteHeader();
	WritePoints();
	WriteCells();
	WritePointData();
	WriteCellData();

	fclose(m_fp);
	m_fp = nullptr;

	return true;
}

//-----------------------------------------------------------------------------
void VTKPlotFile::WriteHeader()
{
	FEModel& fem = *GetFEModel();
	fprintf(m_fp, "%s\n", "# vtk DataFile Version 3.0");
	fprintf(m_fp, "%s %g\n", "time", fem.GetCurrentTime());
	fprintf(m_fp, "%s\n", "ASCII");
	fprintf(m_fp, "%s\n", "DATASET UNSTRUCTURED_GRID");
}

//-----------------------------------------------------------------------------
void VTKPlotFile::WritePoints()
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int nodes = m.Nodes();
	fprintf(m_fp, "POINTS %d float\n", nodes);
	for (int j = 0; j < nodes; j += 3)
	{
		for (int k = 0; k < 3 && j + k < nodes; k++)
		{
			FENode& nd = m.Node(j + k);
			vec3d& r = nd.m_r0;
			fprintf(m_fp, "%lg %lg %lg ", r.x, r.y, r.z);
		}
		fprintf(m_fp, "\n");
	}
	fprintf(m_fp, "%s\n", "");
}

//-----------------------------------------------------------------------------
void VTKPlotFile::WriteCells()
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NE = m.Elements();
    int nsize = 0;
	for (int j = 0; j<NE; ++j)
        nsize += m.Element(j)->Nodes() + 1;

	// Write CELLS
    fprintf(m_fp, "CELLS %d %d\n", NE, nsize);
    for (int j=0; j<NE; ++j)
    {
		FEElement& el = *m.Element(j);
        fprintf(m_fp, "%d ", el.Nodes());
        for (int k=0; k<el.Nodes(); ++k) fprintf(m_fp, "%d ", el.m_node[k]);
        fprintf(m_fp, "\n");
    }
        
	// Write CELL_TYPES
    fprintf(m_fp, "\nCELL_TYPES %d\n", NE);
	for (int j = 0; j<m.Elements(); ++j)
    {
		FEElement& el = *m.Element(j);
		int vtk_type;
        switch (el.Shape()) {
            case ET_HEX8   : vtk_type = VTK_HEXAHEDRON; break;
            case ET_TET4   : vtk_type = VTK_TETRA; break;
            case ET_PENTA6 : vtk_type = VTK_WEDGE; break;
            case ET_PYRA5  : vtk_type = VTK_PYRAMID; break;
            case ET_QUAD4  : vtk_type = VTK_QUAD; break;
            case ET_TRI3   : vtk_type = VTK_TRIANGLE; break;
            case ET_TRUSS2 : vtk_type = VTK_LINE; break;
            case ET_HEX20  : vtk_type = VTK_QUADRATIC_HEXAHEDRON; break;
            case ET_QUAD8  : vtk_type = VTK_QUADRATIC_QUAD; break;
//            case ET_BEAM3  : vtk_type = VTK_QUADRATIC_EDGE; break;
            case ET_TET10  : vtk_type = VTK_QUADRATIC_TETRA; break;
            case ET_TET15  : vtk_type = VTK_QUADRATIC_TETRA; break;
            case ET_PENTA15: vtk_type = VTK_QUADRATIC_WEDGE; break;
            case ET_HEX27  : vtk_type = VTK_QUADRATIC_HEXAHEDRON; break;
            case ET_PYRA13 : vtk_type = VTK_QUADRATIC_PYRAMID; break;
            case ET_TRI6   : vtk_type = VTK_QUADRATIC_TRIANGLE; break;
            case ET_QUAD9  : vtk_type = VTK_QUADRATIC_QUAD; break;
            default: vtk_type = -1; break;
        }
            
        fprintf(m_fp, "%d\n", vtk_type);
    }
}

void VTKPlotFile::WriteScalarData(std::vector<float>& val, const std::string& name)
{
	fprintf(m_fp, "%s %s %s\n", "SCALARS", name.c_str(), "float");
	fprintf(m_fp, "%s %s\n", "LOOKUP_TABLE", "default");
	for (int i = 0; i < val.size(); ++i) fprintf(m_fp, "%g\n", val[i]);
}

void VTKPlotFile::WriteVectorData(std::vector<float>& val, const std::string& name)
{
	fprintf(m_fp, "%s %s %s\n", "VECTORS", name.c_str(), "float");
	for (int i = 0; i < val.size(); i += 3) fprintf(m_fp, "%g %g %g\n", val[i], val[i + 1], val[i + 2]);
}

void VTKPlotFile::WriteMat3FSData(std::vector<float>& val, const std::string& name)
{
	fprintf(m_fp, "%s %s %s\n", "TENSORS", name.c_str(), "float");
	for (int i = 0; i < val.size(); i += 6)
		fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
			val[i    ], val[i + 3], val[i + 5],
			val[i + 3], val[i + 1], val[i + 4],
			val[i + 5], val[i + 4], val[i + 2]);
}

void VTKPlotFile::WriteMat3FDData(std::vector<float>& val, const std::string& name)
{
	fprintf(m_fp, "%s %s %s\n", "TENSORS", name.c_str(), "float");
	for (int i = 0; i < val.size(); i += 3)
		fprintf(m_fp, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
			val[i], 0.f, 0.f,
			0.f, val[i + 1], 0.f,
			0.f, 0.f, val[i + 2]);
}

static void Space2_(string& s)
{
	int n = (int)s.size();
	for (int i = 0; i < n; ++i)
		if (s[i] == ' ') s[i] = '_';
}

void VTKPlotFile::WriteArrayData(std::vector<float>& val, const std::string& name, FEPlotData* pd)
{
	int arraySize = pd->GetArraysize();
	fprintf(m_fp, "FIELD %s %d\n", name.c_str(), arraySize);
	std::vector<string> arrayNames = pd->GetArrayNames();
	int NE = val.size() / arraySize;
	for (int j = 0; j < arraySize; ++j)
	{
		string name = arrayNames[j];
		Space2_(name);
		fprintf(m_fp, "%s %d %d float\n", name.c_str(), 1, NE);
		for (int i = 0; i < NE; ++i)
		{
			float f = val[arraySize*i + j];
			fprintf(m_fp, "%g\n", f);
		}
	}
}

void VTKPlotFile::WriteArrayVec3fData(std::vector<float>& val, const std::string& name, FEPlotData* pd)
{
	int arraySize = pd->GetArraysize();
	fprintf(m_fp, "FIELD %s %d\n", name.c_str(), arraySize);
	std::vector<string> arrayNames = pd->GetArrayNames();
	int NE = val.size() / (3*arraySize);
	for (int j = 0; j < arraySize; ++j)
	{
		string name = arrayNames[j];
		Space2_(name);
		fprintf(m_fp, "%s %d %d float\n", name.c_str(), 3, NE);
		float f[3];
		for (int i = 0; i < NE; ++i)
		{
			f[0] = val[3*arraySize * i + 3*j    ];
			f[1] = val[3*arraySize * i + 3*j + 1];
			f[2] = val[3*arraySize * i + 3*j + 2];
			fprintf(m_fp, "%g %g %g\n", f[0], f[1], f[2]);
		}
	}
}

void VTKPlotFile::WritePointData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nodes = mesh.Nodes();

	// we count nodal variables and element variables that use NODE format
	PlotFile::Dictionary& dic = GetDictionary();
	int nodalVars = dic.NodalVariables();
	int domainVars = 0;
	auto& domainData = dic.DomainVariableList();
	list<DICTIONARY_ITEM>::iterator it = domainData.begin();
	for (int i = 0; i < domainData.size(); ++i, ++it)
	{
		if (it->m_psave)
		{
			FEPlotData* pd = it->m_psave;
			if (pd->StorageFormat() == Storage_Fmt::FMT_NODE) domainVars++;
		}
	}
	if ((nodalVars + domainVars) == 0) return;

	fprintf(m_fp, "\nPOINT_DATA %d\n", nodes);
	auto& nodeData = dic.NodalVariableList();
	it = nodeData.begin();
	for (int n = 0; n < nodeData.size(); ++n, ++it)
	{
		if (it->m_psave)
		{
			FEPlotData* pd = it->m_psave;
			int ndata = pd->VarSize(pd->DataType());

			int N = fem.GetMesh().Nodes();
			FEDataStream a; a.reserve(ndata * N);
			if (pd->Save(fem.GetMesh(), a))
			{
				// pad mismatches
				assert(a.size() == N * ndata);
				if (a.size() != N * ndata) a.resize(N * ndata, 0.f);

				// must remove all whitespace
				string dataName = it->m_szname;
				for (size_t i = 0; i < dataName.size(); ++i)
					if (isspace(dataName[i])) dataName[i] = '_';
				const char* szname = dataName.c_str();

				// write the value array
				std::vector<float>& val = a.data();
				switch (pd->DataType())
				{
				case PLT_FLOAT : WriteScalarData(val, szname); break;
				case PLT_VEC3F : WriteVectorData(val, szname); break;
				case PLT_MAT3FS: WriteMat3FSData(val, szname); break;
				case PLT_MAT3FD: WriteMat3FDData(val, szname); break;
				default:
					assert(false);
				}
			}
		}
	}

	// export all domain data that uses NODE storage format
	it = domainData.begin();
	for (int i = 0; i < domainData.size(); ++i, ++it)
	{
		if (it->m_psave)
		{
			FEPlotData* pd = it->m_psave;
			if (pd->StorageFormat() == Storage_Fmt::FMT_NODE)
			{
				// must remove all whitespace
				string dataName = it->m_szname;
				Space2_(dataName);
				const char* szname = dataName.c_str();

				int ndata = pd->VarSize(pd->DataType());

				// For now, we store all data in a global array
				int N = fem.GetMesh().Nodes();
				std::vector<float> val(ndata * N, 0.f);

				// loop over all domains and fill global val array
				for (int i = 0; i < mesh.Domains(); ++i)
				{
					FEDomain& dom = mesh.Domain(i);
					int NN = dom.Nodes();
					FEDataStream a; a.reserve(ndata * NN);
					pd->Save(dom, a);

					// pad mismatches
					if (a.size() != NN * ndata) a.resize(NN * ndata, 0.f);

					// copy to global array
					for (int j = 0; j < NN; ++j)
					{
						int nj = dom.NodeIndex(j);
						for (int k = 0; k < ndata; ++k) val[nj * ndata + k] = a[ndata * j + k];
					}
				}

				// write the value array
				switch (pd->DataType())
				{
				case PLT_FLOAT : WriteScalarData(val, szname); break;
				case PLT_VEC3F : WriteVectorData(val, szname); break;
				case PLT_MAT3FS: WriteMat3FSData(val, szname); break;
				case PLT_MAT3FD: WriteMat3FDData(val, szname); break;
				case PLT_ARRAY : WriteArrayData (val, szname, pd); break;
				case PLT_ARRAY_VEC3F: WriteArrayVec3fData(val, szname, pd); break;
				default:
					assert(false);
				}
			}
		}
	}
}

void VTKPlotFile::WriteCellData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int totalElements = mesh.Elements();

	PlotFile::Dictionary& dic = GetDictionary();
	if (dic.DomainVariables() == 0) return;

	// write cell data
	fprintf(m_fp, "\nCELL_DATA %d\n", totalElements);

	// write the part IDs first
	fprintf(m_fp, "SCALARS part_id int\n");
	fprintf(m_fp, "LOOKUP_TABLE default\n");
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int n = 0; n < NE; ++n) fprintf(m_fp, "%d\n", i);
	}

	auto& elemData = dic.DomainVariableList();
	list<DICTIONARY_ITEM>::iterator it = elemData.begin();
	for (int n = 0; n < elemData.size(); ++n, ++it)
	{
		if (it->m_psave)
		{
			FEPlotData* pd = it->m_psave;

			// For now, we can only store FE_REGION_DOMAIN/FMT_ITEM
			int nregion = pd->RegionType();
			int nformat = pd->StorageFormat();
			if ((nregion == FE_REGION_DOMAIN) && (nformat == FMT_ITEM))
			{
				// get the number of floats per data value
				int ndata = pd->VarSize(pd->DataType());

				// For now, we store all data in a global array
				std::vector<float> val(ndata * totalElements, 0.f);

				// loop over all domains and fill global val array
				int nc = 0;
				for (int i = 0; i < mesh.Domains(); ++i)
				{
					FEDomain& dom = mesh.Domain(i);
					int NE = dom.Elements();
					FEDataStream a; a.reserve(ndata * NE);
					pd->Save(dom, a);

					// pad mismatches
					if (a.size() != NE * ndata) a.resize(NE * ndata, 0.f);

					// copy into global array
					vector<float>& vi = a.data();
					for (int iel = 0; iel < NE; ++iel)
					{
						for (int k = 0; k < ndata; ++k)
						{
							val[nc++] = vi[iel * ndata + k];
						}
					}
				}

				// must remove all whitespace
				string dataName = it->m_szname;
				for (size_t i = 0; i < dataName.size(); ++i)
					if (isspace(dataName[i])) dataName[i] = '_';
				const char* szname = dataName.c_str();

				// write the value array
				switch (pd->DataType())
				{
				case PLT_FLOAT : WriteScalarData(val, szname); break;
				case PLT_VEC3F : WriteVectorData(val, szname); break;
				case PLT_MAT3FS: WriteMat3FSData(val, szname); break;
				case PLT_MAT3FD: WriteMat3FDData(val, szname); break;
				case PLT_ARRAY : WriteArrayData (val, szname, pd); break;
				case PLT_ARRAY_VEC3F: WriteArrayVec3fData(val, szname, pd); break;
				default:
					assert(false);
				}
			}
		}
	}
}
