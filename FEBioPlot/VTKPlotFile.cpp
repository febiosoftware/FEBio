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

}

//! Open the plot database
bool VTKPlotFile::Open(const char* szfile)
{
	return true;
}

//! Open for appending
bool VTKPlotFile::Append(const char* szfile)
{
	return true;
}

//! see if the plot file is valid
bool VTKPlotFile::IsValid() const
{
	return true;
}

//! Write current FE state to plot database
bool VTKPlotFile::Write(float ftime, int flag)
{
	static int n = 1;
	FEModel& fem = *GetFEModel();

	std::stringstream ss;
	ss << "out" << n++ << ".vtk";
	string fileName = ss.str();
	
	m_fp = fopen(fileName.c_str(), "wt");
	if (m_fp == nullptr) return false;

	// --- H E A D E R ---
	WriteHeader();

	// --- N O D E S ---
	WritePoints();

	// --- E L E M E N T S ---
	WriteCells();

	fclose(m_fp);
	m_fp = nullptr;

	return true;
}

//-----------------------------------------------------------------------------
void VTKPlotFile::WriteHeader()
{
	FEModel& fem = *GetFEModel();
	fprintf(m_fp, "%s\n", "# vtk DataFile Version 3.0");
	fprintf(m_fp, "%s %g\n", "vtk output at time", fem.GetCurrentTime());
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
			vec3d& r = m.Node(j + k).m_rt;
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
