// PlotFile.cpp: implementation of the PlotFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "PlotFile.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlotFile::PlotFile()
{
}

PlotFile::~PlotFile()
{
	Close();
}

//-----------------------------------------------------------------------------
void PlotFile::Close()
{
	m_ar.Close();
}

//-----------------------------------------------------------------------------
void PlotFile::write_displacements()
{
	FEM& fem = *m_pfem;

	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_rt.x;
		xf[1] = (float) node.m_rt.y;
		xf[2] = (float) node.m_rt.z;

		m_ar.write(xf, sizeof(float), 3);
	}
}
