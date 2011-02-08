#include "stdafx.h"
#include "FEPlotNodeData.h"
#include "FEPlotDataFactory.h"

//-----------------------------------------------------------------------------
REGISTER_PLOTDATA(FEPlotNodeDisplacement, "displacement");
REGISTER_PLOTDATA(FEPlotNodeVelocity    , "velocity"    );
REGISTER_PLOTDATA(FEPlotNodeAcceleration, "acceleration");
REGISTER_PLOTDATA(FEPlotNodeTemperature , "temperature" );
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotNodeDisplacement::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) (node.m_rt.x - node.m_r0.x);
		xf[1] = (float) (node.m_rt.y - node.m_r0.y);
		xf[2] = (float) (node.m_rt.z - node.m_r0.z);

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeVelocity::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_vt.x;
		xf[1] = (float) node.m_vt.y;
		xf[2] = (float) node.m_vt.z;

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotNodeAcceleration::Save(FEMesh& m, vector<float>& a)
{
	float xf[3];
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_at.x;
		xf[1] = (float) node.m_at.y;
		xf[2] = (float) node.m_at.z;

		a.push_back(xf[0]);
		a.push_back(xf[1]);
		a.push_back(xf[2]);
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotNodeTemperature::Save(FEMesh& m, vector<float>& a)
{
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		float f = (float) node.m_T;
		a.push_back(f);
	}
	return true;
}
