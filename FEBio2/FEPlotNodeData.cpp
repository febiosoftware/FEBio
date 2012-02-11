#include "stdafx.h"
#include "FEPlotNodeData.h"
#include "FECore/FEMesh.h"
#include "FECore/febio.h"
//-----------------------------------------------------------------------------
REGISTER_FEBIO_CLASS(FEPlotNodeDisplacement  , FEPlotData, "displacement"   );
REGISTER_FEBIO_CLASS(FEPlotNodeVelocity      , FEPlotData, "velocity"       );
REGISTER_FEBIO_CLASS(FEPlotNodeAcceleration  , FEPlotData, "acceleration"   );
REGISTER_FEBIO_CLASS(FEPlotNodeTemperature   , FEPlotData, "temperature"    );
REGISTER_FEBIO_CLASS(FEPlotNodeReactionForces, FEPlotData, "reaction forces");
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

//-----------------------------------------------------------------------------
//! Store nodal reaction forces
bool FEPlotNodeReactionForces::Save(FEMesh& m, vector<float>& a)
{
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		FENode& node = m.Node(i);
		a.push_back((float) node.m_Fr.x);
		a.push_back((float) node.m_Fr.y);
		a.push_back((float) node.m_Fr.z);
	}
	return true;
}
