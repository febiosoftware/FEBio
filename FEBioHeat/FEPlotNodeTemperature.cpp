#include "FEPlotNodeTemperature.h"

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
