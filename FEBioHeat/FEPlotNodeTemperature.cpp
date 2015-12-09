#include "FEPlotNodeTemperature.h"

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotNodeTemperature::Save(FEMesh& m, FEDataStream& a)
{
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.m_T;
	}
	return true;
}
