#include "FEPlotNodeTemperature.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//! Store the nodal temperatures
bool FEPlotNodeTemperature::Save(FEMesh& m, FEDataStream& a)
{
	// get the temperature dof index
	int dof_t = GetFEModel()->GetDOFIndex("t");
	if (dof_t == -1) return false;

	// store the temperatures
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.get(dof_t);
	}
	return true;
}
