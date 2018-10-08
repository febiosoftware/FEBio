#pragma once
#include "FEPlotData.h"

//-----------------------------------------------------------------------------
// class for exporting element specific material parameters to plot file
class FEPlotParameter : public FEPlotData
{
public:
	FEPlotParameter(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool Save(FESurface& dom, FEDataStream& a);
	bool Save(FEMesh& mesh, FEDataStream& a);

	virtual bool SetFilter(const char* sz);

protected:
	FEParam*		m_param;	//!< parameter name
	int				m_index;	//!< index for array parameters
};
