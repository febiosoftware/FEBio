#pragma once
#include "FEPlotData.h"

//-----------------------------------------------------------------------------
// class for exporting element specific material parameters to plot file
class FEPlotParameter : public FEPlotDomainData
{
public:
	FEPlotParameter(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

	virtual bool SetFilter(const char* sz);

protected:
	FEModel*		m_fem;		//!< TODO: Base class also has a model parameter, but it doesn't look like it's set
	FEParam*		m_param;	//!< parameter name
	int				m_index;	//!< index for array parameters
};
