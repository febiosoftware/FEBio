#pragma once
#include "FEPlotData.h"

//-----------------------------------------------------------------------------
// class for exporting element specific material parameters to plot file
class FEPlotMaterialParameter : public FEDomainData
{
public:
	FEPlotMaterialParameter(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

	virtual bool SetFilter(const char* sz);

protected:
	char	m_szmat[128];	// material name
	char	m_szparam[128];	// parameter name
	int		m_index;		// index for array parameters
};
