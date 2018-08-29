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
	std::string		m_matName;		//!< material name
	std::string		m_paramName;	//!< parameter name
	int				m_index;		//!< index for array parameters
};
