#pragma once
#include "FEPlotData.h"

//-----------------------------------------------------------------------------
class FEMaterial;
class FEDomainList;
class FEVectorGenerator;
class FEFacetSet;

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
	FEParamValue	m_param;	//!< parameter
	int				m_index;	//!< index for array parameters

private:
	FEMaterial*			m_mat;
	FEDomainList*		m_dom;
	FEFacetSet*			m_surf;
	FEVectorGenerator*	m_vec;
};
