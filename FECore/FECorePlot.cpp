#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
FEPlotParameter::FEPlotParameter(FEModel* pfem) : FEPlotDomainData(PLT_FLOAT, FMT_MULT) 
{
	m_fem = pfem;
	m_param = 0;
	m_index = 0; 
}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotParameter::SetFilter(const char* sz)
{
	assert(m_fem);
	if (m_fem == 0) return false;

	// find the parameter
	ParamString ps(sz);
	m_param = m_fem->FindParameter(ps);
	if (m_param == 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotParameter::Save(FEDomain& dom, FEDataStream& a)
{
	if (m_param == 0) return false;
	FEParam* param = m_param;

	if (param->type() != FE_PARAM_DOUBLE_MAPPED) return false;
	FEModelParam& map = param->value<FEModelParam>();

	FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

	// loop over all the elements in the domain
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the element and loop over its integration points
		// we only calculate the element's average
		// but since most material parameters can only defined 
		// at the element level, this should get the same answer
		FESolidElement& e = sd.Element(i);
		int nint = e.GaussPoints();
		int neln = e.Nodes();

		vector<double> gv(nint);
		double E = 0.0;
		int nc = 0;
		for (int j=0; j<nint; ++j)
		{
			// get the material point data for this integration point
			FEMaterialPoint& mp = *e.GetMaterialPoint(j);

			gv[j] = map.eval(mp);
		}

		vector<double> nv(neln, 0.0);
		e.project_to_nodes(&gv[0], &nv[0]);

		// store the result
		for (int j=0; j<neln; ++j) a << nv[j];
	}

	return true;
}
