#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
FEPlotMaterialParameter::FEPlotMaterialParameter(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotMaterialParameter::SetFilter(const char* sz)
{
	char szbuf[256] = {0};
	strcpy(szbuf, sz);
	char* ch = strchr(szbuf, '.');
	if (ch) *ch++ = 0; else return false;

	strcpy(m_szmat, szbuf);
	strcpy(m_szparam, ch);

	return true;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotMaterialParameter::Save(FEDomain& dom, FEPlotStream& a)
{
	// First, get the domain material
	FEMaterial* pmat = dom.GetMaterial();
	if (pmat==0) return false;

	// check the name of this material
	if (strcmp(pmat->GetName(), m_szmat) != 0) return false;

	// setup the parameter string which is a way
	// to reference a specific parameter
	ParamString s(m_szparam);

	// loop over all the elements in the domain
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the element and loop over its integration points
		// we only calculate the element's average
		// but since most material parameters can only defined 
		// at the element level, this should get the same answer
		FEElement& e = dom.ElementRef(i);
		int nint = e.GaussPoints();
		double E = 0.0;
		int nc = 0;
		for (int j=0; j<nint; ++j)
		{
			// get the material point data for this integration point
			FEMaterialPoint& mp = *e.GetMaterialPoint(j);

			// extract the parameter
			// Note that for now this only works for double parameters
			FEParam* pv = mp.GetParameter(s);
			if (pv && (pv->m_itype==FE_PARAM_DOUBLE))
			{
				E += pv->value<double>();
				nc++;
			}
		}
		if (nc > 0) E /= (double) nc;

		// store the result
		a << E;
	}

	return true;
}
