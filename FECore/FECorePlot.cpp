#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"

//-----------------------------------------------------------------------------
FEPlotMaterialParameter::FEPlotMaterialParameter(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_MULT) { m_index = 0; }

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotMaterialParameter::SetFilter(const char* sz)
{
	// copy name into temp buffer
	char szbuf[256] = {0};
	strcpy(szbuf, sz);

	// strip the material name
	char* ch = strchr(szbuf, '.');
	if (ch) *ch++ = 0; else return false;

	// get the last component
	char* cd = strrchr(ch, '.');
	if (cd == 0) cd = ch; else cd++;

	// get the optional index
	char* chl = strchr(cd, '[');
	if (chl)
	{
		*chl++ = 0;
		char* chr = strrchr(chl, ']');
		if (chr == 0) return false;
		*chr=0;
		m_index = atoi(chl);
		if (m_index < 0) return false;
	}

	strcpy(m_szmat, szbuf);
	strcpy(m_szparam, ch);

	return true;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
// TODO: use ParamString
bool FEPlotMaterialParameter::Save(FEDomain& dom, FEDataStream& a)
{
	// First, get the domain material
	FEMaterial* pmat = dom.GetMaterial();
	if (pmat==0) return false;

	// check the name of this material
	if (strcmp(pmat->GetName(), m_szmat) != 0) return false;

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

			// extract the parameter
			// Note that for now this only works for double parameters
			FEParam* pv = mp.FindParameter(m_szparam);
			if (pv && (pv->type()==FE_PARAM_DOUBLE) && (m_index < pv->dim()))
			{
				gv[j] = pv->value<double>(m_index);
				nc++;
			}
		}

		vector<double> nv(neln, 0.0);
		if (nc == nint)
		{
			e.project_to_nodes(&gv[0], &nv[0]);
		}

		// store the result
		for (int j=0; j<neln; ++j) a << nv[j];
	}

	return true;
}
