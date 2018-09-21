#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
FEPlotParameter::FEPlotParameter(FEModel* pfem)
{
	m_fem = pfem;
	m_param = 0;
	m_index = 0; 
	SetStorageFormat(FMT_MULT);
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

	switch (m_param->type())
	{
	case FE_PARAM_DOUBLE_MAPPED:
	{
		FEParamDouble& p = m_param->value<FEParamDouble>();
		FEDomainList& domList = p.getDomainList();

		// we assume for now that all domains in this list are of the same type.
		if (domList.IsEmpty()) return false;
		FEDomain* dom = domList.GetDomain(0);

		// If the domain is not set, we assume it's a element domain parameter
		if (dom == 0) SetRegionType(FE_REGION_DOMAIN);
		else
		{
			if (dynamic_cast<FESolidDomain*>(dom)) SetRegionType(FE_REGION_DOMAIN);
			else if (dynamic_cast<FESurface*>(dom)) SetRegionType(FE_REGION_SURFACE);
			else return false;
		}
		SetVarType(PLT_FLOAT);
	}
	break;
	case FE_PARAM_VEC3D_MAPPED:
	{
		FEParamVec3& p = m_param->value<FEParamVec3>();
		FEDomainList& domList = p.getDomainList();

		// we assume for now that all domains in this list are of the same type.
		if (domList.IsEmpty()) return false;
		FEDomain* dom = domList.GetDomain(0);

		// If the domain is not set, we assume it's a element domain parameter
		if (dom == 0) SetRegionType(FE_REGION_DOMAIN);
		else
		{
			if (dynamic_cast<FESolidDomain*>(dom)) SetRegionType(FE_REGION_DOMAIN);
			else if (dynamic_cast<FESurface*>(dom)) SetRegionType(FE_REGION_SURFACE);
			else return false;
		}
		SetVarType(PLT_VEC3F);
	}
	break;
	default:
		assert(false);
		return false;
		break;
	}

	return true;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotParameter::Save(FEDomain& dom, FEDataStream& a)
{
	if (m_param == 0) return false;
	FEParam* param = m_param;

	if (param->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = param->value<FEParamDouble>();
		if (map.getDomainList().IsMember(&dom) == false) return false;
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		double gi[FEElement::MAX_INTPOINTS] = { 0 };
		double gn[FEElement::MAX_NODES] = { 0 };

		// loop over all the elements in the domain
		int NE = dom.Elements();
		for (int i = 0; i < NE; ++i)
		{
			// get the element and loop over its integration points
			// we only calculate the element's average
			// but since most material parameters can only defined 
			// at the element level, this should get the same answer
			FESolidElement& e = sd.Element(i);
			int nint = e.GaussPoints();
			int neln = e.Nodes();

			for (int j = 0; j < nint; ++j)
			{
				// get the material point data for this integration point
				FEMaterialPoint& mp = *e.GetMaterialPoint(j);
				gi[j] = map.eval(mp);
			}
			
			e.project_to_nodes(gi, gn);

			// store the result
			for (int j=0; j<neln; ++j) a << gn[j];
		}
	}
	else if (param->type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = param->value<FEParamVec3>();
		if (map.getDomainList().IsMember(&dom) == false) return false;
		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		vec3d gi[FEElement::MAX_INTPOINTS];
		vec3d gn[FEElement::MAX_NODES];

		// loop over all the elements in the domain
		int NE = dom.Elements();
		for (int i = 0; i < NE; ++i)
		{
			// get the element and loop over its integration points
			// we only calculate the element's average
			// but since most material parameters can only defined 
			// at the element level, this should get the same answer
			FESolidElement& e = sd.Element(i);
			int nint = e.GaussPoints();
			int neln = e.Nodes();

			for (int j = 0; j < nint; ++j)
			{
				// get the material point data for this integration point
				FEMaterialPoint& mp = *e.GetMaterialPoint(j);
				gi[j] = map.eval(mp);
			}

			e.project_to_nodes(gi, gn);

			// store the result
			for (int j = 0; j<neln; ++j) a << gn[j];
		}
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotParameter::Save(FESurface& dom, FEDataStream& a)
{
	if (m_param == 0) return false;
	FEParam* param = m_param;

	if (param->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = param->value<FEParamDouble>();
		if (map.getDomainList().IsMember(&dom) == false) return false;

		double gi[FEElement::MAX_INTPOINTS] = { 0 };
		double gn[FEElement::MAX_NODES] = { 0 };

		// loop over all the elements in the domain
		int NE = dom.Elements();
		for (int i = 0; i < NE; ++i)
		{
			// get the element and loop over its integration points
			// we only calculate the element's average
			// but since most material parameters can only defined 
			// at the element level, this should get the same answer
			FESurfaceElement& e = dom.Element(i);
			int nint = e.GaussPoints();
			int neln = e.Nodes();

			for (int j = 0; j < nint; ++j)
			{
				// get the material point data for this integration point
				FEMaterialPoint& mp = *e.GetMaterialPoint(j);
				gi[j] = map.eval(mp);
			}

			e.project_to_nodes(gi, gn);

			// store the result
			// NOTE: Note that we always need to store 10 entries. This is because of a limitation of the plot file format.
			for (int j = 0; j < 10; ++j) a << gn[j];
		}
	}
	else if (param->type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = param->value<FEParamVec3>();
		if (map.getDomainList().IsMember(&dom) == false) return false;

		vec3d gi[FEElement::MAX_INTPOINTS];
		vec3d gn[FEElement::MAX_NODES];

		// loop over all the elements in the domain
		int NE = dom.Elements();
		for (int i = 0; i < NE; ++i)
		{
			// get the element and loop over its integration points
			// we only calculate the element's average
			// but since most material parameters can only defined 
			// at the element level, this should get the same answer
			FESurfaceElement& e = dom.Element(i);
			int nint = e.GaussPoints();
			int neln = e.Nodes();

			vec3d v(0, 0, 0);
			for (int j = 0; j < nint; ++j)
			{
				// get the material point data for this integration point
				FEMaterialPoint& mp = *e.GetMaterialPoint(j);
				gi[j] = map.eval(mp);
			}

			e.project_to_nodes(gi, gn);

			// store the result
			// NOTE: Note that we always need to store 10 entries. This is because of a limitation of the plot file format.
			for (int j = 0; j < 10; ++j) a << gn[j];
		}
	}
	else return false;

	return true;
}
