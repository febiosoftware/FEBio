#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"
#include "FEBodyLoad.h"
#include "FEPlotData.h"
#include "FESurface.h"

//-----------------------------------------------------------------------------
FEPlotParameter::FEPlotParameter(FEModel* pfem) : FEPlotData(pfem)
{
	m_param = 0;
	m_index = 0; 
	SetStorageFormat(FMT_MULT);
}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotParameter::SetFilter(const char* sz)
{
	// find the parameter
	ParamString ps(sz);
	m_param = GetFEModel()->FindParameter(ps);
	if (m_param == 0) return false;

	switch (m_param->type())
	{
	case FE_PARAM_DOUBLE_MAPPED:
	{
		FEParamDouble& p = m_param->value<FEParamDouble>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if (dynamic_cast<FEMaterial*>(m_param->parent())) SetRegionType(FE_REGION_DOMAIN);
			else if (dynamic_cast<FEBodyLoad*>(m_param->parent())) SetRegionType(FE_REGION_DOMAIN);
			else return false;
		}
		else
		{
			if (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) SetRegionType(FE_REGION_SURFACE);
			else if (dynamic_cast<FEElementSet*>(itemList)) SetRegionType(FE_REGION_DOMAIN);
			else return false;
		}

		SetVarType(PLT_FLOAT);
	}
	break;
	case FE_PARAM_VEC3D_MAPPED:
	{
		FEParamVec3& p = m_param->value<FEParamVec3>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if (dynamic_cast<FEMaterial*>(m_param->parent())) SetRegionType(FE_REGION_DOMAIN);
			else if (dynamic_cast<FEBodyLoad*>(m_param->parent())) SetRegionType(FE_REGION_DOMAIN);
			else return false;
		}
		else
		{
			if (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) SetRegionType(FE_REGION_SURFACE);
			else if (dynamic_cast<FEElementSet*>(itemList)) SetRegionType(FE_REGION_DOMAIN);
			else return false;
		}

		SetVarType(PLT_VEC3F);
	}
	break;
	case FE_PARAM_DOUBLE:
	{
		FEParamContainer* pc = m_param->parent();
		if (pc == 0) return false;

		if (dynamic_cast<FEMaterial*>(pc))
		{
			SetRegionType(FE_REGION_DOMAIN);
			SetStorageFormat(FMT_REGION);
		}
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

	if ((param->type() == FE_PARAM_DOUBLE_MAPPED) ||
		(param->type() == FE_PARAM_VEC3D_MAPPED))
	{
		FEModelParam& map = param->value<FEModelParam>();

		FEDomainList* domList = 0;

		FEElementSet* elset = dynamic_cast<FEElementSet*>(map.GetItemList());
		if (elset == 0)
		{
			FEMaterial* mat = dynamic_cast<FEMaterial*>(param->parent());
			if (mat) domList = &mat->GetDomainList();
			else
			{
				FEBodyLoad* bl = dynamic_cast<FEBodyLoad*>(param->parent());
				if (bl) domList = &bl->GetDomainList();
				else return false;
			}
		}
		else domList = &elset->GetDomainList();
		if (domList->IsMember(&dom) == false) return false;

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		if (param->type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& mapDouble = dynamic_cast<FEParamDouble&>(map);
			writeNodalProjectedElementValues(sd, a, mapDouble);
		}
		else if (param->type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& mapVec3 = dynamic_cast<FEParamVec3&>(map);
			writeNodalProjectedElementValues(sd, a, mapVec3);
		}
	
		return true;
	}
	else if (param->type() == FE_PARAM_DOUBLE)
	{
		FEParamContainer* pc = param->parent();
		if (dynamic_cast<FEMaterial*>(pc))
		{
			FEMaterial* mat = dynamic_cast<FEMaterial*>(pc);
			if (dom.GetMaterial() == mat)
			{
				double val = param->value<double>();
				a << val;
				return true;
			}
		}
	}

	return false;
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
		FEFacetSet* surf = dynamic_cast<FEFacetSet*>(map.GetItemList());
		if (surf != dom.GetFacetSet()) return false;

		writeNodalProjectedElementValues(dom, a, map);
	}
	else if (param->type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = param->value<FEParamVec3>();

		FEFacetSet* surf = dynamic_cast<FEFacetSet*>(map.GetItemList());
		if (surf != dom.GetFacetSet()) return false;

		writeNodalProjectedElementValues(dom, a, map);
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotParameter::Save(FEMesh& mesh, FEDataStream& a)
{
	if (m_param == 0) return false;
	FEParam* param = m_param;

	if (param->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = param->value<FEParamDouble>();
		FENodeSet* nset = dynamic_cast<FENodeSet*>(map.GetItemList());
		if (nset == 0) return false;

		// write the nodal values
		writeNodalValues(*nset, a, map);

		return true;
	}

	return false;
}
