#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"
#include "FEBodyLoad.h"
#include "FEPlotData.h"
#include "FESurface.h"
#include "FEVectorGenerator.h"
#include "FEMaterialPointMember.h"
#include "writeplot.h"

//-----------------------------------------------------------------------------
FEPlotParameter::FEPlotParameter(FEModel* pfem) : FEPlotData(pfem)
{
	m_index = 0; 
	SetStorageFormat(FMT_MULT);

	m_mat = 0;
	m_dom = 0;
	m_vec = 0;
}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotParameter::SetFilter(const char* sz)
{
	// find the parameter
	ParamString ps(sz);
	m_param = GetFEModel()->GetParameterValue(ps);
	if (m_param.isValid() == false)
	{
		// the string does not reference a parameter, let's see if it is a value property
		FECoreBase* pp = GetFEModel()->FindComponent(ps);
		if (pp == nullptr) return false;

		m_vec = dynamic_cast<FEVectorGenerator*>(pp);
		if (m_vec == 0) return false;

		m_mat = dynamic_cast<FEMaterial*>(m_vec->GetAncestor());
		if (m_mat == 0) return false;

		SetRegionType(FE_REGION_DOMAIN);
		SetVarType(PLT_VEC3F);

		return true;
	}

	FEParam* param = m_param.param();
	if (param == 0) return false;

	FECoreBase* pc = dynamic_cast<FECoreBase*>(param->parent());
	if (pc == nullptr) return false;

	switch (m_param.type())
	{
	case FE_PARAM_DOUBLE_MAPPED:
	{
		FEParamDouble& p = m_param.value<FEParamDouble>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if (     dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
			else if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
			else return false;
		}
		else
		{
			if (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) SetRegionType(FE_REGION_SURFACE);
			else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
			else return false;
		}

		SetVarType(PLT_FLOAT);
	}
	break;
	case FE_PARAM_VEC3D_MAPPED:
	{
		FEParamVec3& p = m_param.value<FEParamVec3>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if (     dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
			else if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
			else return false;
		}
		else
		{
			if (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) SetRegionType(FE_REGION_SURFACE);
			else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
			else return false;
		}

		SetVarType(PLT_VEC3F);
	}
	break;
	case FE_PARAM_DOUBLE:
	{
		FEParamContainer* pc = param->parent();
		if (pc == 0) return false;

		if (dynamic_cast<FEMaterial*>(pc))
		{
			m_mat = dynamic_cast<FEMaterial*>(pc);
			SetRegionType(FE_REGION_DOMAIN);
			SetStorageFormat(FMT_REGION);
		}
		else return false;
	}
	break;
	case FE_PARAM_MATERIALPOINT:
	{
		FEMaterialPointMember& prop = param->value<FEMaterialPointMember>();
		m_mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
		if (m_mat == nullptr) return false;

		SetRegionType(FE_REGION_DOMAIN);

		switch (prop.dataType())
		{
		case FE_DOUBLE: SetVarType(PLT_FLOAT); break;
		case FE_VEC3D : SetVarType(PLT_VEC3F); break;
		case FE_MAT3D : SetVarType(PLT_MAT3F); break;
		default:
			return false;
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
	if (m_vec)
	{
		if (m_mat == 0) return false;
		if (dom.GetMaterial() != m_mat) return false;

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		writeNodalProjectedElementValues<vec3d>(sd, a, [=](const FEMaterialPoint& mp) {
			return m_vec->GetVector(mp);
		});

		return true;
	}

	if (m_param.isValid() == false) return false;

	if ((m_param.type() == FE_PARAM_DOUBLE_MAPPED) ||
		(m_param.type() == FE_PARAM_VEC3D_MAPPED))
	{
		FEModelParam& map = m_param.value<FEModelParam>();

		if (m_dom == nullptr)
		{
			if ((m_mat == nullptr) || (dom.GetMaterial() != m_mat))
				return false;
		}
		else if (m_dom->IsMember(&dom) == false) return false;

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		if (m_param.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& mapDouble = dynamic_cast<FEParamDouble&>(map);
			writeNodalProjectedElementValues<double>(sd, a, mapDouble);
		}
		else if (m_param.type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& mapVec3 = dynamic_cast<FEParamVec3&>(map);
			writeNodalProjectedElementValues<vec3d>(sd, a, mapVec3);
		}
	
		return true;
	}
	else if (m_param.type() == FE_PARAM_DOUBLE)
	{
		if (dom.GetMaterial() == m_mat)
		{
			double val = m_param.value<double>();
			a << val;
			return true;
		}
	}
	else if (m_param.type() == FE_PARAM_MATERIALPOINT)
	{
		if (dom.GetMaterial() == m_mat)
		{
			FEMaterialPointMember& prop = m_param.value<FEMaterialPointMember>();

			switch (prop.dataType())
			{
			case FE_DOUBLE: writeNodalProjectedElementValues<double>(dom, a, [&](const FEMaterialPoint& mp) {
					FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
					double d;
					prop.get(mp_noconst, d);
					return d;
				});	
				break;
			case FE_VEC3D: writeNodalProjectedElementValues<vec3d>(dom, a, [&](const FEMaterialPoint& mp) {
					FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
					vec3d d;
					prop.get(mp_noconst, d);
					return d;
				});	
				break;
			case FE_MAT3D: writeNodalProjectedElementValues<mat3d>(dom, a, [&](const FEMaterialPoint& mp) {
					FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
					mat3d d;
					prop.get(mp_noconst, d);
					return d;
				});	
				break;

			default:
				return false;
			}
			return true;
		}
	}

	return false;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotParameter::Save(FESurface& dom, FEDataStream& a)
{
	if (m_param.isValid() == false) return false;

	if (m_param.type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = m_param.value<FEParamDouble>();
		FEFacetSet* surf = dynamic_cast<FEFacetSet*>(map.GetItemList());
		if (surf != dom.GetFacetSet()) return false;

		writeNodalProjectedElementValues<double>(dom, a, map);
	}
	else if (m_param.type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = m_param.value<FEParamVec3>();

		FEFacetSet* surf = dynamic_cast<FEFacetSet*>(map.GetItemList());
		if (surf != dom.GetFacetSet()) return false;

		writeNodalProjectedElementValues<vec3d>(dom, a, map);
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotParameter::Save(FEMesh& mesh, FEDataStream& a)
{
	if (m_param.isValid() == false) return false;

	if (m_param.type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = m_param.value<FEParamDouble>();
		FENodeSet* nset = dynamic_cast<FENodeSet*>(map.GetItemList());
		if (nset == 0) return false;

		// write the nodal values
		writeNodalValues<double>(*nset, a, map);

		return true;
	}

	return false;
}
