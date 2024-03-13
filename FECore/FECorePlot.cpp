/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"
#include "FEBodyLoad.h"
#include "FENodalLoad.h"
#include "FEPlotData.h"
#include "FESurface.h"
#include "FEMaterialPointProperty.h"
#include "writeplot.h"
#include "FESurfaceLoad.h"
#include "FEDomainMap.h"
#include "FEModel.h"
#include "FEPIDController.h"
#include "FEDataMap.h"

//-----------------------------------------------------------------------------
FEPlotParameter::FEPlotParameter(FEModel* pfem) : FEPlotData(pfem)
{
	m_index = 0; 
	SetStorageFormat(FMT_MULT);

	m_mat = 0;
	m_dom = 0;
	m_surf = 0;
}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotParameter::SetFilter(const char* sz)
{
	// store the filter for serialization
	m_filter = sz;

	// find the parameter
	ParamString ps(sz);
	m_param = GetFEModel()->GetParameterValue(ps);
	if (m_param.isValid() == false) return false;

	FEParam* param = m_param.param();
	if (param == 0) return false;

	FECoreBase* pc = dynamic_cast<FECoreBase*>(param->parent());
	if (pc == nullptr) return false;

	switch (m_param.type())
	{
	case FE_PARAM_DOUBLE_MAPPED:
	{
		FEParamDouble& p = m_param.value<FEParamDouble>();

		// check for materials first
		if (dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
		else 
		{
			FEItemList* itemList = p.GetItemList();
			if (itemList == 0)
			{
				// for some classes, the item list can be empty
				if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
	//			else if (dynamic_cast<FESurfaceLoad*>(pc)) { SetRegionType(FE_REGION_SURFACE); m_surf = dynamic_cast<FESurfaceLoad*>(pc)->GetSurface().GetFacetSet(); }
				else if (dynamic_cast<FENodalLoad*>(pc)) SetRegionType(FE_REGION_NODE);
				else return false;
			}
			else
			{
				if      (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
				else if (dynamic_cast<FEFacetSet*>(itemList)) { SetRegionType(FE_REGION_SURFACE); m_surf = dynamic_cast<FEFacetSet*>(itemList); }
				else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
				else return false;
			}
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
			if      (dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
			else if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
			else return false;
		}
		else
		{
			if      (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) { SetRegionType(FE_REGION_SURFACE); m_surf = dynamic_cast<FEFacetSet*>(itemList); }
			else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
			else return false;
		}

		SetVarType(PLT_VEC3F);
	}
	break;
	case FE_PARAM_MAT3D_MAPPED:
	{
		FEParamMat3d& p = m_param.value<FEParamMat3d>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if      (dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
			else if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
			else return false;
		}
		else
		{
			if      (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) { SetRegionType(FE_REGION_SURFACE); m_surf = dynamic_cast<FEFacetSet*>(itemList); }
			else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
			else return false;
		}

		SetStorageFormat(FMT_ITEM);
		SetVarType(PLT_MAT3F);
	}
	break;
	case FE_PARAM_MAT3DS_MAPPED:
	{
		FEParamMat3ds& p = m_param.value<FEParamMat3ds>();

		FEItemList* itemList = p.GetItemList();
		if (itemList == 0)
		{
			// for material parameters, the item list can be empty
			if      (dynamic_cast<FEMaterial*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_mat = dynamic_cast<FEMaterial*>(pc); m_mat = dynamic_cast<FEMaterial*>(m_mat->GetAncestor()); }
			else if (dynamic_cast<FEBodyLoad*>(pc)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList(); }
			else return false;
		}
		else
		{
			if      (dynamic_cast<FENodeSet*>(itemList)) SetRegionType(FE_REGION_NODE);
			else if (dynamic_cast<FEFacetSet*>(itemList)) { SetRegionType(FE_REGION_SURFACE); m_surf = dynamic_cast<FEFacetSet*>(itemList); }
			else if (dynamic_cast<FEElementSet*>(itemList)) { SetRegionType(FE_REGION_DOMAIN); m_dom = &(dynamic_cast<FEElementSet*>(itemList))->GetDomainList(); }
			else return false;
		}

		SetStorageFormat(FMT_ITEM);
		SetVarType(PLT_MAT3FS);
	}
	break;
	case FE_PARAM_DOUBLE:
	{
		if (dynamic_cast<FEMaterial*>(pc))
		{
			m_mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
			SetRegionType(FE_REGION_DOMAIN);
			SetStorageFormat(FMT_REGION);
		}
		else if (dynamic_cast<FESurfaceLoad*>(pc))
		{
			FESurfaceLoad* psl = dynamic_cast<FESurfaceLoad*>(pc);
			SetRegionType(FE_REGION_SURFACE);
			SetStorageFormat(FMT_REGION);
			m_surf = psl->GetSurface().GetFacetSet();
		}
		else return false;
	}
	break;
	case FE_PARAM_VEC3D:
	{
		if (dynamic_cast<FEMaterial*>(pc))
		{
			m_mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
			SetRegionType(FE_REGION_DOMAIN);
			SetStorageFormat(FMT_REGION);
		}
		else if (dynamic_cast<FEBodyLoad*>(pc))
		{
			SetRegionType(FE_REGION_DOMAIN);
			SetStorageFormat(FMT_REGION);
			m_dom = &(dynamic_cast<FEBodyLoad*>(pc))->GetDomainList();
		}
		else if (dynamic_cast<FESurfaceLoad*>(pc))
		{
			FESurfaceLoad* psl = dynamic_cast<FESurfaceLoad*>(pc);
			SetRegionType(FE_REGION_SURFACE);
			SetStorageFormat(FMT_REGION);
			m_surf = psl->GetSurface().GetFacetSet();
		}
		else return false;

		SetVarType(PLT_VEC3F);
	}
	break;
	case FE_PARAM_MATERIALPOINT:
	{
		FEMaterialPointProperty& prop = m_param.value<FEMaterialPointProperty>();
		m_mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
		if (m_mat == nullptr) return false;

		SetRegionType(FE_REGION_DOMAIN);

		switch (prop.dataType())
		{
		case FE_DOUBLE: SetVarType(PLT_FLOAT); break;
		case FE_VEC3D: SetVarType(PLT_VEC3F); break;
		case FE_MAT3D: SetVarType(PLT_MAT3F); break;
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
void FEPlotParameter::Serialize(DumpStream& ar)
{
	FEPlotData::Serialize(ar);
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
		ar << m_filter;
	else
	{
		string filter;
		ar >> filter;
		SetFilter(filter.c_str());
	}
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotParameter::Save(FEDomain& dom, FEDataStream& a)
{
	if (m_param.isValid() == false) return false;

	FEParam* param = m_param.param();

	if ((m_param.type() == FE_PARAM_DOUBLE_MAPPED) ||
		(m_param.type() == FE_PARAM_VEC3D_MAPPED ) ||
		(m_param.type() == FE_PARAM_MAT3D_MAPPED ) ||
		(m_param.type() == FE_PARAM_MAT3DS_MAPPED))
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

			FEMappedValue* val = dynamic_cast<FEMappedValue*>(mapDouble.valuator());
			if (val)
			{
				FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap()); assert(map);
				if (map->StorageFormat() == FMT_MULT)
				{
					// loop over all elements
					int NE = dom.Elements();
					for (int i = 0; i < NE; ++i)
					{
						FEElement& e = dom.ElementRef(i);
						int ne = e.Nodes();

						vector<double> sn(ne);
						for (int j = 0; j < ne; ++j)
						{
							sn[j] = map->value<double>(i, j);
						}

						// push data to archive
						for (int j = 0; j < ne; ++j) a << sn[j];
					}

					return true;
				}
			}

			writeNodalProjectedElementValues<double>(sd, a, mapDouble);
		}
		else if (m_param.type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& mapVec3 = dynamic_cast<FEParamVec3&>(map);
			writeNodalProjectedElementValues<vec3d>(sd, a, mapVec3);
		}
		else if (m_param.type() == FE_PARAM_MAT3D_MAPPED)
		{
			FEParamMat3d& mapMat3 = dynamic_cast<FEParamMat3d&>(map);
			writeElementValue<mat3d>(sd, a, mapMat3);
		}
		else if (m_param.type() == FE_PARAM_MAT3DS_MAPPED)
		{
			FEParamMat3ds& mapMat3 = dynamic_cast<FEParamMat3ds&>(map);
			writeElementValue<mat3ds>(sd, a, mapMat3);
		}
		else return false;

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
	else if (m_param.type() == FE_PARAM_VEC3D)
	{
		if (m_dom && (m_dom->IsMember(&dom)))
		{
			vec3d val = m_param.value<vec3d>();
			a << val;
			return true;
		}
		else if (dom.GetMaterial() == m_mat)
		{
			vec3d val = m_param.value<vec3d>();
			a << val;
			return true;
		}
		else return false;
	}
	else if (m_param.type() == FE_PARAM_MATERIALPOINT)
	{
		if (dom.GetMaterial() == m_mat)
		{
			FEMaterialPointProperty& prop = m_param.value<FEMaterialPointProperty>();

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

	FEFacetSet* surf = dom.GetFacetSet();
	if (m_surf != surf) return false;

	if (m_param.type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = m_param.value<FEParamDouble>();
		writeNodalProjectedElementValues<double>(dom, a, map);
	}
	else if (m_param.type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = m_param.value<FEParamVec3>();
		writeNodalProjectedElementValues<vec3d>(dom, a, map);
	}
	else if (m_param.type() == FE_PARAM_DOUBLE)
	{
		a << m_param.value<double>();
	}
	else if (m_param.type() == FE_PARAM_VEC3D)
	{
		a << m_param.value<vec3d>();
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
		if (nset)
			writeNodalValues<double>(*nset, a, map);
		else
		{
			writeNodalValues<double>(mesh, a, [&](const FENode& node) {
				FEMaterialPoint mp;
				mp.m_r0 = node.m_r0;
				mp.m_index = -1;
				double v = map(mp);
				return v;
			});
		}

		return true;
	}
	else if (m_param.type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = m_param.value<FEParamVec3>();
		FENodeSet* nset = dynamic_cast<FENodeSet*>(map.GetItemList());
		if (nset == 0) return false;

		// write the nodal values
		writeNodalValues<vec3d>(*nset, a, map);

		return true;
	}


	return false;
}

//-----------------------------------------------------------------------------
FEPlotPIDController::FEPlotPIDController(FEModel* pfem) : FEPlotGlobalData(pfem, PLT_ARRAY)
{
	m_pid = nullptr;
	SetArraySize(3);
	std::vector<std::string> names;
	names.push_back("measurement");
	names.push_back("error");
	names.push_back("output");
	SetArrayNames(names);
}

bool FEPlotPIDController::SetFilter(const char* sz)
{
	if (sz == nullptr) return false;
	FEModel* fem = GetFEModel(); assert(fem);
	if (fem == nullptr) return false;

	for (int i = 0; i < fem->LoadControllers(); ++i)
	{
		FEPIDController* pid = dynamic_cast<FEPIDController*>(fem->GetLoadController(i));
		if (pid && (pid->GetName() == string(sz)))
		{
			m_pid = pid;
			return true;
		}
	}

	m_pid = nullptr;
	return false;
}

bool FEPlotPIDController::Save(FEDataStream& a)
{
	if (m_pid == nullptr) return false;
	a << m_pid->GetParameterValue();
	a << m_pid->GetError();
	a << m_pid->Value();
	return true;
}

//=============================================================================

//-----------------------------------------------------------------------------
FEPlotMeshData::FEPlotMeshData(FEModel* pfem) : FEPlotData(pfem)
{
	m_map = nullptr;
	m_dom = nullptr;
}

//-----------------------------------------------------------------------------
// This plot field requires a filter which defines the material name and 
// the material parameter in the format [materialname.parametername].
bool FEPlotMeshData::SetFilter(const char* sz)
{
	// store the filter for serialization
	m_filter = sz;

	// find the map
	m_map = GetFEModel()->GetMesh().FindDataMap(sz);
	if (m_map == nullptr) return false;

	switch (m_map->DataMapType())
	{
	case FE_DOMAIN_MAP:
	{
		FEDomainMap* map = dynamic_cast<FEDomainMap*>(m_map); assert(map);
		FEDataType dataType = map->DataType();

		SetRegionType(FE_REGION_DOMAIN);
		SetStorageFormat((Storage_Fmt)map->StorageFormat());

		if (map->StorageFormat() == FMT_MATPOINTS) SetStorageFormat(FMT_MULT);
		
		switch (dataType)
		{
		case FE_DOUBLE: SetVarType(PLT_FLOAT); break;
		default:
			assert(false);
			return false;
		}

		const FEElementSet* set = map->GetElementSet();
		if (set == nullptr) return false;

		m_dom = &(set->GetDomainList());
	}
	break;
	default:
		assert(false);
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPlotMeshData::Serialize(DumpStream& ar)
{
	FEPlotData::Serialize(ar);
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
		ar << m_filter;
	else
	{
		string filter;
		ar >> filter;
		SetFilter(filter.c_str());
	}
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotMeshData::Save(FEDomain& dom, FEDataStream& a)
{
	FEDomainMap* map = dynamic_cast<FEDomainMap*>(m_map);
	if (map == nullptr) return false;
	if (m_dom == nullptr) return false;

	if (m_dom->IsMember(&dom))
	{
		if (StorageFormat() == FMT_NODE)
		{
			if (DataType() == PLT_FLOAT)
			{
				int n = dom.Nodes();
				for (int i = 0; i < n; ++i)
				{
					int m = dom.NodeIndex(i);
					double v = map->NodalValue(m);
					a << v;
				}
				return true;
			}
		}
		else if (StorageFormat() == FMT_ITEM)
		{
			if (DataType() == PLT_FLOAT)
			{
				int n = dom.Elements();
				for (int i = 0; i < n; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					FEMaterialPoint mp;
					mp.m_elem = &el;
					mp.m_index = 0;
					double v = map->value(mp);
					a << v;
				}
				return true;
			}
		}
		else if (StorageFormat() == FMT_MULT)
		{
			if (DataType() == PLT_FLOAT)
			{
				int n = dom.Elements();
				for (int i = 0; i < n; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					for (int j = 0; j < el.Nodes(); ++j)
					{
						double v = map->value<double>(i, j);
						a << v;
					}
				}
				return true;
			}
			else if (DataType() == PLT_VEC3F)
			{
				int n = dom.Elements();
				for (int i = 0; i < n; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					for (int j = 0; j < el.Nodes(); ++j)
					{
						vec3d v = map->value<vec3d>(i, j);
						a << v;
					}
				}
				return true;
			}
		}
		else if (map->StorageFormat() == FMT_MATPOINTS)
		{
			// Note that for this map type, the plot format was changed to FMT_MULT
			if (DataType() == PLT_FLOAT)
			{
				int NE = dom.Elements();
				double vi[FEElement::MAX_INTPOINTS] = { 0 };
				double vn[FEElement::MAX_NODES] = { 0 };
				for (int i = 0; i < NE; ++i)
				{
					FEElement& el = dom.ElementRef(i);
					for (int j = 0; j < el.GaussPoints(); ++j)
					{
						FEMaterialPoint mp;
						mp.m_elem = &el;
						mp.m_index = j;
						vi[j] = map->value(mp);
					}

					el.project_to_nodes(vi, vn);
					for (int j=0; j<el.Nodes(); ++j)
						a << vn[j];
				}
				return true;
			}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
// The Save function stores the material parameter data to the plot file.
bool FEPlotMeshData::Save(FESurface& dom, FEDataStream& a)
{
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotMeshData::Save(FEMesh& mesh, FEDataStream& a)
{
	return false;
}

FEPlotFieldVariable::FEPlotFieldVariable(FEModel* pfem) : FEPlotNodeData(pfem, Var_Type::PLT_FLOAT, Storage_Fmt::FMT_ITEM)
{

}

bool FEPlotFieldVariable::SetFilter(const char* sz)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int nvar = dofs.GetVariableIndex(sz);
	if (nvar == -1) return false;
	dofs.GetDOFList(sz, m_dofs);
	int vartype = dofs.GetVariableType(nvar);
	switch (vartype)
	{
	case VAR_SCALAR: SetVarType(PLT_FLOAT); break;
	case VAR_VEC3  : SetVarType(PLT_VEC3F); break;
	case VAR_ARRAY :
	{
		SetVarType(Var_Type::PLT_ARRAY);
		SetArraySize(m_dofs.size());
		vector<string> dofNames;
		for (int i = 0; i < m_dofs.size(); ++i) dofNames.push_back(dofs.GetDOFName(nvar, i));
		SetArrayNames(dofNames);
	}
	break;
	default:
		return false;
	}
	return true;
}

bool FEPlotFieldVariable::Save(FEMesh& mesh, FEDataStream& a)
{
	if (m_dofs.empty()) return false;
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < m_dofs.size(); ++j)
		{
			a << node.get(m_dofs[j]);
		}
	}
	return true;
}

void FEPlotFieldVariable::Serialize(DumpStream& ar)
{
	FEPlotNodeData::Serialize(ar);
	ar & m_dofs;
}
