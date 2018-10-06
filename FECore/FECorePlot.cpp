#include "stdafx.h"
#include "FECorePlot.h"
#include "FEMaterial.h"
#include "FESolidDomain.h"
#include "FEModelParam.h"
#include "FEBodyLoad.h"
#include "FEPlotData.h"

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
				if (bl) domList = &bl->GetDomaintList();
				else return false;
			}
		}
		else domList = &elset->GetDomainList();
		if (domList->IsMember(&dom) == false) return false;

		FESolidDomain& sd = dynamic_cast<FESolidDomain&>(dom);

		if (param->type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& mapDouble = dynamic_cast<FEParamDouble&>(map);
			writeNodalProjectedElementValues(sd, mapDouble, a);
		}
		else if (param->type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& mapVec3 = dynamic_cast<FEParamVec3&>(map);
			writeNodalProjectedElementValues(sd, mapVec3, a);
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
				gi[j] = map(mp);
			}

			e.FEElement::project_to_nodes(gi, gn);

			// store the result
			// NOTE: Note that we always need to store 10 entries. This is because of a limitation of the plot file format.
			for (int j = 0; j < 10; ++j) a << gn[j];
		}
	}
	else if (param->type() == FE_PARAM_VEC3D_MAPPED)
	{
		FEParamVec3& map = param->value<FEParamVec3>();

		FEFacetSet* surf = dynamic_cast<FEFacetSet*>(map.GetItemList());
		if (surf != dom.GetFacetSet()) return false;

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
				gi[j] = map(mp);
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

bool FEPlotParameter::Save(FEMesh& mesh, FEDataStream& a)
{
	if (m_param == 0) return false;
	FEParam* param = m_param;

	if (param->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		FEParamDouble& map = param->value<FEParamDouble>();

		int NN = mesh.Nodes();
		vector<double> data(NN, 0.0);

		FENodeSet* nset = dynamic_cast<FENodeSet*>(map.GetItemList());
		if (nset == 0) return false;

		const std::vector<int> nodeList = nset->GetNodeList();
		FEMaterialPoint mp;
		for (int i = 0; i < nset->size(); ++i)
		{
			int nodeId = nodeList[i];
			FENode& node = mesh.Node(nodeId);
			mp.m_r0 = node.m_r0;
			mp.m_index = i;

			double vi = map(mp);

			data[nodeId] = vi;
		}
		a << data;

		return true;
	}

	return false;
}
