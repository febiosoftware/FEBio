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
#include "FEBioPlotFile.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEDataExport.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include <FECore/FESurface.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/log.h>

FEBioPlotFile::DICTIONARY_ITEM::DICTIONARY_ITEM()
{
	m_psave = 0;
	m_ntype = 0;
	m_nfmt = 0;
	m_arraySize = 0;
	m_szname[0] = 0;
	m_szunit[0] = 0;
}

FEBioPlotFile::DICTIONARY_ITEM::DICTIONARY_ITEM(const FEBioPlotFile::DICTIONARY_ITEM& item)
{
	m_psave = item.m_psave;
	m_ntype = item.m_ntype;
	m_nfmt = item.m_nfmt;
	m_arraySize = item.m_arraySize;
	m_arrayNames = item.m_arrayNames;
	m_szname[0] = 0;
	m_szunit[0] = 0;
	if (item.m_szname[0]) strcpy(m_szname, item.m_szname);
	if (item.m_szunit[0]) strcpy(m_szunit, item.m_szunit);
}

class FEPlotSurfaceDataExport : public FEPlotData
{
public:
	FEPlotSurfaceDataExport(FEModel* fem, const char* szname, Var_Type itype, Storage_Fmt fmt) : FEPlotData(fem, FE_REGION_SURFACE, itype, fmt) { m_szname = szname; }
	void Save(FEModel& fem, PltArchive& ar)
	{
		FEMesh& mesh = fem.GetMesh();
		int NS = mesh.Surfaces();
		for (int i = 0; i<NS; ++i)
		{
			FESurface& s = mesh.Surface(i);
			int ND = s.DataExports();
			if (ND > 0)
			{
				for (int j = 0; j<ND; ++j)
				{
					FEDataExport* pd = s.GetDataExport(j);
					if (strcmp(pd->m_szname, m_szname) == 0)
					{
						FEDataStream d;
						pd->Serialize(d);
						ar.WriteData(i + 1, d.data());
						break;
					}
				}
			}
		}
	}

private:
	const char*		m_szname;
};

class FEPlotDomainDataExport : public FEPlotData
{
public:
	FEPlotDomainDataExport(FEModel* fem, const char* szname, Var_Type itype, Storage_Fmt fmt) : FEPlotData(fem, FE_REGION_DOMAIN, itype, fmt) { m_szname = szname; }
	void Save(FEModel& fem, PltArchive& ar)
	{
		FEMesh& mesh = fem.GetMesh();
		int NDOMS = mesh.Domains();
		for (int i = 0; i<NDOMS; ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			int NDATA = dom.DataExports();
			if (NDATA > 0)
			{
				for (int j = 0; j<NDATA; ++j)
				{
					FEDataExport* pd = dom.GetDataExport(j);
					if (strcmp(pd->m_szname, m_szname) == 0)
					{
						FEDataStream d;
						pd->Serialize(d);
						ar.WriteData(i + 1, d.data());
						break;
					}
				}
			}
		}
	}

private:
	const char*		m_szname;
};


class FEBioPlotVariable : public FEPlotNodeData
{
public:
	FEBioPlotVariable(FEModel* fem, const char* szname, Var_Type itype, Storage_Fmt fmt) : FEPlotNodeData(fem, itype, fmt) { strcpy(m_szname, szname); }
	bool Save(FEMesh& mesh, FEDataStream& str)
	{
		// get the DOFS
		FEModel& fem = *GetFEModel();
		DOFS& dofs = fem.GetDOFS();

		// see if this variable exists
		int nvar = dofs.GetVariableIndex(m_szname);
		if (nvar < 0) return false;

		// get the size of the variable
		int n = dofs.GetVariableSize(nvar);
		if (n == 0) return false;

		// get the start index of the DOFS
		int ndof = dofs.GetDOF(nvar, 0);
		if (ndof < 0) return false;

		// store the nodal data
		int NN = mesh.Nodes();
		for (int i = 0; i<NN; ++i)
		{
			FENode& node = mesh.Node(i);
			for (int j = 0; j<n; ++j) str << node.get(ndof + j);
		}

		return true;
	}

private:
	char	m_szname[256];
};

class FEPlotArrayVariable : public FEPlotDomainData
{
public:
	FEPlotArrayVariable(FEModel* fem, const char* szname, int index) : FEPlotDomainData(fem, PLT_FLOAT, FMT_NODE) 
	{ 
		strcpy(m_szname, szname); 
		m_index = index;
		if (index == -1)
		{
			// get the DOFS
			FEModel& fem = *GetFEModel();
			DOFS& dofs = fem.GetDOFS();

			// see if this variable exists
			int nvar = dofs.GetVariableIndex(m_szname);
			if (nvar >= 0)
			{
				// get the size of the variable
				int n = dofs.GetVariableSize(nvar);
				if (n > 0)
				{
					SetVarType(PLT_ARRAY);
					SetArraySize(n);
				}
			}
		}
	}
	bool Save(FEDomain& D, FEDataStream& a)
	{
		// get the DOFS
		FEModel& fem = *GetFEModel();
		DOFS& dofs = fem.GetDOFS();

		// see if this variable exists
		int nvar = dofs.GetVariableIndex(m_szname);
		if (nvar < 0) return false;

		// get the size of the variable
		int n = dofs.GetVariableSize(nvar);
		if (n == 0) return false;

		// get the start index of the DOFS
		if (m_index >= 0)
		{
			int ndof = dofs.GetDOF(nvar, m_index);
			if (ndof < 0) return false;

			// see if this domain contains this dof
			const FEDofList& domDofs = D.GetDOFList();
			bool bfound = false;
			for (int i = 0; i < (int)domDofs.Size(); ++i)
			{
				if (domDofs[i] == ndof)
				{
					bfound = true;
					break;
				}
			}
			if (bfound == false) return false;

			// store the nodal data
			int NN = D.Nodes();
			for (int i = 0; i < NN; ++i)
			{
				FENode& node = D.Node(i);
				a << node.get(ndof);
			}
		}
		else
		{
			const FEDofList& domDofs = D.GetDOFList();

			vector<int> dofList(n, -1);
			for (int dof_i = 0; dof_i < n; dof_i++)
			{
				int ndof = dofs.GetDOF(nvar, dof_i);
				// see if this domain contains this dof
				bool bfound = false;
				for (int i = 0; i < (int)domDofs.Size(); ++i)
				{
					if (domDofs[i] == ndof)
					{
						bfound = true;
						dofList[dof_i] = i;
						break;
					}
				}
			}

			// store the nodal data
			int NN = D.Nodes();
			for (int i = 0; i < NN; ++i)
			{
				FENode& node = D.Node(i);
				for (int dof_i = 0; dof_i < n; dof_i++)
				{
					int ndof = dofList[dof_i];
					if (ndof < 0)
					{
						a << 0.0;
					}
					else
					{
						// store the nodal data
						a << node.get(ndof);
					}
				}
			}
		}

		return true;
	}

private:
	char	m_szname[256];
	int		m_index;
};

//-----------------------------------------------------------------------------
//! Adds a variable to the plot file. 
//! 
//! The name of the filter can be composed of three parts and in general takes on
//! the following format.
//!
//! szname = "field_name[filter]=alias". 
//!
//! field_name = This is the actual filter naeme as it is registered with the framework.
//! filter     = This is a filter that is used to resolve ambiguities.
//! alias      = This is an alternative name for the field variable.
//!
//! The alias is optional but can be used by post-processing software to present an alternative
//! (often simpler) name for the field variable than the default field_name + filter combo.
//! 
//! Whether a filter is required depends entirely on the field variable. Most field variables don't
//! require it, but some do in order to resolve an ambiguity. For instance, the "parameter" field
//! allows users to plot the spatially varying value of a material parameter. The filter is used to
//! specify the material and parameter name. 
//!
//! The filter can be a numerical value or a string. If it's a string then it must be enclosed in
//! single quotes. 
//!
//! szname = "field_name[12]"  \\ example of a numerical filter
//! szname = "field_name['val'] \\ example of a string filter
//!
//! The interpretation of these filters is entirely left up to the field variable. 
//!
bool FEBioPlotFile::Dictionary::AddVariable(FEModel* pfem, const char* szname, vector<int>& item, const char* szdom)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// create a copy so we can strip the alias and the filter from the name
	char sz[1024] = { 0 };
	strcpy(sz, szname);

	// This is the name that will be stored in the plot file
	const char* szfield = sz;

	// see if there is an alias defined
	char* ch = strchr(sz, '=');
	if (ch)
	{
		// replace the equal sign with a null character.
		*ch++ = 0;

		// make sure there is an alias
		if (ch == 0) return false;

		// store the alias instead of the actual name
		szfield = ch;
	}

	// extract the filter
	char* szflt = strchr(sz, '[');
	int index = 0;
	int fltType = 0;
	if (szflt)
	{
		*szflt++ = 0;
		char* ch = strrchr(szflt, ']');
		if (ch == 0) return false;
		*ch = 0;

		// see if the filter is a number or a string
		ch = strchr(szflt, '\'');
		if (ch)
		{
			*szflt++ = 0;
			// find the end character
			char* ch2 = strrchr(szflt, '\'');
			if (ch2 == 0) return false;
			*ch2 = 0;
		}
		else
		{
			fltType = 1;
			index = atoi(szflt);
		}
	}

	// create the plot variable
	FEPlotData* ps = fecore_new<FEPlotData>(sz, pfem);
	if (ps)
	{
		// set the optional item list and filter
		ps->SetItemList(item);
		if (szflt)
		{
			if (fltType == 0)
			{
				if (ps->SetFilter(szflt) == false) return false;
			}
			else if (fltType == 1)
			{
				if (ps->SetFilter(index) == false) return false;
			}
		}

		// add the field to the plot file
		ps->SetDomainName(szdom);
		switch (ps->RegionType())
		{
		case FE_REGION_NODE: return AddNodalVariable(ps, szfield, item);
		case FE_REGION_DOMAIN: return AddDomainVariable(ps, szfield, item);
		case FE_REGION_SURFACE: return AddSurfaceVariable(ps, szfield, item);
		default:
			assert(false);
			return false;
		}
	}
	else
	{
		// If we get here then this variable is not a plot field.
		// But let's see if it is an export variable from a domain
		// Check the surfaces first
		FEMesh& mesh = pfem->GetMesh();
		for (int i = 0; i<mesh.Surfaces(); ++i)
		{
			FESurface& s = mesh.Surface(i);
			int ND = s.DataExports();
			for (int j = 0; j<ND; ++j)
			{
				FEDataExport* pd = s.GetDataExport(j);
				if (strcmp(pd->m_szname, szname) == 0)
				{
					// We have a match. Create a plot field for this export
					ps = new FEPlotSurfaceDataExport(pfem, pd->m_szname, pd->m_type, pd->m_fmt);
					return AddSurfaceVariable(ps, szname, item);
				}
			}
		}

		// now the domains.
		for (int i = 0; i<mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			int ND = dom.DataExports();
			for (int j = 0; j<ND; ++j)
			{
				FEDataExport* pd = dom.GetDataExport(j);
				if (strcmp(pd->m_szname, szname) == 0)
				{
					// We have a match. Create a plot field for this export
					ps = new FEPlotDomainDataExport(pfem, pd->m_szname, pd->m_type, pd->m_fmt);
					return AddDomainVariable(ps, szname, item);
				}
			}
		}

		// If we still didn't find it, maybe it's a model variable.
		DOFS& dofs = pfem->GetDOFS();
		int nvar = dofs.GetVariableIndex(sz);
		if (nvar >= 0)
		{
			int vartype = dofs.GetVariableType(nvar);
			if (vartype == VAR_SCALAR)
			{
				ps = new FEBioPlotVariable(pfem, sz, PLT_FLOAT, FMT_NODE);
				return AddNodalVariable(ps, szname, item);
			}
			else if (vartype == VAR_VEC3)
			{
				ps = new FEBioPlotVariable(pfem, sz, PLT_VEC3F, FMT_NODE);
				return AddNodalVariable(ps, szname, item);
			}
			else if (vartype == VAR_ARRAY)
			{
				int ndofs = dofs.GetVariableSize(sz);
				if (fltType == 0)
				{
					index = -1;
					if (szflt)
					{

						index = dofs.GetIndex(sz, szflt);
						if ((index < 0) || (index >= ndofs)) return false;
					}
				}

				ps = new FEPlotArrayVariable(pfem, sz, index);
				return AddDomainVariable(ps, szname, item);
			}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Dictionary::AddGlobalVariable(FEPlotData* ps, const char* szname)
{
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Dictionary::AddMaterialVariable(FEPlotData* ps, const char* szname)
{
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Dictionary::AddNodalVariable(FEPlotData* ps, const char* szname, vector<int>& item)
{
	assert(ps->RegionType()==FE_REGION_NODE);
	if (ps->RegionType()==FE_REGION_NODE)
	{
		DICTIONARY_ITEM it;
		it.m_ntype = ps->DataType();
		it.m_nfmt  = ps->StorageFormat();
		it.m_psave = ps;
		it.m_arraySize = ps->GetArraysize();
		it.m_arrayNames = ps->GetArrayNames();
		strcpy(it.m_szname, szname);
		if (ps->GetUnits())
		{
			strcpy(it.m_szunit, ps->GetUnits());
		}
		m_Node.push_back(it);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Dictionary::AddDomainVariable(FEPlotData* ps, const char* szname, vector<int>& item)
{
	assert(ps->RegionType()==FE_REGION_DOMAIN);
	if (ps->RegionType()==FE_REGION_DOMAIN)
	{
		DICTIONARY_ITEM it;
		it.m_ntype = ps->DataType();
		it.m_nfmt  = ps->StorageFormat();
		it.m_psave = ps;
		it.m_arraySize = ps->GetArraysize();
		it.m_arrayNames = ps->GetArrayNames();
		strcpy(it.m_szname, szname);
		if (ps->GetUnits())
		{
			strcpy(it.m_szunit, ps->GetUnits());
		}
		m_Elem.push_back(it);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Dictionary::AddSurfaceVariable(FEPlotData* ps, const char* szname, vector<int>& item)
{
	assert(ps->RegionType()==FE_REGION_SURFACE);
	if (ps->RegionType()==FE_REGION_SURFACE)
	{
		DICTIONARY_ITEM it;
		it.m_ntype = ps->DataType();
		it.m_nfmt  = ps->StorageFormat();
		it.m_psave = ps;
		it.m_arraySize = ps->GetArraysize();
		it.m_arrayNames = ps->GetArrayNames();
		strcpy(it.m_szname, szname);
		if (ps->GetUnits())
		{
			strcpy(it.m_szunit, ps->GetUnits());
		}
		m_Face.push_back(it);
		return true;
	}
	return false;
}


//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::Defaults(FEModel& fem)
{
	// First we build the dictionary
	// get the mesh
	FEMesh& m = fem.GetMesh();

	// Define default variables
	if (m_Node.empty() && m_Elem.empty() && m_Face.empty())
	{
		vector<int> l; // empty list
		AddVariable(&fem, "displacement", l);
		AddVariable(&fem, "stress", l);
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::Clear()
{
	m_Glob.clear();
	m_Mat.clear();
	m_Node.clear();
	m_Elem.clear();
	m_Face.clear();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::PlotObject::AddData(const char* szname, Var_Type type, FEPlotData* psave)
{
	DICTIONARY_ITEM item;
	item.m_psave = nullptr;

	item.m_szname[0] = 0;
	int l = strlen(szname);
	if (l >= STR_SIZE) l = STR_SIZE - 1;
	strncpy(item.m_szname, szname, l);
	item.m_szname[l] = 0;

	item.m_ntype = type;
	item.m_nfmt = FMT_ITEM;
	item.m_psave = psave;

	m_data.push_back(item);
}


//=============================================================================

FEBioPlotFile::FEBioPlotFile(FEModel* fem) : PlotFile(fem)
{
	m_ncompress = 0;
	m_meshesWritten = 0;
	m_exportUnitsFlag = false;
}

//-----------------------------------------------------------------------------
FEBioPlotFile::~FEBioPlotFile(void)
{
	// close the archive
	Close();

	// clear all arrays
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Glob.begin();
	for (int i=0; i<(int) m_dic.m_Glob.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Mat.begin();
	for (int i=0; i<(int) m_dic.m_Mat.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Node.begin();
	for (int i=0; i<(int) m_dic.m_Node.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Elem.begin();
	for (int i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Face.begin();
	for (int i=0; i<(int) m_dic.m_Face.size(); ++i, ++it) delete it->m_psave;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::AddVariable(FEPlotData* ps, const char* szname)
{
	vector<int> dummy;
	switch (ps->RegionType())
	{
	case FE_REGION_NODE: return m_dic.AddNodalVariable(ps, szname, dummy);
	case FE_REGION_DOMAIN: return m_dic.AddDomainVariable(ps, szname, dummy);
	case FE_REGION_SURFACE: return m_dic.AddSurfaceVariable(ps, szname, dummy);
	default:
		assert(false);
		return false;
	}
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::AddVariable(const char* sz)
{
	vector<int> dummy;
	return AddVariable(sz, dummy);
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::AddVariable(const char* sz, vector<int>& item, const char* szdom)
{ 
	return m_dic.AddVariable(GetFEModel(), sz, item, szdom); 
}

//-----------------------------------------------------------------------------
int FEBioPlotFile::PointObjects()
{
	return (int) m_Points.size();
}

//-----------------------------------------------------------------------------
FEBioPlotFile::PointObject* FEBioPlotFile::GetPointObject(int i)
{
	if ((i >= 0) && (i < PointObjects())) return m_Points[i];
	return nullptr;
}

//-----------------------------------------------------------------------------
FEBioPlotFile::PointObject* FEBioPlotFile::AddPointObject(const std::string& name)
{
	PointObject* po = new PointObject;
	m_Points.push_back(po);
	po->m_name = name;
	po->m_id = m_Points.size();
	return po;
}

//-----------------------------------------------------------------------------
int FEBioPlotFile::LineObjects()
{
	return (int)m_Lines.size();
}

//-----------------------------------------------------------------------------
FEBioPlotFile::LineObject* FEBioPlotFile::GetLineObject(int i)
{
	if ((i >= 0) && (i < LineObjects())) return m_Lines[i];
	return nullptr;
}

//-----------------------------------------------------------------------------
FEBioPlotFile::LineObject* FEBioPlotFile::AddLineObject(const std::string& name)
{
	LineObject* po = new LineObject;
	m_Lines.push_back(po);
	po->m_name = name;
	po->m_id = m_Lines.size();
	return po;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::SetCompression(int n)
{
	m_ncompress = n;
}

//-----------------------------------------------------------------------------
//! set the version string
void FEBioPlotFile::SetSoftwareString(const std::string& softwareString)
{
	m_softwareString = softwareString;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::IsValid() const
{
	return m_ar.IsValid();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Close()
{
	m_ar.Close();
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Open(const char *szfile)
{
	FEModel* fem = GetFEModel();

	// open the archive
	m_ar.Create(szfile);

	// set compression
	FEPlotDataStore& pltData = fem->GetPlotDataStore();
	SetCompression(pltData.GetPlotCompression());

	// add plot variables
	for (int n = 0; n < pltData.PlotVariables(); ++n)
	{
		FEPlotVariable& vi = pltData.GetPlotVariable(n);
		const std::string& varName = vi.Name();
		const std::string& domName = vi.DomainName();

		// add the plot output variable
		if (AddVariable(varName.c_str(), vi.m_item, domName.c_str()) == false)
		{
			feLog("FATAL ERROR: Output variable \"%s\" is not defined\n", varName.c_str());
			throw "FATAL ERROR";
		}
	}

	try
	{
		// write the root element
		if (WriteRoot(*fem) == false) return false;

		// write the mesh section
		if (WriteMeshSection(*fem) == false) return false;
	}
	catch (...)
	{
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteRoot(FEModel& fem)
{
	// write the root element
	// (don't compress this section)
	m_ar.SetCompression(0);
	m_ar.BeginChunk(PLT_ROOT);
	{
		// --- save the header file ---
		m_ar.BeginChunk(PLT_HEADER);
		{
			if (WriteHeader(fem) == false) return false;
		}
		m_ar.EndChunk();

		// --- save the dictionary ---
		m_ar.BeginChunk(PLT_DICTIONARY);
		{
			if (WriteDictionary(fem) == false) return false;
		}
		m_ar.EndChunk();
	}
	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteHeader(FEModel& fem)
{
	// setup the header
	unsigned int nversion = PLT_VERSION;

	// output header
	m_ar.WriteChunk(PLT_HDR_VERSION, nversion);

	// compression flag
	m_ar.WriteChunk(PLT_HDR_COMPRESSION, m_ncompress);

	// software flag
	if (m_softwareString.empty() == false)
	{
		const char* sz = m_softwareString.c_str();
		m_ar.WriteChunk(PLT_HDR_SOFTWARE, sz);
	}

	// units flag
	m_exportUnitsFlag = false;
	const char* szunits = fem.GetUnits();
	if (szunits != nullptr)
	{
		m_exportUnitsFlag = true;
		m_ar.WriteChunk(PLT_HDR_UNITS, szunits);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteDictionary(FEModel& fem)
{
	// setup defaults for the dictionary
	m_dic.Defaults(fem);

	// Next, we save the dictionary
	// Global variables
	if (!m_dic.m_Glob.empty())
	{
		m_ar.BeginChunk(PLT_DIC_GLOBAL);
		{
			WriteDicList(m_dic.m_Glob);
		}
		m_ar.EndChunk();
	}

	// store nodal variables
	if (!m_dic.m_Node.empty())
	{
		m_ar.BeginChunk(PLT_DIC_NODAL);
		{
			WriteDicList(m_dic.m_Node);
		}
		m_ar.EndChunk();
	}

	// store element variables
	if (!m_dic.m_Elem.empty())
	{
		m_ar.BeginChunk(PLT_DIC_DOMAIN);
		{
			WriteDicList(m_dic.m_Elem);
		}
		m_ar.EndChunk();
	}

	// store surface data
	if (!m_dic.m_Face.empty())
	{
		m_ar.BeginChunk(PLT_DIC_SURFACE);
		{
			WriteDicList(m_dic.m_Face);
		}
		m_ar.EndChunk();
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDicList(list<FEBioPlotFile::DICTIONARY_ITEM>& dic)
{
	int N = (int) dic.size();
	list<DICTIONARY_ITEM>::iterator pi = dic.begin();
	for (int i=0; i<N; ++i, ++pi)
	{
		m_ar.BeginChunk(PLT_DIC_ITEM);
		{
			WriteDictionaryItem(*pi);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDictionaryItem(DICTIONARY_ITEM& it)
{
	m_ar.WriteChunk(PLT_DIC_ITEM_TYPE, it.m_ntype);
	m_ar.WriteChunk(PLT_DIC_ITEM_FMT, it.m_nfmt);
	m_ar.WriteChunk(PLT_DIC_ITEM_ARRAYSIZE, it.m_arraySize);
	if ((it.m_arraySize > 0) && (it.m_arrayNames.size() == it.m_arraySize))
	{
		for (int i = 0; i < (int)it.m_arraySize; ++i)
		{
			string& si = it.m_arrayNames[i];
			const char* c = si.c_str();
			m_ar.WriteChunk(PLT_DIC_ITEM_ARRAYNAME, (char*)c, STR_SIZE);
		}
	}
	m_ar.WriteChunk(PLT_DIC_ITEM_NAME, it.m_szname, STR_SIZE);

	if (m_exportUnitsFlag && it.m_szunit && it.m_szunit[0])
	{
		m_ar.WriteChunk(PLT_DIC_ITEM_UNITS, it.m_szunit, STR_SIZE);
	}
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteMeshSection(FEModel& fem)
{
	// get the mesh
	FEMesh& m = fem.GetMesh();

	m_ar.BeginChunk(PLT_MESH);
	{
		// node section
		m_ar.BeginChunk(PLT_NODE_SECTION);
		{
			WriteNodeSection(m);
		}
		m_ar.EndChunk();

		// domain section
		m_ar.BeginChunk(PLT_DOMAIN_SECTION);
		{
			WriteDomainSection(m);
		}
		m_ar.EndChunk();

		// surface section
		if (m.Surfaces() > 0)
		{
			BuildSurfaceTable();
			m_ar.BeginChunk(PLT_SURFACE_SECTION);
			{
				WriteSurfaceSection(m);
			}
			m_ar.EndChunk();
		}

		// node sets
		if (m.NodeSets() > 0)
		{
			m_ar.BeginChunk(PLT_NODESET_SECTION);
			{
				WriteNodeSetSection(m);
			}
			m_ar.EndChunk();
		}

		// element sets
		if (m.ElementSets() > 0)
		{
			m_ar.BeginChunk(PLT_ELEMENTSET_SECTION);
			{
				WriteElementSetSection(m);
			}
			m_ar.EndChunk();
		}

		// facet sets
		if (m.FacetSets() > 0)
		{
			m_ar.BeginChunk(PLT_FACETSET_SECTION);
			{
				WriteFacetSetSection(m);
			}
			m_ar.EndChunk();
		}
		// parts
		// (we write the materials as parts)
		if (fem.Materials() > 0)
		{
			m_ar.BeginChunk(PLT_PARTS_SECTION);
			{
				WritePartsSection(fem);
			}
			m_ar.EndChunk();
		}

		// additional objects
		if (m_meshesWritten == 0)
		{
			if ((m_Points.size() > 0) || (m_Lines.size() > 0))
			{
				m_ar.BeginChunk(PLT_OBJECTS_SECTION);
				{
					WriteObjectsSection();
				}
				m_ar.EndChunk();
			}
		}
	}
	m_ar.EndChunk();

	m_meshesWritten++;

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeSection(FEMesh& m)
{
	// write the node header
	m_ar.BeginChunk(PLT_NODE_HEADER);
	{
		int NN = m.Nodes();
		int dim = 3;
		m_ar.WriteChunk(PLT_NODE_SIZE, NN);
		m_ar.WriteChunk(PLT_NODE_DIM , dim);
//		m_ar.WriteChunk(PLT_NODE_NAME, "AllNodes");
	}
	m_ar.EndChunk();

	// write the reference coordinates
	int NN = m.Nodes();
	vector<float> X(4*NN);
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		*((int*) (&X[0] + 4*i)) = i;
		X[4*i+1] = (float) node.m_r0.x;
		X[4*i+2] = (float) node.m_r0.y;
		X[4*i+3] = (float) node.m_r0.z;
	}
	m_ar.WriteChunk(PLT_NODE_COORDS, X);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDomainSection(FEMesh& m)
{
	// write all domains
	for (int nd = 0; nd<m.Domains(); ++nd)
	{
		FEDomain& dom = m.Domain(nd);
		m_ar.BeginChunk(PLT_DOMAIN);
		{
			switch (dom.Class())
			{
			case FE_DOMAIN_SOLID   : WriteSolidDomain   (static_cast<FESolidDomain&   >(dom)); break;
			case FE_DOMAIN_SHELL   : WriteShellDomain   (static_cast<FEShellDomain&   >(dom)); break;
			case FE_DOMAIN_BEAM    : WriteBeamDomain    (static_cast<FEBeamDomain&    >(dom)); break;
			case FE_DOMAIN_DISCRETE: WriteDiscreteDomain(static_cast<FEDiscreteDomain&>(dom)); break;
            case FE_DOMAIN_2D      : WriteDomain2D      (static_cast<FEDomain2D&      >(dom)); break;
			}
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSolidDomain(FESolidDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	assert(mid > 0);
	int eshape = dom.GetElementShape();

	int i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 0;
	int dtype = 0;
	switch (eshape)
	{
		case ET_HEX8   : ne =  8; dtype = PLT_ELEM_HEX; break;
		case ET_PENTA6 : ne =  6; dtype = PLT_ELEM_PENTA; break;
		case ET_TET4   : ne =  4; dtype = PLT_ELEM_TET4; break;
		case ET_TET5   : ne =  5; dtype = PLT_ELEM_TET5; break;
		case ET_TET10  : ne = 10; dtype = PLT_ELEM_TET10; break;
		case ET_TET15  : ne = 15; dtype = PLT_ELEM_TET15; break;
		case ET_HEX20  : ne = 20; dtype = PLT_ELEM_HEX20; break;
		case ET_HEX27  : ne = 27; dtype = PLT_ELEM_HEX27; break;
		case ET_TET20  : ne = 20; dtype = PLT_ELEM_TET20; break;
        case ET_PENTA15: ne = 15; dtype = PLT_ELEM_PENTA15; break;
		case ET_PYRA5  : ne =  5; dtype = PLT_ELEM_PYRA5; break;
        case ET_PYRA13 : ne = 13; dtype = PLT_ELEM_PYRA13; break;
        default:
			assert(false);
	}

	// write the header
	m_ar.BeginChunk(PLT_DOMAIN_HDR);
	{
		m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
		m_ar.WriteChunk(PLT_DOM_PART_ID  ,   mid);
		m_ar.WriteChunk(PLT_DOM_ELEMS    ,    NE);
		m_ar.WriteChunk(PLT_DOM_NAME     , dom.GetName());
	}
	m_ar.EndChunk();

	// write the element list
	int n[FEElement::MAX_NODES + 1];
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST);
	{
		for (i=0; i<NE; ++i)
		{
			FESolidElement& el = dom.Element(i);
			n[0] = el.GetID();
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j];
			m_ar.WriteChunk(PLT_ELEMENT, n, ne+1);
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteShellDomain(FEShellDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	assert(mid > 0);
	int etype = dom.GetElementType();

	int i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 0;
	int dtype = 0;
	switch (etype)
	{
        case FE_SHELL_QUAD4G8 :
        case FE_SHELL_QUAD4G12 :
            ne = 4; dtype = PLT_ELEM_QUAD; break;
        case FE_SHELL_TRI3G6  :
        case FE_SHELL_TRI3G9  :
            ne = 3; dtype = PLT_ELEM_TRI; break;
        case FE_SHELL_QUAD8G18:
        case FE_SHELL_QUAD8G27:
            ne = 8; dtype = PLT_ELEM_QUAD8; break;
        case FE_SHELL_TRI6G14 :
        case FE_SHELL_TRI6G21 :
            ne = 6; dtype = PLT_ELEM_TRI6; break;
        default:
            assert(false);
	}

	// write the header
	m_ar.BeginChunk(PLT_DOMAIN_HDR);
	{
		m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
		m_ar.WriteChunk(PLT_DOM_PART_ID  ,   mid);
		m_ar.WriteChunk(PLT_DOM_ELEMS    ,    NE);
	}
	m_ar.EndChunk();

	// write the element list
	int n[9];
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST);
	{
		for (i=0; i<NE; ++i)
		{
			FEShellElement& el = dom.Element(i);
			n[0] = el.GetID();
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j];
			m_ar.WriteChunk(PLT_ELEMENT, n, ne+1);
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteBeamDomain(FEBeamDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	assert(mid > 0);

	int i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 2;
	int dtype = PLT_ELEM_LINE2;

	// write the header
	m_ar.BeginChunk(PLT_DOMAIN_HDR);
	{
		m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
		m_ar.WriteChunk(PLT_DOM_PART_ID  ,   mid);
		m_ar.WriteChunk(PLT_DOM_ELEMS    ,    NE);
	}
	m_ar.EndChunk();

	// write the element list
	int n[5];
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST);
	{
		for (i=0; i<NE; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			n[0] = el.GetID();
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j];
			m_ar.WriteChunk(PLT_ELEMENT, n, ne+1);
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDiscreteDomain(FEDiscreteDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	assert(mid > 0);

	int i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 2;
	int dtype = PLT_ELEM_LINE2;

	// write the header
	m_ar.BeginChunk(PLT_DOMAIN_HDR);
	{
		m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
		m_ar.WriteChunk(PLT_DOM_PART_ID  ,   mid);
		m_ar.WriteChunk(PLT_DOM_ELEMS    ,    NE);
	}
	m_ar.EndChunk();

	// write the element list
	int n[5];
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST);
	{
		for (i=0; i<NE; ++i)
		{
			FEElement& el = dom.ElementRef(i);
			n[0] = el.GetID();
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j];
			m_ar.WriteChunk(PLT_ELEMENT, n, ne+1);
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDomain2D(FEDomain2D& dom)
{
    int mid = dom.GetMaterial()->GetID();
    assert(mid > 0);
    int etype = dom.GetElementType();
    
    int i, j;
    int NE = dom.Elements();
    
    // figure out element type
    int ne = 0;
    int dtype = 0;
    switch (etype)
    {
        case FE2D_TRI3G1 : ne = 3; dtype = PLT_ELEM_TRI; break;
        case FE2D_TRI6G3 : ne = 6; dtype = PLT_ELEM_TRI6; break;
        case FE2D_QUAD4G4: ne = 4; dtype = PLT_ELEM_QUAD; break;
        case FE2D_QUAD8G9: ne = 8; dtype = PLT_ELEM_QUAD8; break;
        case FE2D_QUAD9G9: ne = 9; dtype = PLT_ELEM_QUAD9; break;
		default:
			assert(false);
    }
    
    // write the header
    m_ar.BeginChunk(PLT_DOMAIN_HDR);
    {
        m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
        m_ar.WriteChunk(PLT_DOM_PART_ID  ,   mid);
        m_ar.WriteChunk(PLT_DOM_ELEMS    ,    NE);
    }
    m_ar.EndChunk();
    
    // write the element list
    int n[10];
    m_ar.BeginChunk(PLT_DOM_ELEM_LIST);
    {
        for (i=0; i<NE; ++i)
        {
            FEElement2D& el = dom.Element(i);
            n[0] = el.GetID();
            for (j=0; j<ne; ++j) n[j+1] = el.m_node[j];
            m_ar.WriteChunk(PLT_ELEMENT, n, ne+1);
        }
    }
    m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::BuildSurfaceTable()
{
	FEModel& fem = *GetFEModel();

	FEMesh& mesh = fem.GetMesh();
	m_Surf.clear();
	for (int ns = 0; ns < mesh.Surfaces(); ++ns)
	{
		FESurface& s = mesh.Surface(ns);
		int NF = s.Elements();

		// find the max nodes
		int maxNodes = 0;
		for (int i = 0; i < NF; ++i)
		{
			FESurfaceElement& el = s.Element(i);
			if (el.Nodes() > maxNodes) maxNodes = el.Nodes();
		}

		Surface surf;
		surf.maxNodes = maxNodes;
		surf.surf = &s;
		m_Surf.push_back(surf);
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSurfaceSection(FEMesh& m)
{
	for (int ns = 0; ns<m.Surfaces(); ++ns)
	{
		Surface& surf = m_Surf[ns];
		FESurface& s = *surf.surf;
		int NF = s.Elements();
		int maxNodes = surf.maxNodes;

		m_ar.BeginChunk(PLT_SURFACE);
		{
			m_ar.BeginChunk(PLT_SURFACE_HDR);
			{
				int sid = ns+1;
				m_ar.WriteChunk(PLT_SURFACE_ID, sid);
				m_ar.WriteChunk(PLT_SURFACE_FACES, NF);
				m_ar.WriteChunk(PLT_SURFACE_NAME, s.GetName());
				m_ar.WriteChunk(PLT_SURFACE_MAX_FACET_NODES, maxNodes);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_FACE_LIST);
			{
				int n[FEElement::MAX_NODES + 2];
				for (int i=0; i<NF; ++i)
				{
					FESurfaceElement& f = s.Element(i);
					int nf = f.Nodes();
					n[0] = i+1;
					n[1] = nf;
					for (int i=0; i<nf; ++i) n[i+2] = f.m_node[i];
					m_ar.WriteChunk(PLT_FACE, n, maxNodes + 2);
				}
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeSetSection(FEMesh& m)
{
	for (int ns = 0; ns < m.NodeSets(); ++ns)
	{
		FENodeSet& l = *m.NodeSet(ns);
		int nodes = l.Size();
		m_ar.BeginChunk(PLT_NODESET);
		{
			m_ar.BeginChunk(PLT_NODESET_HDR);
			{
				int nid = ns+1;
				m_ar.WriteChunk(PLT_NODESET_ID, nid);
				m_ar.WriteChunk(PLT_NODESET_SIZE, nodes);
				m_ar.WriteChunk(PLT_NODESET_NAME, l.GetName());
			}
			m_ar.EndChunk();

			std::vector<int> nodeList(nodes);
			for (int i = 0; i < nodes; ++i) nodeList[i] = l[i];
			m_ar.WriteChunk(PLT_NODESET_LIST, nodeList);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteElementSetSection(FEMesh& m)
{
	for (int n = 0; n < m.ElementSets(); ++n)
	{
		FEElementSet& l = m.ElementSet(n);
		int elems = l.Elements();
		m_ar.BeginChunk(PLT_ELEMENTSET);
		{
			m_ar.BeginChunk(PLT_ELEMENTSET_HDR);
			{
				int nid = n + 1;
				m_ar.WriteChunk(PLT_ELEMENTSET_ID, nid);
				m_ar.WriteChunk(PLT_ELEMENTSET_SIZE, elems);
				m_ar.WriteChunk(PLT_ELEMENTSET_NAME, l.GetName());
			}
			m_ar.EndChunk();

			std::vector<int> elemList(elems);
			for (int i = 0; i < elems; ++i) elemList[i] = l[i];
			m_ar.WriteChunk(PLT_ELEMENTSET_LIST, elemList);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteFacetSetSection(FEMesh& m)
{
	for (int ns = 0; ns < m.FacetSets(); ++ns)
	{
		FEFacetSet& surf = m.FacetSet(ns);
		int NF = surf.Faces();

		m_ar.BeginChunk(PLT_FACETSET);
		{
			int maxNodes = FEFacetSet::FACET::MAX_NODES;

			m_ar.BeginChunk(PLT_FACETSET_HDR);
			{
				int sid = ns + 1;
				m_ar.WriteChunk(PLT_FACETSET_ID, sid);
				m_ar.WriteChunk(PLT_FACETSET_SIZE, NF);
				m_ar.WriteChunk(PLT_FACETSET_NAME, surf.GetName());
				m_ar.WriteChunk(PLT_FACETSET_MAXNODES, maxNodes);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_FACETSET_LIST);
			{
				int n[FEElement::MAX_NODES + 2];
				for (int i = 0; i < NF; ++i)
				{
					FEFacetSet::FACET& f = surf.Face(i);
					int nf = f.ntype;	// this is the type and the nr. of nodes!
					n[0] = i + 1;
					n[1] = nf;
					for (int i = 0; i < nf; ++i) n[i + 2] = f.node[i];
					m_ar.WriteChunk(PLT_FACET, n, maxNodes + 2);
				}
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WritePartsSection(FEModel& fem)
{
	int NMAT = fem.Materials();
	for (int i=0; i<NMAT; ++i)
	{
		FEMaterial* pm = fem.GetMaterial(i);
		m_ar.BeginChunk(PLT_PART);
		{
			unsigned int nid = (unsigned int) pm->GetID();
			char szname[STR_SIZE] = {0};

			// Make sure that the material name fits in the buffer
			std::string name = pm->GetName();
			const char* sz = name.c_str();
			int l = (int)strlen(sz);
			if (l >= STR_SIZE) l = STR_SIZE - 1;
			strncpy(szname, sz, l);

			// write the material data
			m_ar.WriteChunk(PLT_PART_ID, nid);
			m_ar.WriteChunk(PLT_PART_NAME, szname, STR_SIZE);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteObjectsSection()
{
	for (int i = 0; i < m_Points.size(); ++i)
	{
		PointObject* po = GetPointObject(i);
		m_ar.BeginChunk(PLT_POINT_OBJECT);
		{
			WriteObject(po);
		}
		m_ar.EndChunk();
	}

	for (int i = 0; i < m_Lines.size(); ++i)
	{
		LineObject* po = GetLineObject(i);
		m_ar.BeginChunk(PLT_LINE_OBJECT);
		{
			WriteObject(po);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteObject(PlotObject* po)
{
	m_ar.WriteChunk(PLT_OBJECT_ID, po->m_id);
	m_ar.WriteChunk(PLT_OBJECT_NAME, po->m_name.c_str());
	m_ar.WriteChunk(PLT_OBJECT_TAG, po->m_tag);

	vec3d& r = po->m_pos;
	float f[3] = { (float)r.x, (float)r.y, (float)r.z };
	m_ar.WriteChunk(PLT_OBJECT_POS, f, 3);

	quatd q = po->m_rot;
	float a[4] = { (float)q.x, (float)q.y, (float)q.z, (float)q.w };
	m_ar.WriteChunk(PLT_OBJECT_ROT, a, 4);

	list<DICTIONARY_ITEM>::iterator it = po->m_data.begin();
	for (int j = 0; j < po->m_data.size(); ++j, ++it)
	{
		m_ar.BeginChunk(PLT_OBJECT_DATA);
		{
			WriteDictionaryItem(*it);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Write(float ftime, int flag)
{
	FEModel& fem = *GetFEModel();

	// compress these sections if requested
	m_ar.SetCompression(m_ncompress);
	m_ar.BeginChunk(PLT_STATE);
	{
		// state header
		m_ar.BeginChunk(PLT_STATE_HEADER);
		{
			m_ar.WriteChunk(PLT_STATE_HDR_TIME, ftime);
			m_ar.WriteChunk(PLT_STATE_STATUS, flag);
		}
		m_ar.EndChunk();

		// write the state flags of the mesh
		m_ar.BeginChunk(PLT_MESH_STATE);
		{
			WriteMeshState(fem.GetMesh());
		}
		m_ar.EndChunk();

		m_ar.BeginChunk(PLT_STATE_DATA);
		{
			// Global Data
			if (!m_dic.m_Glob.empty())
			{
				m_ar.BeginChunk(PLT_GLOBAL_DATA);
				{
					WriteGlobalData(fem);
				}
				m_ar.EndChunk();
			}

			// Node Data
			if (!m_dic.m_Node.empty())
			{
				m_ar.BeginChunk(PLT_NODE_DATA);
				{
					WriteNodeData(fem);
				}
				m_ar.EndChunk();
			}

			// Element Data
			if (!m_dic.m_Elem.empty())
			{
				m_ar.BeginChunk(PLT_ELEMENT_DATA);
				{
					WriteDomainData(fem);
				}
				m_ar.EndChunk();
			}

			// surface data
			if (!m_dic.m_Face.empty())
			{
				m_ar.BeginChunk(PLT_FACE_DATA);
				{
					WriteSurfaceData(fem);
				}
				m_ar.EndChunk();
			}
		}
		m_ar.EndChunk();

		if (m_Points.size() > 0)
		{
			m_ar.BeginChunk(PLT_OBJECTS_STATE);
			{
				WriteObjectsState();
			}
			m_ar.EndChunk();
		}
	}
	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteGlobalData(FEModel& fem)
{

}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeData(FEModel& fem)
{
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Node.begin();
	for (int i=0; i<(int) m_dic.m_Node.size(); ++i, ++it)
	{
		m_ar.BeginChunk(PLT_STATE_VARIABLE);
		{
			unsigned int nid = i+1;
			m_ar.WriteChunk(PLT_STATE_VAR_ID, nid);
			m_ar.BeginChunk(PLT_STATE_VAR_DATA);
			{
				if (it->m_psave) WriteNodeDataField(fem, it->m_psave);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDomainData(FEModel& fem)
{
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Elem.begin();
	for (int i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it)
	{
		m_ar.BeginChunk(PLT_STATE_VARIABLE);
		{
			unsigned int nid = i+1;
			m_ar.WriteChunk(PLT_STATE_VAR_ID, nid);
			m_ar.BeginChunk(PLT_STATE_VAR_DATA);
			{
				if (it->m_psave) WriteDomainDataField(fem, it->m_psave);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSurfaceData(FEModel& fem)
{
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Face.begin();
	for (int i=0; i<(int) m_dic.m_Face.size(); ++i, ++it)
	{
		m_ar.BeginChunk(PLT_STATE_VARIABLE);
		{
			unsigned int nid = i+1;
			m_ar.WriteChunk(PLT_STATE_VAR_ID, nid);
			m_ar.BeginChunk(PLT_STATE_VAR_DATA);
			{
				if (it->m_psave) WriteSurfaceDataField(fem, it->m_psave);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeDataField(FEModel &fem, FEPlotData* pd)
{
	// loop over all node sets
	// right now there is only one, namely the node set of all mesh nodes
	// so we just pass the mesh
	int ndata = pd->VarSize(pd->DataType());

	int N = fem.GetMesh().Nodes();
	FEDataStream a; a.reserve(ndata*N);
	if (pd->Save(fem.GetMesh(), a))
	{
		// pad mismatches
		assert(a.size() == N*ndata);
		if (a.size() != N * ndata) a.resize(N*ndata, 0.f);
		m_ar.WriteData(0, a.data());
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSurfaceDataField(FEModel& fem, FEPlotData* pd)
{
	// get the domain name (if any)
	string domName;
	const char* szdom = pd->GetDomainName();
	if (szdom) domName = szdom;

	// loop over all surfaces
	FEMesh& m = fem.GetMesh();
	int NS = m.Surfaces();
	for (int i = 0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);

		if (domName.empty() || (domName == S.GetName()))
		{
			Surface& surf = m_Surf[i];
			assert(surf.surf == &S);

			// Determine data size.
			// Note that for the FMT_MULT case we are 
			// assuming 9 data entries per facet
			// regardless of the nr of nodes a facet really has
			// this is because for surfaces, all elements are not
			// necessarily of the same type
			// TODO: Fix the assumption of the FMT_MULT
			int datasize = pd->VarSize(pd->DataType());
			int nsize = datasize;
			switch (pd->StorageFormat())
			{
			case FMT_NODE: nsize *= S.Nodes(); break;
			case FMT_ITEM: nsize *= S.Elements(); break;
			case FMT_MULT: nsize *= surf.maxNodes * S.Elements(); break;
			case FMT_REGION:
				// one value per surface so nsize remains unchanged
				break;
			default:
				assert(false);
			}

			// save data
			FEDataStream a; a.reserve(nsize);
			if (pd->Save(S, a))
			{
				// in FEBio 3.0, the data streams are assumed to have no padding, but for now we still need to pad 
				// the data stream before we write it to the file
				if (a.size() == nsize)
				{
					// assumed padding is already there, or not needed
					m_ar.WriteData(i + 1, a.data());
				}
				else
				{
					// this is only needed for FMT_MULT storage
					assert(pd->StorageFormat() == FMT_MULT);

					// add padding
					const int M = surf.maxNodes;
					int m = 0;
					FEDataStream b; b.assign(nsize, 0.f);
					for (int n = 0; n < S.Elements(); ++n)
					{
						FESurfaceElement& el = S.Element(n);
						int ne = el.Nodes();
						for (int j = 0; j < ne; ++j)
						{
							for (int k = 0; k < datasize; ++k) b[n * M * datasize + j * datasize + k] = a[m++];
						}
					}

					// write the padded data
					m_ar.WriteData(i + 1, b.data());
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDomainDataField(FEModel &fem, FEPlotData* pd)
{
	FEMesh& m = fem.GetMesh();
	int ND = m.Domains();

	// if the item list is empty, store all domains
	vector<int> item = pd->GetItemList();
	if (item.empty())
	{
		for (int i = 0; i<ND; ++i) item.push_back(i);
	}

	// allow plot data to prepare for save
	if (pd->PreSave() == false)
	{
		assert(false);
		return;
	}

	// get the domain name (if any)
	string domName;
	const char* szdom = pd->GetDomainName();
	if (szdom) domName = szdom;

	// loop over all domains in the item list
	for (int i = 0; i<ND; ++i)
	{
		// get the domain
		FEDomain& D = m.Domain(item[i]);

		if (domName.empty() || (D.GetName() == domName))
		{
			// calculate the size of the data vector
			int nsize = pd->VarSize(pd->DataType());
			switch (pd->StorageFormat())
			{
			case FMT_NODE: nsize *= D.Nodes(); break;
			case FMT_ITEM: nsize *= D.Elements(); break;
			case FMT_MULT:
			{
				// since all elements have the same type within a domain
				// we just grab the number of nodes of the first element 
				// to figure out how much storage we need
				FEElement& e = D.ElementRef(0);
				int n = e.Nodes();
				nsize *= n * D.Elements();
			}
			break;
			case FMT_REGION:
				// one value for this domain so nsize remains unchanged
				break;
			default:
				assert(false);
			}
			assert(nsize > 0);

			// fill data vector and save
			FEDataStream a;
			a.reserve(nsize);
			if (pd->Save(D, a))
			{
				assert(a.size() == nsize);
				m_ar.WriteData(item[i] + 1, a.data());
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Append(const char *szfile)
{
	// try to open the file
	if (m_ar.Open(szfile) == false) return false;

	FEModel* fem = GetFEModel();
	FEPlotDataStore& pltData = fem->GetPlotDataStore();
	SetCompression(pltData.GetPlotCompression());

	// add plot variables
	for (int n = 0; n < pltData.PlotVariables(); ++n)
	{
		FEPlotVariable& vi = pltData.GetPlotVariable(n);
		const std::string& varName = vi.Name();
		const std::string& domName = vi.DomainName();

		// add the plot output variable
		if (AddVariable(varName.c_str(), vi.m_item, domName.c_str()) == false)
		{
			feLog("FATAL ERROR: Output variable \"%s\" is not defined\n", varName.c_str());
			throw "FATAL ERROR";
		}
	}

	// NOTE: Reading the dictionary rebuilds the plot variables too, but
	//       there is not enough data in the dictionary to do that correctly. 
	//       So, we should probably rely on the data store, which gets serialized to the dump file.
	// open the root element
	bool bok = true;
/*	m_ar.OpenChunk();
	unsigned int nid = m_ar.GetChunkID();
	if (nid != PLT_ROOT) return false;

	bok = false;
	while (m_ar.OpenChunk() == IO_OK)
	{
		nid = m_ar.GetChunkID();
		if (nid == PLT_DICTIONARY)
		{
			// read the dictionary
			bok = ReadDictionary();
			break;
		}
		m_ar.CloseChunk();
	}
*/
	// close it again ...
	m_ar.Close();

	// rebuild the surface table
	BuildSurfaceTable();

	// ... and open for appending
	if (bok) return m_ar.Append(szfile);

	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::ReadDictionary()
{
	m_dic.Clear();

	while (m_ar.OpenChunk() == IO_OK)
	{
		unsigned int nid = m_ar.GetChunkID();
		switch (nid)
		{
		case PLT_DIC_GLOBAL: assert(false); return false;
		case PLT_DIC_NODAL  : ReadDicList(); break;
		case PLT_DIC_DOMAIN : ReadDicList(); break;
		case PLT_DIC_SURFACE: ReadDicList(); break;
		default:
			assert(false);
			return false;
		}
		m_ar.CloseChunk();
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::ReadDicList()
{
	vector<int> l; // empty item list
	while (m_ar.OpenChunk() == IO_OK)
	{
		unsigned int nid = m_ar.GetChunkID();
		if (nid == PLT_DIC_ITEM)
		{
			while (m_ar.OpenChunk() == IO_OK)
			{
				unsigned int nid = m_ar.GetChunkID();
				if (nid == PLT_DIC_ITEM_NAME)
				{
					char sz[STR_SIZE];
					m_ar.read(sz, STR_SIZE);
					AddVariable(sz, l);
				}
				m_ar.CloseChunk();
			}
		}
		else return false;
		m_ar.CloseChunk();
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteMeshState(FEMesh& mesh)
{
	vector<unsigned int> flags;
	flags.reserve(mesh.Elements());
	int NDOM = mesh.Domains();
	for (int i = 0; i < NDOM; ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			unsigned int status = el.status();
			flags.push_back(status);
		}
	}

	m_ar.WriteChunk(PLT_ELEMENT_STATE, flags);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteObjectsState()
{
	for (int i = 0; i < PointObjects(); ++i)
	{
		PointObject* po = m_Points[i];
		m_ar.BeginChunk(PLT_POINT_OBJECT);
		{
			m_ar.WriteChunk(PLT_OBJECT_ID, po->m_id);

			vec3d r = po->m_pos;
			float f[3] = { (float)r.x, (float)r.y, (float)r.z };
			m_ar.WriteChunk(PLT_OBJECT_POS, f, 3);

			quatd q = po->m_rot;
			float a[4] = { (float)q.x, (float)q.y, (float)q.z, (float)q.w };
			m_ar.WriteChunk(PLT_OBJECT_ROT, a, 4);

			r = po->m_r;
			float c[3] = { (float)r.x, (float)r.y, (float)r.z };
			m_ar.WriteChunk(PLT_POINT_COORD, c, 3);

			m_ar.BeginChunk(PLT_OBJECT_DATA);
			{
				WriteObjectData(po);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}

	for (int i = 0; i < m_Lines.size(); ++i)
	{
		LineObject* po = GetLineObject(i);
		m_ar.BeginChunk(PLT_LINE_OBJECT);
		{
			m_ar.WriteChunk(PLT_OBJECT_ID, po->m_id);

			vec3d r = po->m_pos;
			float f[3] = { (float)r.x, (float)r.y, (float)r.z };
			m_ar.WriteChunk(PLT_OBJECT_POS, f, 3);

			quatd q = po->m_rot;
			float a[4] = { (float)q.x, (float)q.y, (float)q.z, (float)q.w };
			m_ar.WriteChunk(PLT_OBJECT_ROT, a, 4);

			vec3d r1 = po->m_r1;
			vec3d r2 = po->m_r2;
			float c[6] = { (float)r1.x, (float)r1.y, (float)r1.z, (float)r2.x, (float)r2.y, (float)r2.z };
			m_ar.WriteChunk(PLT_LINE_COORDS, c, 6);

			m_ar.BeginChunk(PLT_OBJECT_DATA);
			{
				WriteObjectData(po);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteObjectData(PlotObject* po)
{
	list<DICTIONARY_ITEM>::iterator it = po->m_data.begin();
	for (int j = 0; j < po->m_data.size(); ++j, ++it)
	{
		assert(it->m_psave);

		FEPlotObjectData* pd = dynamic_cast<FEPlotObjectData*>(it->m_psave); assert(pd);

		FEDataStream a;
		if (pd->Save(po, a))
		{
			m_ar.WriteData(j, a.data());
		}
	}
}
