#include "stdafx.h"
#include "FEBioPlotFile.h"
#include "fem.h"
#include "FETransverselyIsotropic.h"

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddGlobalVariable(FEPlotData* ps, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ps->DataType();
	it.m_nfmt  = ps->StorageFormat();
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddMaterialVariable(FEPlotData* ps, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ps->DataType();
	it.m_nfmt  = ps->StorageFormat();
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddNodalVariable(FEPlotData* ps, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ps->DataType();
	it.m_nfmt  = ps->StorageFormat();
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Node.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddElementVariable(FEPlotData* ps, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ps->DataType();
	it.m_nfmt  = ps->StorageFormat();
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Elem.push_back(it);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

FEBioPlotFile::FEBioPlotFile(void)
{
}

FEBioPlotFile::~FEBioPlotFile(void)
{
	m_ar.Close();

	int i;
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Glob.begin();
	for (i=0; i<(int) m_dic.m_Glob.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Mat.begin();
	for (i=0; i<(int) m_dic.m_Mat.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Node.begin();
	for (i=0; i<(int) m_dic.m_Node.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Elem.begin();
	for (i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Face.begin();
	for (i=0; i<(int) m_dic.m_Face.size(); ++i, ++it) delete it->m_psave;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Open(FEM &fem, const char *szfile)
{
	// open the archive
	if (m_ar.Create(szfile) == false) return false;

	// write the root element
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

		// --- save the materials
		m_ar.BeginChunk(PLT_MATERIALS);
		{
			if (WriteMaterials(fem) == false) return false;
		}
		m_ar.EndChunk();

		// --- save the geometry ---
		m_ar.BeginChunk(PLT_GEOMETRY);
		{
			if (WriteGeometry(fem) == false) return false;
		}
		m_ar.EndChunk();
	}
	// Don't call EndChunk yet since we still 
	// need to write state data
//	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Append(FEM &fem, const char *szfile)
{
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteHeader(FEM& fem)
{
	// setup the header
	unsigned int nversion = PLT_VERSION;

	// output header
	m_ar.WriteChunk(PLT_HDR_VERSION, nversion);

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteDictionary(FEM& fem)
{
	// First we build the dictionary
	// get the mesh
	FEMesh& m = fem.m_mesh;

	int nmode = fem.m_pStep->m_nModule;
	int ntype = fem.m_pStep->m_nanalysis;

	// setup the dictionary
	m_dic.AddNodalVariable  (new FEPlotNodeDisplacement, "Displacement");
	m_dic.AddElementVariable(new FEPlotElementStress   , "Stress");

	// store dynamic analysis data
	if ((ntype == FE_DYNAMIC) || (nmode == FE_POROELASTIC)) m_dic.AddNodalVariable(new FEPlotNodeVelocity    , "Velocity");
	if (ntype == FE_DYNAMIC) m_dic.AddNodalVariable(new FEPlotNodeAcceleration, "Acceleration");

	// store contact data
	if (fem.m_CI.size() > 0)
	{
		m_dic.AddNodalVariable(new FEPlotContactGap     , "Contact gap");
		m_dic.AddNodalVariable(new FEPlotContactTraction, "Contact traction");
	}

	// store poro data
	if (nmode == FE_POROELASTIC)
	{
		m_dic.AddNodalVariable  (new FEPlotFluidPressure, "Fluid Pressure");
		m_dic.AddElementVariable(new FEPlotFluidFlux    , "Fluid Flux");
	}

	// if any material is trans-iso we store material fibers and strain
	int ntiso = 0;
	for (int i=0; i<fem.Materials(); ++i)
	{
		FEElasticMaterial* pm = fem.GetElasticMaterial(i);
		if (dynamic_cast<FETransverselyIsotropic*>(pm)) ntiso++;
	}
	if (ntiso)
	{
		m_dic.AddElementVariable(new FEPlotFiberVector, "Fiber vector");
	}

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

	// store material variables
	if (!m_dic.m_Mat.empty())
	{
		m_ar.BeginChunk(PLT_DIC_MATERIAL);
		{
			WriteDicList(m_dic.m_Mat);
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
			m_ar.WriteChunk(PLT_DIC_ITEM_TYPE, pi->m_ntype);
			m_ar.WriteChunk(PLT_DIC_ITEM_FMT , pi->m_nfmt);
			m_ar.WriteChunk(PLT_DIC_ITEM_NAME, pi->m_szname, STR_SIZE);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteMaterials(FEM& fem)
{
	int NMAT = fem.Materials();
	for (int i=0; i<NMAT; ++i)
	{
		FEMaterial* pm = fem.GetMaterial(i);
		m_ar.BeginChunk(PLT_DIC_MATERIAL);
		{
			unsigned int nid = (unsigned int) pm->GetID();
			char szname[STR_SIZE] = {0};
			strcpy(szname, pm->GetName());
			m_ar.WriteChunk(PLT_MAT_ID, nid);
			m_ar.WriteChunk(PLT_MAT_NAME, szname, STR_SIZE);
		}
		m_ar.EndChunk();
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteGeometry(FEM& fem)
{
	// get the mesh
	FEMesh& m = fem.m_mesh;

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
		m_ar.BeginChunk(PLT_SURFACE_SECTION);
		{
			WriteSurfaceSection(m);
		}
		m_ar.EndChunk();
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeSection(FEMesh& m)
{
	// write the material coordinates
	int NN = m.Nodes();
	vector<float> X(3*NN);
	for (int i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		X[3*i  ] = (float) node.m_r0.x;
		X[3*i+1] = (float) node.m_r0.y;
		X[3*i+2] = (float) node.m_r0.z;
	}

	m_ar.WriteChunk(PLT_NODE_COORDS, X);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDomainSection(FEMesh& m)
{
	// write all domains
	for (int nd = 0; nd<m.Domains(); ++nd)
	{
		m_ar.BeginChunk(PLT_DOMAIN);
		{
			FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&m.Domain(nd));
			if (pbd) WriteSolidDomain(*pbd);

			FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
			if (psd) WriteShellDomain(*psd);

			FETrussDomain* ptd = dynamic_cast<FETrussDomain*>(&m.Domain(nd));
			if (ptd) WriteTrussDomain(*ptd);
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSolidDomain(FESolidDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	int etype = dom.GetElementType();

	int n[9], i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 0;
	int dtype = 0;
	switch (etype)
	{
		case FE_HEX:
		case FE_RIHEX:
		case FE_UDGHEX: ne = 8; dtype = PLT_ELEM_HEX; break;
		case FE_PENTA: ne = 6; dtype = PLT_ELEM_PENTA; break;
		case FE_TET:
		case FE_TETG1: ne = 4; dtype = PLT_ELEM_TET; break;
	}

	// write the header
	m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
	m_ar.WriteChunk(PLT_DOM_MAT_ID   ,   mid);

	// write the element list
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST, NE*(ne+1));
	{
		for (i=0; i<NE; ++i)
		{
			FESolidElement& el = dom.Element(i);
			n[0] = el.m_nID;
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j]+1;

			for (j=0; j<=ne; ++j) m_ar << n[j];
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteShellDomain(FEShellDomain& dom)
{
	int mid = dom.GetMaterial()->GetID();
	int etype = dom.GetElementType();

	int n[5], i, j;
	int NE = dom.Elements();

	// figure out element type
	int ne = 0;
	int dtype = 0;
	switch (etype)
	{
		case FE_SHELL_QUAD: ne = 4; dtype = PLT_ELEM_QUAD; break;
		case FE_SHELL_TRI : ne = 3; dtype = PLT_ELEM_TRI; break;
	}

	// write the header
	m_ar.WriteChunk(PLT_DOM_ELEM_TYPE, dtype);
	m_ar.WriteChunk(PLT_DOM_MAT_ID   ,   mid);

	// write the element list
	m_ar.BeginChunk(PLT_DOM_ELEM_LIST, NE*(ne+1));
	{
		for (i=0; i<NE; ++i)
		{
			FEShellElement& el = dom.Element(i);
			n[0] = el.m_nID;
			for (j=0; j<ne; ++j) n[j+1] = el.m_node[j]+1;
			for (j=0; j<=ne; ++j) m_ar << n[j];
		}
	}
	m_ar.EndChunk();
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteTrussDomain(FETrussDomain& dom)
{
	assert(false);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSurfaceSection(FEMesh& m)
{
	assert(false);
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Write(FEM &fem)
{
	// store the fem pointer
	m_pfem = &fem;

	m_ar.BeginChunk(PLT_STATE);
	{
		// state header
		m_ar.BeginChunk(PLT_STATE_HEADER);
		{
			m_ar.WriteChunk(PLT_STATE_HDR_TIME, fem.m_ftime);
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

			// Material Data
			if (!m_dic.m_Mat.empty())
			{
				m_ar.BeginChunk(PLT_MATERIAL_DATA);
				{
					WriteMaterialData(fem);
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
					WriteElementData(fem);
				}
				m_ar.EndChunk();
			}

			// surface data
			if (!m_dic.m_Face.empty())
			{
				m_ar.BeginChunk(PLT_FACE_DATA);
				{
					WriteFaceData(fem);
				}
				m_ar.EndChunk();
			}
		}
		m_ar.EndChunk();
	}
	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteGlobalData(FEM& fem)
{

}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteMaterialData(FEM& fem)
{

}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteNodeData(FEM& fem)
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
				if (it->m_psave) (it->m_psave)->Save(fem, m_ar);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteElementData(FEM& fem)
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
				if (it->m_psave) (it->m_psave)->Save(fem, m_ar);
			}
			m_ar.EndChunk();
		}
		m_ar.EndChunk();
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteFaceData(FEM& fem)
{

}
