#include "stdafx.h"
#include "FEBioPlotFile.h"
#include "fem.h"
#include "FETransverselyIsotropic.h"

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddGlobalVariable(FEPlotData* ps, unsigned int ntype, unsigned int nfmt, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_nfmt  = nfmt;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddMaterialVariable(FEPlotData* ps, unsigned int ntype, unsigned int nfmt, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_nfmt  = nfmt;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::Dictionary::AddNodalVariable(FEPlotData* ps, unsigned int ntype, unsigned int nfmt, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_nfmt  = nfmt;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Node.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddElementVariable(FEPlotData* ps, unsigned int ntype, unsigned int nfmt, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_nfmt  = nfmt;
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
		if (WriteHeader(fem) == false) return false;

		// --- save the dictionary ---
		if (WriteDictionary(fem) == false) return false;

		// --- save the materials
		if (WriteMaterials(fem) == false) return false;

		// --- save the geometry ---
		if (WriteGeometry(fem) == false) return false;
	}
	// Don't call EndChunk yet since we still 
	// need to write state data
//	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteHeader(FEM& fem)
{
	// setup the header
	m_hdr.nversion = PLT_VERSION;

	// output header
	m_ar.BeginChunk(PLT_HEADER);
	{
		m_ar.WriteChunk(PLT_HDR_VERSION, m_hdr.nversion);
	}
	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::WriteDictionary(FEM& fem)
{
	int i;

	// First we build the dictionary

	// get the mesh
	FEMesh& m = fem.m_mesh;

	int nmode = fem.m_pStep->m_nModule;
	int ntype = fem.m_pStep->m_nanalysis;

	// setup the dictionary
	m_dic.AddNodalVariable  (new FEPlotNodeDisplacement,  VEC3F, ITEM_DATA, "Displacement");
	m_dic.AddElementVariable(new FEPlotElementStress   , MAT3FS, ITEM_DATA, "Stress");

	// store dynamic analysis data
	if ((ntype == FE_DYNAMIC) || (nmode == FE_POROELASTIC)) m_dic.AddNodalVariable(new FEPlotNodeVelocity    , VEC3F, ITEM_DATA, "Velocity");
	if (ntype == FE_DYNAMIC) m_dic.AddNodalVariable(new FEPlotNodeAcceleration, VEC3F, ITEM_DATA, "Acceleration");

	// store contact data
	if (fem.m_CI.size() > 0)
	{
		m_dic.AddNodalVariable(new FEPlotContactGap     , FLOAT, ITEM_DATA, "Contact gap");
		m_dic.AddNodalVariable(new FEPlotContactTraction, VEC3F, ITEM_DATA, "Contact traction");
	}

	// store poro data
	if (nmode == FE_POROELASTIC)
	{
		m_dic.AddNodalVariable  (new FEPlotFluidPressure, FLOAT, ITEM_DATA, "Fluid Pressure");
		m_dic.AddElementVariable(new FEPlotFluidFlux    , VEC3F, ITEM_DATA, "Fluid Flux");
	}

	// if any material is trans-iso we store material fibers and strain
	int ntiso = 0;
	for (i=0; i<fem.Materials(); ++i)
	{
		FEElasticMaterial* pm = fem.GetElasticMaterial(i);
		if (dynamic_cast<FETransverselyIsotropic*>(pm)) ntiso++;
	}
	if (ntiso)
	{
		m_dic.AddElementVariable(new FEPlotFiberVector, VEC3F, ITEM_DATA, "Fiber vector");
	}

	// Next, we save the dictionary
	m_ar.BeginChunk(PLT_DICTIONARY);
	{
		// Global variables
		m_ar.BeginChunk(PLT_DIC_GLOBAL);
		{
			WriteDicList(m_dic.m_Glob);
		}
		m_ar.EndChunk();

		// store material variables
		m_ar.BeginChunk(PLT_DIC_MATERIAL);
		{
			WriteDicList(m_dic.m_Mat);
		}
		m_ar.EndChunk();

		// store nodal variables
		m_ar.BeginChunk(PLT_DIC_NODAL);
		{
			WriteDicList(m_dic.m_Node);
		}
		m_ar.EndChunk();

		// store element variables
		m_ar.BeginChunk(PLT_DIC_DOMAIN);
		{
			WriteDicList(m_dic.m_Elem);
		}
		m_ar.EndChunk();

		// store surface data
		m_ar.BeginChunk(PLT_DIC_SURFACE);
		{
			WriteDicList(m_dic.m_Face);
		}
		m_ar.EndChunk();
	}
	m_ar.EndChunk();

	return true;
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteDicList(list<FEBioPlotFile::DICTIONARY_ITEM>& dic)
{
	int N = (int) dic.size();
	m_ar.WriteChunk(PLT_DIC_ENTRIES, N);
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
			unsigned int nid = i+1;
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
	m_ar.BeginChunk(PLT_SURFACE_SECTION);
	{
		WriteSurfaceSection(m);
	}
	m_ar.EndChunk();

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
	int n[9], i, j;
	int NE = dom.Elements();
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		n[0] = el.m_nID;
		switch (el.Type())
		{
		case FE_HEX:
		case FE_RIHEX:
		case FE_UDGHEX:
			for (j=0; j<8; ++j) n[j+1] = el.m_node[j]+1;
			break;
		case FE_PENTA:
			n[1] = el.m_node[0]+1;
			n[2] = el.m_node[1]+1;
			n[3] = el.m_node[2]+1;
			n[4] = el.m_node[2]+1;
			n[5] = el.m_node[3]+1;
			n[6] = el.m_node[4]+1;
			n[7] = el.m_node[5]+1;
			n[8] = el.m_node[5]+1;
			break;
		case FE_TET:
		case FE_TETG1:
			n[1] = el.m_node[0]+1;
			n[2] = el.m_node[1]+1;
			n[3] = el.m_node[2]+1;
			n[4] = el.m_node[2]+1;
			n[5] = n[6] = n[7] = n[8] = el.m_node[3]+1;
			break;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteShellDomain(FEShellDomain& dom)
{
	int n[5];
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FEShellElement& el = dom.Element(i);
		n[0] = el.m_nID;
		switch (el.Type())
		{
		case FE_SHELL_QUAD:
			n[1] = el.m_node[0]+1;
			n[2] = el.m_node[1]+1;
			n[3] = el.m_node[2]+1;
			n[4] = el.m_node[3]+1;
			break;
		case FE_SHELL_TRI:
			n[1] = el.m_node[0]+1;
			n[2] = el.m_node[1]+1;
			n[3] = el.m_node[2]+1;
			n[4] = el.m_node[2]+1;
			break;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteTrussDomain(FETrussDomain& dom)
{
	int n[3];
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FETrussElement& el = dom.Element(i);
		n[0] = el.m_nID;
		n[1] = el.m_node[0]+1;
		n[2] = el.m_node[1]+1;
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteSurfaceSection(FEMesh& m)
{

}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Append(FEM &fem, const char *szfile)
{
	return false;
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

		}
		m_ar.EndChunk();

		m_ar.BeginChunk(PLT_STATE_DATA);
		{
			m_ar.BeginChunk(PLT_GLOBAL_DATA);
			{
				WriteGlobalData(fem);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_MATERIAL_DATA);
			{
				WriteMaterialData(fem);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_NODE_DATA);
			{
				WriteNodeData(fem);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_ELEMENT_DATA);
			{
				WriteElementData(fem);
			}
			m_ar.EndChunk();

			m_ar.BeginChunk(PLT_FACE_DATA);
			{
				WriteFaceData(fem);
			}
			m_ar.EndChunk();
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
	// save the nodal variables
	if (m_dic.m_Node.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator it = m_dic.m_Node.begin();
		for (int i=0; i<(int) m_dic.m_Node.size(); ++i, ++it)
			if (it->m_psave) (it->m_psave)->Save(fem, m_ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteElementData(FEM& fem)
{
	// save the solid variables
	if (m_dic.m_Elem.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator it = m_dic.m_Elem.begin();
		for (int i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it)
			if (it->m_psave) (it->m_psave)->Save(fem, m_ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioPlotFile::WriteFaceData(FEM& fem)
{

}
