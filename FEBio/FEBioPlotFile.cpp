#include "StdAfx.h"
#include "FEBioPlotFile.h"
#include "fem.h"

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddGlobalVariable(unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddNodalVariable(unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	strcpy(it.m_szname, szname);
	m_Node.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddSolidVariable(unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	strcpy(it.m_szname, szname);
	m_Elem.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddShellVariable(unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	strcpy(it.m_szname, szname);
	m_Shell.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddBeamVariable(unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	strcpy(it.m_szname, szname);
	m_Beam.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::Save(Archive &ar)
{
	// store global variables
	if (m_Glob.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Glob.begin();
		for (int i=0; i<(int) m_Glob.size(); ++i, ++pi)
		{
			ar << pi->m_ntype;
			ar.write(pi->m_szname, sizeof(char), DI_NAME_SIZE);
		}
	}

	// store nodal variables
	if (m_Node.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Node.begin();
		for (int i=0; i<(int) m_Node.size(); ++i, ++pi)
		{
			ar << pi->m_ntype;
			ar.write(pi->m_szname, sizeof(char), DI_NAME_SIZE);
		}
	}

	// store solid variables
	if (m_Elem.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Elem.begin();
		for (int i=0; i<(int) m_Elem.size(); ++i, ++pi)
		{
			ar << pi->m_ntype;
			ar.write(pi->m_szname, sizeof(char), DI_NAME_SIZE);
		}
	}

	// store shell variables
	if (m_Shell.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Shell.begin();
		for (int i=0; i<(int) m_Shell.size(); ++i, ++pi)
		{
			ar << pi->m_ntype;
			ar.write(pi->m_szname, sizeof(char), DI_NAME_SIZE);
		}
	}

	// store beam variables
	if (m_Beam.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Beam.begin();
		for (int i=0; i<(int) m_Beam.size(); ++i, ++pi)
		{
			ar << pi->m_ntype;
			ar.write(pi->m_szname, sizeof(char), DI_NAME_SIZE);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

FEBioPlotFile::FEBioPlotFile(void)
{
	// set the default export data
	m_dic.AddNodalVariable(VEC3F, "Displacement");
}

FEBioPlotFile::~FEBioPlotFile(void)
{
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Open(FEM &fem, const char *szfile)
{
	// get the mesh
	FEMesh& m = fem.m_mesh;

	// setup the header
	m_hdr.nversion = 0;
	m_hdr.nnodes = m.Nodes();
	m_hdr.n3d    = m.SolidElements();
	m_hdr.n2d    = m.ShellElements();
	m_hdr.n1d    = fem.m_DE.size();

	m_hdr.nglv = m_dic.m_Glob.size();
	m_hdr.nnv  = m_dic.m_Node.size();
	m_hdr.nv3d = m_dic.m_Elem.size();
	m_hdr.nv2d = m_dic.m_Shell.size();
	m_hdr.nv1d = m_dic.m_Beam.size();

	// open the archive
	if (m_ar.Create(szfile) == false) return false;

	// --- save the header file ---
	m_ar.write(&m_hdr, sizeof(HEADER), 1);

	// --- save the dictionary ---
	m_dic.Save(m_ar);

	// --- save the geometry ---
	int i, j, N;
	
	// write the material coordinates
	float xf[3];
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		xf[0] = (float) node.m_r0.x;
		xf[1] = (float) node.m_r0.y;
		xf[2] = (float) node.m_r0.z;

		m_ar.write(xf, sizeof(float), 3);
	}

	// write the connectivity and material number
	// Note that we increment all numbers by 1 since
	// the plot database expects 1-based arrays
	int n[9];

	int nid = 1;

	// write solid element data
	// note that we reindex all elements so that the ID
	// corresponds to the nr in the plot file
	for (i=0; i<m.SolidElements(); ++i)
	{
		FESolidElement& el = m.SolidElement(i);

		el.m_nID = nid++;

		N = el.Nodes();
		switch (el.Type())
		{
		case FE_HEX:
		case FE_RIHEX:
		case FE_UDGHEX:
			for (j=0; j<N; ++j) n[j] = el.m_node[j]+1;
			break;
		case FE_PENTA:
			n[0] = el.m_node[0]+1;
			n[1] = el.m_node[1]+1;
			n[2] = el.m_node[2]+1;
			n[3] = el.m_node[2]+1;
			n[4] = el.m_node[3]+1;
			n[5] = el.m_node[4]+1;
			n[6] = el.m_node[5]+1;
			n[7] = el.m_node[5]+1;
			break;
		case FE_TET:
			n[0] = el.m_node[0]+1;
			n[1] = el.m_node[1]+1;
			n[2] = el.m_node[2]+1;
			n[3] = el.m_node[2]+1;
			n[4] = n[5] = n[6] = n[7] = el.m_node[3]+1;
			break;
		}

		n[8] = el.GetMatID()+1;

		m_ar.write(n, sizeof(int), 9);
	}

	// write shell element data
	for (i=0; i<m.ShellElements(); ++i)
	{
		FEShellElement& el = m.ShellElement(i);

		el.m_nID = nid++;

		N = el.Nodes();
		switch (el.Type())
		{
		case FE_SHELL_QUAD:
			n[0] = el.m_node[0]+1;
			n[1] = el.m_node[1]+1;
			n[2] = el.m_node[2]+1;
			n[3] = el.m_node[3]+1;
			break;
		case FE_SHELL_TRI:
			n[0] = el.m_node[0]+1;
			n[1] = el.m_node[1]+1;
			n[2] = el.m_node[2]+1;
			n[3] = el.m_node[2]+1;
			break;
		}

		n[4] = el.GetMatID()+1;

		m_ar.write(n, sizeof(int), 5);
	}

	// write truss element data
	for (i=0; i<m.TrussElements(); ++i)
	{
		FETrussElement& el = m.TrussElement(i);
		el.m_nID = nid++;
		n[0] = el.m_node[0]+1;
		n[1] = el.m_node[1]+1;
		n[2] = 0;
		n[3] = 0;
		n[4] = 0;
		n[5] = el.GetMatID()+1;

		m_ar.write(n, sizeof(int), 6);
	}

	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Append(FEM &fem, const char *szfile)
{
	return false;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Write(FEM &fem)
{
	// make sure the archive is opened
	if (!m_ar.IsValid()) return false;

	// store the fem pointer
	m_pfem = &fem;

	// get the mesh
	FEMesh& mesh = fem.m_mesh;

	// save the time stamp
	float time = (float) fem.m_ftime;
	m_ar << time;

	// save the global variables
	if (m_dic.m_Glob.size() > 0)
	{
	}

	// save the nodal variables
	if (m_dic.m_Node.size() > 0)
	{
		write_displacements();
	}

	// save the solid variables
	if (m_dic.m_Elem.size() > 0)
	{
	}

	// save the shell variables
	if (m_dic.m_Shell.size() > 0)
	{
	}

	// save the beam variables
	if (m_dic.m_Beam.size() > 0)
	{
	}

	return false;
}
