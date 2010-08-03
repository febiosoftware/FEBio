#include "stdafx.h"
#include "FEBioPlotFile.h"
#include "fem.h"

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddGlobalVariable(FEPlotData* ps, unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Glob.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddNodalVariable(FEPlotData* ps, unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Node.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddSolidVariable(FEPlotData* ps, unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Elem.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddShellVariable(FEPlotData* ps, unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Shell.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::AddBeamVariable(FEPlotData* ps, unsigned int ntype, const char* szname)
{
	DICTIONARY_ITEM it;
	it.m_ntype = ntype;
	it.m_psave = ps;
	strcpy(it.m_szname, szname);
	m_Beam.push_back(it);
}

//-----------------------------------------------------------------------------

void FEBioPlotFile::Dictionary::Save(FILE* fp)
{
	// store global variables
	if (m_Glob.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Glob.begin();
		for (int i=0; i<(int) m_Glob.size(); ++i, ++pi)
		{
			fwrite(&pi->m_ntype, sizeof(pi->m_ntype), 1, fp);
			fwrite(pi->m_szname, sizeof(char), DI_NAME_SIZE, fp);
		}
	}

	// store nodal variables
	if (m_Node.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Node.begin();
		for (int i=0; i<(int) m_Node.size(); ++i, ++pi)
		{
			fwrite(&pi->m_ntype, sizeof(pi->m_ntype), 1, fp);
			fwrite(pi->m_szname, sizeof(char), DI_NAME_SIZE, fp);
		}
	}

	// store solid variables
	if (m_Elem.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Elem.begin();
		for (int i=0; i<(int) m_Elem.size(); ++i, ++pi)
		{
			fwrite(&pi->m_ntype, sizeof(pi->m_ntype), 1, fp);
			fwrite(pi->m_szname, sizeof(char), DI_NAME_SIZE, fp);
		}
	}

	// store shell variables
	if (m_Shell.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Shell.begin();
		for (int i=0; i<(int) m_Shell.size(); ++i, ++pi)
		{
			fwrite(&pi->m_ntype, sizeof(pi->m_ntype), 1, fp);
			fwrite(pi->m_szname, sizeof(char), DI_NAME_SIZE, fp);
		}
	}

	// store beam variables
	if (m_Beam.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator pi = m_Beam.begin();
		for (int i=0; i<(int) m_Beam.size(); ++i, ++pi)
		{
			fwrite(&pi->m_ntype, sizeof(pi->m_ntype), 1, fp);
			fwrite(pi->m_szname, sizeof(char), DI_NAME_SIZE, fp);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

FEBioPlotFile::FEBioPlotFile(void)
{
	m_fp = 0;
}

FEBioPlotFile::~FEBioPlotFile(void)
{
	if (m_fp) fclose(m_fp);

	int i;
	list<DICTIONARY_ITEM>::iterator it = m_dic.m_Glob.begin();
	for (i=0; i<(int) m_dic.m_Glob.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Node.begin();
	for (i=0; i<(int) m_dic.m_Node.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Elem.begin();
	for (i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Shell.begin();
	for (i=0; i<(int) m_dic.m_Shell.size(); ++i, ++it) delete it->m_psave;

	it = m_dic.m_Beam.begin();
	for (i=0; i<(int) m_dic.m_Beam.size(); ++i, ++it) delete it->m_psave;
}

//-----------------------------------------------------------------------------
bool FEBioPlotFile::Open(FEM &fem, const char *szfile)
{
	int i, j, N;

	// get the mesh
	FEMesh& m = fem.m_mesh;

	int nmode = fem.m_pStep->m_nModule;
	int ntype = fem.m_pStep->m_nanalysis;

	// setup the dictionary
	m_dic.AddNodalVariable(new FEPlotNodeDisplacement,  VEC3F, "Displacement");
	m_dic.AddSolidVariable(new FEPlotElementStress   , MAT3FS, "Stress");

	// store dynamic analysis data
	if ((ntype == FE_DYNAMIC) || (nmode == FE_POROELASTIC)) m_dic.AddNodalVariable(new FEPlotNodeVelocity    , VEC3F, "Velocity");
	if (ntype == FE_DYNAMIC) m_dic.AddNodalVariable(new FEPlotNodeAcceleration, VEC3F, "Acceleration");

	// store contact data
	if (fem.m_CI.size() > 0)
	{
		m_dic.AddNodalVariable(new FEPlotContactGap     , SCALAR, "Contact gap");
		m_dic.AddNodalVariable(new FEPlotContactTraction,  VEC3F, "Contact traction");
	}

	// store poro data
	if (nmode == FE_POROELASTIC)
	{
		m_dic.AddNodalVariable(new FEPlotFluidPressure, SCALAR, "Fluid Pressure");
		m_dic.AddSolidVariable(new FEPlotFluidFlux    ,  VEC3F, "Fluid Flux");
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
		m_dic.AddSolidVariable(new FEPlotFiberVector, VEC3F, "Fiber vector");
	}

	// setup the header
	m_hdr.nsize = sizeof(HEADER);
	m_hdr.nversion = 1;
	m_hdr.nnodes = m.Nodes();
	m_hdr.n3d    = m.SolidElements();
	m_hdr.n2d    = m.ShellElements();
	m_hdr.n1d    = fem.m_DE.size();
	m_hdr.nmat   = fem.Materials();

	m_hdr.nglv = m_dic.m_Glob.size();
	m_hdr.nnv  = m_dic.m_Node.size();
	m_hdr.nv3d = m_dic.m_Elem.size();
	m_hdr.nv2d = m_dic.m_Shell.size();
	m_hdr.nv1d = m_dic.m_Beam.size();

	// open the archive
	if ((m_fp = fopen(szfile, "wb")) == 0) return false;

	// write the tag
	unsigned int tag = FEBIO_TAG;
	fwrite(&tag, sizeof(tag), 1, m_fp);

	// --- save the header file ---
	fwrite(&m_hdr, sizeof(HEADER), 1, m_fp);

	// --- save the dictionary ---
	m_dic.Save(m_fp);

	// --- save the geometry ---
	
	// write the material coordinates
	float xf[3];
	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);

		xf[0] = (float) node.m_r0.x;
		xf[1] = (float) node.m_r0.y;
		xf[2] = (float) node.m_r0.z;

		fwrite(xf, sizeof(float), 3, m_fp);
	}

	// write the connectivity and material number
	// Note that we increment all numbers by 1 since
	// the plot database expects 1-based arrays
	int n[9], nd;

	int nid = 1;

	// write solid element data
	// note that we reindex all elements so that the ID
	// corresponds to the nr in the plot file
	for (nd = 0; nd<m.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&m.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				el.m_nID = nid++;

				// store material number
				n[0] = el.GetMatID()+1;

				N = el.Nodes();
				switch (el.Type())
				{
				case FE_HEX:
				case FE_RIHEX:
				case FE_UDGHEX:
					for (j=0; j<N; ++j) n[j+1] = el.m_node[j]+1;
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

				fwrite(n, sizeof(int), 9, m_fp);
			}
		}
	}

	// write shell element data
	for (nd = 0; nd<m.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);

				el.m_nID = nid++;

				// save material ID
				n[0] = el.GetMatID()+1;

				N = el.Nodes();
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

				fwrite(n, sizeof(int), 5, m_fp);
			}
		}
	}

	// write truss element data
	for (nd = 0; nd < m.Domains(); ++nd)
	{
		FETrussDomain* ptd = dynamic_cast<FETrussDomain*>(&m.Domain(nd));
		if (ptd)
		{
			for (i=0; i<ptd->Elements(); ++i)
			{
				FETrussElement& el = ptd->Element(i);
				el.m_nID = nid++;
				n[0] = el.GetMatID()+1;
				n[1] = el.m_node[0]+1;
				n[2] = el.m_node[1]+1;

				fwrite(n, sizeof(int), 3, m_fp);
			}
		}
	}

	return true;
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
	if (m_fp == 0) return false;

	// store the fem pointer
	m_pfem = &fem;

	// get the mesh
	FEMesh& mesh = fem.m_mesh;

	// save the time stamp
	float time = (float) fem.m_ftime;
	fwrite(&time, sizeof(float), 1, m_fp);

	// save the global variables
	if (m_dic.m_Glob.size() > 0)
	{
	}

	// save the nodal variables
	if (m_dic.m_Node.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator it = m_dic.m_Node.begin();
		for (int i=0; i<(int) m_dic.m_Node.size(); ++i, ++it)
			if (it->m_psave) (it->m_psave)->Save(fem, m_fp);
	}

	// save the solid variables
	if (m_dic.m_Elem.size() > 0)
	{
		list<DICTIONARY_ITEM>::iterator it = m_dic.m_Elem.begin();
		for (int i=0; i<(int) m_dic.m_Elem.size(); ++i, ++it)
			if (it->m_psave) (it->m_psave)->Save(fem, m_fp);
	}

	// save the shell variables
	if (m_dic.m_Shell.size() > 0)
	{
	}

	// save the beam variables
	if (m_dic.m_Beam.size() > 0)
	{
	}

	return true;
}
