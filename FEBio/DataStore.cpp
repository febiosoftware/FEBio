// DataStore.cpp: implementation of the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DataStore.h"
#include "fem.h"
#include "FERigid.h"
#include "FESolidSolver.h"
#include "log.h"

//////////////////////////////////////////////////////////////////////
// DataStore
//////////////////////////////////////////////////////////////////////

DataStore::DataStore()
{
}

DataStore::~DataStore()
{
	for (size_t i=0; i<m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

//-----------------------------------------------------------------------------

void DataStore::Write()
{
	for (size_t i=0; i<m_data.size(); ++i)
	{
		DataRecord& DR = *m_data[i];
		DR.Write();
	}
}

//-----------------------------------------------------------------------------

void DataStore::AddRecord(DataRecord* prec)
{
	static int nid = 1;
	prec->m_nid = nid++;
	m_data.push_back(prec);
}

//////////////////////////////////////////////////////////////////////
// DataRecord
//////////////////////////////////////////////////////////////////////

DataRecord::DataRecord(FEM* pfem, const char* szfile)
{
	m_pfem = pfem;
	m_fp = 0;
	m_szfile[0] = 0;
	m_bcomm = true;

	if (szfile)
	{
		strcpy(m_szfile, szfile);
		m_fp = fopen(szfile, "wt");
		fprintf(m_fp, "*Title:%s\n", pfem->GetTitle());
	}

	m_nid = 0;
	m_szdata[0] = 0;
	m_szname[0] = 0;
	m_szfile[0] = 0;
	strcpy(m_szdelim," ");
}

DataRecord::~DataRecord()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}
}

bool DataRecord::Write()
{
	FEM& fem = *m_pfem;

	int nstep = fem.m_pStep->m_ntimesteps;
	double ftime = fem.m_ftime;
	char* sz, *ch;
	if (strlen(m_szname) == 0) sz = m_szdata; else sz = m_szname;
	double val;

	FILE* fplog = (FILE*) clog;

	// make a note in the log file
	fprintf(fplog, "\nData Record #%d\n", m_nid);
	fprintf(fplog, "===========================================================================\n");
	fprintf(fplog, "Step = %d\n", nstep);
	fprintf(fplog, "Time = %.9lg\n", ftime);
	fprintf(fplog, "Data = %s\n", sz);

	// see if we are saving the data to the logfile or to a 
	// seperate data file
	FILE* fp = m_fp;
	if (fp == 0)
	{
		// we store the data in the logfile
		fp = fplog;
	}
	else if (m_bcomm)
	{
		// we save the data in a seperate file
		fprintf(fplog, "File = %s\n", m_szfile);

		// make a note in the data file
		fprintf(fp,"*Step  = %d\n", nstep);
		fprintf(fp,"*Time  = %.9lg\n", ftime);
		fprintf(fp,"*Data  = %s\n", sz);
	}

	// save the data
	for (size_t i=0; i<m_item.size(); ++i)
	{
		fprintf(fp, "%d%s", m_item[i], m_szdelim);
		sz = m_szdata;
		do
		{
			ch = strchr(sz, ';');
			if (ch) *ch = 0;
			val = Evaluate(m_item[i], sz);
			fprintf(fp, "%lg", val);
			if (ch) 
			{
				*ch = ';';
				sz = ch+1;
				fprintf(fp, "%s", m_szdelim);
			}
			else fprintf(fp, "\n");
		}
		while (ch != 0);
	}

	return true;
}

//-----------------------------------------------------------------------------

void DataRecord::SetItemList(const char* szlist)
{
	int i, n = 0, n0, n1, nn;
	char* ch;
	char* sz = (char*) szlist;
	int nread;
	do
	{
		ch = strchr(sz, ',');
		if (ch) *ch = 0;
		nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		case 3:
			break;
		default:
			n0 = 0;
			n1 = -1;
			nn = 1;
		}

		for (i=n0; i<=n1; i += nn) ++n;

		if (ch) *ch = ',';
		sz = ch+1;
	}
	while (ch != 0);

	if (n != 0)
	{
		m_item.resize(n);

		sz = (char*) szlist;
		n = 0;
		do
		{
			ch = strchr(sz, ',');
			if (ch) *ch = 0;
			nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
			switch (nread)
			{
			case 1:
				n1 = n0;
				nn = 1;
				break;
			case 2:
				nn = 1;
			}

			for (i=n0; i<=n1; i += nn) m_item[n++] = i;
			assert(n <= (int) m_item.size());

			if (ch) *ch = ',';
			sz = ch+1;
		}
		while (ch != 0);
	}
	else SelectAllItems();
}

//-----------------------------------------------------------------------------

double NodeDataRecord::Evaluate(int item, const char* szexpr)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*fem.m_pStep->m_psolver);
	vector<double>& Fr = solver.m_Fr;
	int nnode = item - 1;
	double val = 0;
	int ierr;
	if ((nnode >= 0) && (nnode < mesh.Nodes()))
	{
		FENode& node = mesh.Node(nnode);
		int* id = node.m_ID;
		m_calc.SetVariable("x", node.m_rt.x);
		m_calc.SetVariable("y", node.m_rt.y);
		m_calc.SetVariable("z", node.m_rt.z);
		m_calc.SetVariable("ux", node.m_rt.x - node.m_r0.x);
		m_calc.SetVariable("uy", node.m_rt.y - node.m_r0.y);
		m_calc.SetVariable("uz", node.m_rt.z - node.m_r0.z);
		m_calc.SetVariable("Rx", (-id[0] - 2 >= 0 ? Fr[-id[0]-2] : 0));
		m_calc.SetVariable("Ry", (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0));
		m_calc.SetVariable("Rz", (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0));
		if (fem.m_pStep->m_nModule == FE_POROELASTIC)
		{
			m_calc.SetVariable("p", node.m_pt);
			m_calc.SetVariable("vx", node.m_vt.x);
			m_calc.SetVariable("vy", node.m_vt.y);
			m_calc.SetVariable("vz", node.m_vt.z);
		}
		else if (fem.m_pStep->m_nModule == FE_POROSOLUTE)
		{
			m_calc.SetVariable("p", node.m_pt);
			m_calc.SetVariable("c", node.m_ct);
			m_calc.SetVariable("vx", node.m_vt.x);
			m_calc.SetVariable("vy", node.m_vt.y);
			m_calc.SetVariable("vz", node.m_vt.z);
		}
		val = m_calc.eval(szexpr, ierr);
	}
	return val;
}

void NodeDataRecord::SelectAllItems()
{
	int n = m_pfem->m_mesh.Nodes();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

void NodeDataRecord::SetItemList(FENodeSet* pns)
{
	int n = pns->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pns)[i];
}

//-----------------------------------------------------------------------------
// TODO: this function will become very slow for a while
//       I need to fix this as soon as possible
double ElementDataRecord::Evaluate(int item, const char* szexpr)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	double val = 0;
	int ierr;

	// find the element
	// note that currently item is an index into the element array. When shells
	// are combined with solid elements, this may not result in the correct element as in the input file.
	// However, in any case it should correspond to the correct element in the plot file
	// TODO: THIS IS UNACCEPTABLY SLOW !!!
	FEElement* pe = mesh.FindElementFromID(item);
	mat3ds E;
	if (dynamic_cast<FESolidElement*>(pe)) 
	{
		// this is a solid element
		FESolidElement& el = dynamic_cast<FESolidElement&>(*pe);
		// TODO: do I need to unpack this element?
//		bd.UnpackElement(el);

		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();

			E = pt.Strain();
			m_calc.SetVariable("sx", pt.s.xx());
			m_calc.SetVariable("sy", pt.s.yy());
			m_calc.SetVariable("sz", pt.s.zz());
			m_calc.SetVariable("sxy", pt.s.xy());
			m_calc.SetVariable("syz", pt.s.yz());
			m_calc.SetVariable("sxz", pt.s.xz());
			m_calc.SetVariable("Ex", E.xx());
			m_calc.SetVariable("Ey", E.yy());
			m_calc.SetVariable("Ez", E.zz());
			m_calc.SetVariable("Exy", E.xy());
			m_calc.SetVariable("Eyz", E.yz());
			m_calc.SetVariable("Exz", E.xz());
			val += m_calc.eval(szexpr, ierr);
			if (fem.m_pStep->m_nModule == FE_POROELASTIC)
			{
				FEPoroElasticMaterialPoint* ppt = el.m_State[i]->ExtractData<FEPoroElasticMaterialPoint>();
				if (ppt)
				{
					m_calc.SetVariable("wx", ppt->m_w.x);
					m_calc.SetVariable("wy", ppt->m_w.y);
					m_calc.SetVariable("wz", ppt->m_w.z);
				}
			}
			else if (fem.m_pStep->m_nModule == FE_POROSOLUTE)
			{
				FESolutePoroElasticMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutePoroElasticMaterialPoint>();
				if (ppt)
				{
					m_calc.SetVariable("wx", ppt->m_w.x);
					m_calc.SetVariable("wy", ppt->m_w.y);
					m_calc.SetVariable("wz", ppt->m_w.z);
					m_calc.SetVariable("jx", ppt->m_j.x);
					m_calc.SetVariable("jy", ppt->m_j.y);
					m_calc.SetVariable("jz", ppt->m_j.z);
				}
			}
		}
		val /= nint;
	}
	else if (dynamic_cast<FEShellElement*>(pe))
	{
		// this is a shell element
		FEShellElement& el = dynamic_cast<FEShellElement&>(*pe);
		// TODO: do I need to unpack this element?
//		sd.UnpackElement(el);

		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();

			E = pt.Strain();
			m_calc.SetVariable("sx", pt.s.xx());
			m_calc.SetVariable("sy", pt.s.yy());
			m_calc.SetVariable("sz", pt.s.zz());
			m_calc.SetVariable("sxy", pt.s.xy());
			m_calc.SetVariable("syz", pt.s.yz());
			m_calc.SetVariable("sxz", pt.s.xz());
			m_calc.SetVariable("Ex", E.xx());
			m_calc.SetVariable("Ey", E.yy());
			m_calc.SetVariable("Ez", E.zz());
			m_calc.SetVariable("Exy", E.xy());
			m_calc.SetVariable("Eyz", E.yz());
			m_calc.SetVariable("Exz", E.xz());
			val += m_calc.eval(szexpr, ierr);
		}
		val /= nint;
	}

	return val;
}

void ElementDataRecord::SelectAllItems()
{
	FEMesh& m = m_pfem->m_mesh;
	int n = m.Elements();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//-----------------------------------------------------------------------------

double RigidBodyDataRecord::Evaluate(int item, const char* szexpr)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;
	int nrb = item - 1;
	double val = 0;
	int ierr;
	if ((nrb >= 0) && (nrb < fem.Materials()))
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(nrb));
		if (pm)
		{
			// find the rigid body that has this material
			for (int i=0; i<fem.m_nrb; ++i)
			{
				FERigidBody& RB = fem.m_RB[i];
				if (RB.m_mat == nrb)
				{
					m_calc.SetVariable("x", RB.m_rt.x);
					m_calc.SetVariable("y", RB.m_rt.y);
					m_calc.SetVariable("z", RB.m_rt.z);
					m_calc.SetVariable("qx", RB.m_qt.x);
					m_calc.SetVariable("qy", RB.m_qt.y);
					m_calc.SetVariable("qz", RB.m_qt.z);
					m_calc.SetVariable("qw", RB.m_qt.w);
					m_calc.SetVariable("Fx", RB.m_Fr.x);
					m_calc.SetVariable("Fy", RB.m_Fr.y);
					m_calc.SetVariable("Fz", RB.m_Fr.z);
					m_calc.SetVariable("Mx", RB.m_Mr.x);
					m_calc.SetVariable("My", RB.m_Mr.y);
					m_calc.SetVariable("Mz", RB.m_Mr.z);
					val = m_calc.eval(szexpr, ierr);
					break;
				}
			}
		}
	}

	return val;
}

void RigidBodyDataRecord::SelectAllItems()
{
	int n = 0, i;
	for (i=0; i<m_pfem->Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(i));
		if (pm) ++n;
	}

	if (n > 0)
	{
		m_item.resize(n);
		n = 0;
		for (i=0; i<m_pfem->Materials(); ++i)
		{
			FERigidMaterial* pm  = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(i));
			if (pm)
			{
				m_item[n++] = i+1;
			}
		}
	}
}
