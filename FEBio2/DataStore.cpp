// DataStore.cpp: implementation of the DataStore class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DataStore.h"
#include "fem.h"
#include "FEBioLib/FERigid.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
UnknownDataField::UnknownDataField(const char* sz)
{
	m_szdata[0] = 0;
	int l = strlen(sz);
	if (l > 63) l = 63;
	if (l>0) strncpy(m_szdata, sz, l);
}

//////////////////////////////////////////////////////////////////////
// DataStore
//////////////////////////////////////////////////////////////////////

DataStore::DataStore()
{
}

DataStore::~DataStore()
{
}

void DataStore::Clear()
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

//-----------------------------------------------------------------------------

void DataStore::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << (int) m_data.size();
		for (int i=0; i<(int) m_data.size(); ++i)
		{
			DataRecord* pd = m_data[i];

			int ntype = -1;
			if (dynamic_cast<NodeDataRecord*>(pd)) ntype = FE_DATA_NODE;
			if (dynamic_cast<ElementDataRecord*>(pd)) ntype = FE_DATA_ELEM;
			if (dynamic_cast<RigidBodyDataRecord*>(pd)) ntype = FE_DATA_RB;
			assert(ntype != -1);
			ar << ntype;
			pd->Serialize(ar);
		}
	}
	else
	{
		FEM* pfem = dynamic_cast<FEM*>(ar.GetFEM());

		int ndr;
		Clear();
		ar >> ndr;
		for (int i=0; i<ndr; ++i)
		{
			int ntype;
			ar >> ntype;

			DataRecord* pd = 0;
			switch(ntype)
			{
			case FE_DATA_NODE: pd = new NodeDataRecord(pfem, 0); break;
			case FE_DATA_ELEM: pd = new ElementDataRecord(pfem, 0); break;
			case FE_DATA_RB  : pd = new RigidBodyDataRecord(pfem, 0); break;
			}
			assert(pd);
			pd->Serialize(ar);
			AddRecord(pd);
		}
	}
}

//////////////////////////////////////////////////////////////////////
// DataRecord
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
DataRecord::DataRecord(FEM* pfem, const char* szfile)
{
	m_pfem = pfem;
	m_nid = 0;
	m_szname[0] = 0;

	strcpy(m_szdelim, " ");
	
	m_bcomm = true;

	m_fp = 0;
	m_szfile[0] = 0;

	if (szfile)
	{
		strcpy(m_szfile, szfile);
		m_fp = fopen(szfile, "wt");
		fprintf(m_fp, "*Title:%s\n", pfem->GetTitle());
	}
}

//-----------------------------------------------------------------------------
DataRecord::~DataRecord()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}
}

//-----------------------------------------------------------------------------
void DataRecord::SetName(const char* sz)
{
	strcpy(m_szname, sz);
}

//-----------------------------------------------------------------------------
void DataRecord::SetDelim(const char* sz)
{
	strcpy(m_szdelim, sz);
}

//-----------------------------------------------------------------------------
bool DataRecord::Write()
{
	FEM& fem = *m_pfem;

	int nstep = fem.m_pStep->m_ntimesteps;
	double ftime = fem.m_ftime;
	double val;

	FILE* fplog = (FILE*) clog;

	// make a note in the log file
	fprintf(fplog, "\nData Record #%d\n", m_nid);
	fprintf(fplog, "===========================================================================\n");
	fprintf(fplog, "Step = %d\n", nstep);
	fprintf(fplog, "Time = %.9lg\n", ftime);
	fprintf(fplog, "Data = %s\n", m_szname);

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
		fprintf(fp,"*Data  = %s\n", m_szname);
	}

	// save the data
	for (size_t i=0; i<m_item.size(); ++i)
	{
		fprintf(fp, "%d%s", m_item[i], m_szdelim);
		int nd = (int) m_data.size();
		for (int j=0; j<nd; ++j)
		{
			val = Evaluate(m_item[i], m_data[j]);
			fprintf(fp, "%lg", val);
			if (j!=nd-1) fprintf(fp, "%s", m_szdelim);
			else fprintf(fp, "\n");
		}
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

void DataRecord::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_nid;
		ar << m_szname;
		ar << m_szdelim;
		ar << m_szfile;
		ar << m_bcomm;
		ar << m_item;
	}
	else
	{
		ar >> m_nid;
		ar >> m_szname;
		ar >> m_szdelim;
		ar >> m_szfile;
		ar >> m_bcomm;
		ar >> m_item;

		if (m_fp) fclose(m_fp);
		m_fp = 0;
		if (m_szfile[0] != 0)
		{
			// reopen data file for appending
			m_fp = fopen(m_szfile, "a+");
		}
	}
}

//=============================================================================
// NodeDataRecord
//-----------------------------------------------------------------------------
void NodeDataRecord::Parse(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x" ) == 0) m_data.push_back(X );
		else if (strcmp(sz, "y" ) == 0) m_data.push_back(Y );
		else if (strcmp(sz, "z" ) == 0) m_data.push_back(Z );
		else if (strcmp(sz, "ux") == 0) m_data.push_back(UX);
		else if (strcmp(sz, "uy") == 0) m_data.push_back(UY);
		else if (strcmp(sz, "uz") == 0) m_data.push_back(UZ);
		else if (strcmp(sz, "vx") == 0) m_data.push_back(VX);
		else if (strcmp(sz, "vy") == 0) m_data.push_back(VY);
		else if (strcmp(sz, "vz") == 0) m_data.push_back(VZ);
		else if (strcmp(sz, "Rx") == 0) m_data.push_back(RX);
		else if (strcmp(sz, "Ry") == 0) m_data.push_back(RY);
		else if (strcmp(sz, "Rz") == 0) m_data.push_back(RZ);
		else if (strcmp(sz, "p" ) == 0) m_data.push_back(P );
		else if (strcmp(sz, "c" ) == 0) m_data.push_back(C );
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double NodeDataRecord::Evaluate(int item, int ndata)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*fem.m_pStep->m_psolver);
	vector<double>& Fr = solver.m_Fr;
	int nnode = item - 1;
	if ((nnode < 0) || (nnode >= mesh.Nodes())) return 0;
	FENode& node = mesh.Node(nnode);
	int* id = node.m_ID;

	double val = 0;
	switch (ndata)
	{
	case X : val = node.m_rt.x; break;
	case Y : val = node.m_rt.y; break;
	case Z : val = node.m_rt.z; break;
	case UX: val = node.m_rt.x - node.m_r0.x; break;
	case UY: val = node.m_rt.y - node.m_r0.y; break;
	case UZ: val = node.m_rt.z - node.m_r0.z; break;
	case VX: val = node.m_vt.x; break;
	case VY: val = node.m_vt.y; break;
	case VZ: val = node.m_vt.z; break;
	case RX: val = (-id[0] - 2 >= 0 ? Fr[-id[0]-2] : 0);  break;
	case RY: val = (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);  break;
	case RZ: val = (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);  break;
	case P : val = node.m_pt; break;
	case C : val = node.m_ct; break;
	}
	return val;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SelectAllItems()
{
	int n = m_pfem->m_mesh.Nodes();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SetItemList(FENodeSet* pns)
{
	int n = pns->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pns)[i];
}

//=============================================================================
// ElementDataRecord
//-----------------------------------------------------------------------------
void ElementDataRecord::Parse(const char *szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x"  ) == 0) m_data.push_back(X  );
		else if (strcmp(sz, "y"  ) == 0) m_data.push_back(Y  );
		else if (strcmp(sz, "z"  ) == 0) m_data.push_back(Z  );
		else if (strcmp(sz, "J"  ) == 0) m_data.push_back(J  );
		else if (strcmp(sz, "Ex" ) == 0) m_data.push_back(EX );
		else if (strcmp(sz, "Ey" ) == 0) m_data.push_back(EY );
		else if (strcmp(sz, "Ez" ) == 0) m_data.push_back(EZ );
		else if (strcmp(sz, "Exy") == 0) m_data.push_back(EXY);
		else if (strcmp(sz, "Eyz") == 0) m_data.push_back(EYZ);
		else if (strcmp(sz, "Exz") == 0) m_data.push_back(EXZ);
		else if (strcmp(sz, "sx" ) == 0) m_data.push_back(SX );
		else if (strcmp(sz, "sy" ) == 0) m_data.push_back(SY );
		else if (strcmp(sz, "sz" ) == 0) m_data.push_back(SZ );
		else if (strcmp(sz, "sxy") == 0) m_data.push_back(SXY);
		else if (strcmp(sz, "syz") == 0) m_data.push_back(SYZ);
		else if (strcmp(sz, "sxz") == 0) m_data.push_back(SXZ);
		else if (strcmp(sz, "Fx" ) == 0) m_data.push_back(FX );
		else if (strcmp(sz, "Fy" ) == 0) m_data.push_back(FY );
		else if (strcmp(sz, "Fz" ) == 0) m_data.push_back(FZ );
		else if (strcmp(sz, "Fyz") == 0) m_data.push_back(FYZ);
		else if (strcmp(sz, "Fzx") == 0) m_data.push_back(FZX);
		else if (strcmp(sz, "Fxy") == 0) m_data.push_back(FXY);
		else if (strcmp(sz, "Fyx") == 0) m_data.push_back(FYX);
		else if (strcmp(sz, "Fxz") == 0) m_data.push_back(FXZ);
		else if (strcmp(sz, "Fzy") == 0) m_data.push_back(FZY);
		else if (strcmp(sz, "p"  ) == 0) m_data.push_back(P  );
		else if (strcmp(sz, "wx" ) == 0) m_data.push_back(WX );
		else if (strcmp(sz, "wy" ) == 0) m_data.push_back(WY );
		else if (strcmp(sz, "wz" ) == 0) m_data.push_back(WZ );
		else if (strcmp(sz, "c"  ) == 0) m_data.push_back(C  );
		else if (strcmp(sz, "jx" ) == 0) m_data.push_back(JX );
		else if (strcmp(sz, "jy" ) == 0) m_data.push_back(JY );
		else if (strcmp(sz, "jz" ) == 0) m_data.push_back(JZ );
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double ElementDataRecord::Evaluate(int item, int ndata)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	// make sure we have an ELT
	if (m_ELT.empty()) BuildELT();

	// find the element
	assert((item >= 1) && (item <= mesh.Elements()));
	ELEMREF e = m_ELT[item-1];
	assert((e.ndom != -1) && (e.nid != -1));
	FEElement* pe = &mesh.Domain(e.ndom).ElementRef(e.nid);

	// calculate the return val
	double val = 0;
	mat3ds E;
	if (dynamic_cast<FESolidElement*>(pe)) 
	{
		// this is a solid element
		FESolidElement& el = dynamic_cast<FESolidElement&>(*pe);

		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
			E = pt.Strain();

			switch (ndata)
			{
			case X: val += pt.rt.x; break;
			case Y: val += pt.rt.y; break;
			case Z: val += pt.rt.z; break;
			case J: val += pt.J; break;
			case EX: val += E.xx(); break;
			case EY: val += E.yy(); break;
			case EZ: val += E.zz(); break;
			case EXY: val += E.xy(); break;
			case EYZ: val += E.yz(); break;
			case EXZ: val += E.xz(); break;
			case SX: val += pt.s.xx(); break;
			case SY: val += pt.s.yy(); break;
			case SZ: val += pt.s.zz(); break;
			case SXY: val += pt.s.xy(); break;
			case SYZ: val += pt.s.yz(); break;
			case SXZ: val += pt.s.xz(); break;
			case FX: val += pt.F(0,0); break;
			case FY: val += pt.F(1,1); break;
			case FZ: val += pt.F(2,2); break;
			case FYZ: val += pt.F(1,2); break;
			case FZX: val += pt.F(2,0); break;
			case FXY: val += pt.F(0,1); break;
			case FYX: val += pt.F(1,0); break;
			case FXZ: val += pt.F(0,2); break;
			case FZY: val += pt.F(2,1); break;
			}

			if (fem.m_pStep->m_nModule == FE_POROELASTIC)
			{
				FEPoroElasticMaterialPoint* ppt = el.m_State[i]->ExtractData<FEPoroElasticMaterialPoint>();
				if (ppt)
				{
					switch (ndata)
					{
					case P: val += ppt->m_pa; break;
					case WX: val += ppt->m_w.x; break;
					case WY: val += ppt->m_w.y; break;
					case WZ: val += ppt->m_w.z; break;
					}
				}
			}
			else if (fem.m_pStep->m_nModule == FE_POROSOLUTE)
			{
				FESolutePoroElasticMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutePoroElasticMaterialPoint>();
				if (ppt)
				{
					switch (ndata)
					{
					case P: val += ppt->m_pa; break;
					case WX: val += ppt->m_w.x; break;
					case WY: val += ppt->m_w.y; break;
					case WZ: val += ppt->m_w.z; break;
					case C: val += ppt->m_ca; break;
					case JX: val += ppt->m_j.x; break;
					case JY: val += ppt->m_j.y; break;
					case JZ: val += ppt->m_j.z; break;
					}
				}
			}
		}
		val /= nint;
	}
	else if (dynamic_cast<FEShellElement*>(pe))
	{
		// this is a shell element
		FEShellElement& el = dynamic_cast<FEShellElement&>(*pe);
		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();

			E = pt.Strain();
			switch (ndata)
			{
			case X: val += pt.rt.x; break;
			case Y: val += pt.rt.y; break;
			case Z: val += pt.rt.z; break;
			case J: val += pt.J; break;
			case EX: val += E.xx(); break;
			case EY: val += E.yy(); break;
			case EZ: val += E.zz(); break;
			case EXY: val += E.xy(); break;
			case EYZ: val += E.yz(); break;
			case EXZ: val += E.xz(); break;
			case SX: val += pt.s.xx(); break;
			case SY: val += pt.s.yy(); break;
			case SZ: val += pt.s.zz(); break;
			case SXY: val += pt.s.xy(); break;
			case SYZ: val += pt.s.yz(); break;
			case SXZ: val += pt.s.xz(); break;
			case FX: val += pt.F(0,0); break;
			case FY: val += pt.F(1,1); break;
			case FZ: val += pt.F(2,2); break;
			case FYZ: val += pt.F(1,2); break;
			case FZX: val += pt.F(2,0); break;
			case FXY: val += pt.F(0,1); break;
			case FYX: val += pt.F(1,0); break;
			case FXZ: val += pt.F(0,2); break;
			case FZY: val += pt.F(2,1); break;
			}
		}
		val /= nint;
	}
	return val;
}

//-----------------------------------------------------------------------------
void ElementDataRecord::BuildELT()
{
	int i, j;
	m_ELT.clear();
	FEMesh& m = m_pfem->m_mesh;
	int NE = m.Elements();
	m_ELT.resize(NE);
	for (i=0; i<NE; ++i) 
	{
		m_ELT[i].ndom = -1;
		m_ELT[i].nid  = -1;
	}

	for (i=0; i<m.Domains(); ++i)
	{
		FEDomain& d = m.Domain(i);
		int ne = d.Elements();
		for (j=0; j<ne; ++j)
		{
			FEElement& el = d.ElementRef(j);
			m_ELT[el.m_nID-1].ndom = i;
			m_ELT[el.m_nID-1].nid  = j;
		}
	}
}

//-----------------------------------------------------------------------------
void ElementDataRecord::SelectAllItems()
{
	FEMesh& m = m_pfem->m_mesh;
	int n = m.Elements();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//=============================================================================
// RigidBodyDataRecord
//-----------------------------------------------------------------------------
void RigidBodyDataRecord::Parse(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x" ) == 0) m_data.push_back(X );
		else if (strcmp(sz, "y" ) == 0) m_data.push_back(Y );
		else if (strcmp(sz, "z" ) == 0) m_data.push_back(Z );
		else if (strcmp(sz, "qx") == 0) m_data.push_back(QX);
		else if (strcmp(sz, "qy") == 0) m_data.push_back(QY);
		else if (strcmp(sz, "qz") == 0) m_data.push_back(QZ);
		else if (strcmp(sz, "qw") == 0) m_data.push_back(QW);
		else if (strcmp(sz, "Fx") == 0) m_data.push_back(FX);
		else if (strcmp(sz, "Fy") == 0) m_data.push_back(FY);
		else if (strcmp(sz, "Fz") == 0) m_data.push_back(FZ);
		else if (strcmp(sz, "Mx") == 0) m_data.push_back(MX);
		else if (strcmp(sz, "My") == 0) m_data.push_back(MY);
		else if (strcmp(sz, "Mz") == 0) m_data.push_back(MZ);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}
//-----------------------------------------------------------------------------
double RigidBodyDataRecord::Evaluate(int item, int ndata)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;
	int nrb = item - 1;
	if ((nrb < 0) || (nrb >= fem.Materials())) return 0;

	double val = 0;
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(nrb));
	assert(pm);
	if (pm == 0) return 0;

	// find the rigid body that has this material
	for (int i=0; i<fem.m_nrb; ++i)
	{
		FERigidBody& RB = fem.m_RB[i];
		if (RB.m_mat == nrb)
		{
			switch (ndata)
			{
			case X: val = RB.m_rt.x; break;
			case Y: val = RB.m_rt.y; break;
			case Z: val = RB.m_rt.z; break;
			case QX: val = RB.m_qt.x; break;
			case QY: val = RB.m_qt.y; break;
			case QZ: val = RB.m_qt.z; break;
			case QW: val = RB.m_qt.w; break;
			case FX: val = RB.m_Fr.x; break;
			case FY: val = RB.m_Fr.y; break;
			case FZ: val = RB.m_Fr.z; break;
			case MX: val = RB.m_Mr.x; break;
			case MY: val = RB.m_Mr.y; break;
			case MZ: val = RB.m_Mr.z; break;
			}
		}
	}

	return val;
}

//-----------------------------------------------------------------------------
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
