#include "stdafx.h"
#include "ElementDataRecord.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"

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
		else if (strcmp(sz, "E1" ) == 0) m_data.push_back(E1);
		else if (strcmp(sz, "E2" ) == 0) m_data.push_back(E2);
		else if (strcmp(sz, "E3" ) == 0) m_data.push_back(E3);
		else if (strcmp(sz, "sx" ) == 0) m_data.push_back(SX );
		else if (strcmp(sz, "sy" ) == 0) m_data.push_back(SY );
		else if (strcmp(sz, "sz" ) == 0) m_data.push_back(SZ );
		else if (strcmp(sz, "sxy") == 0) m_data.push_back(SXY);
		else if (strcmp(sz, "syz") == 0) m_data.push_back(SYZ);
		else if (strcmp(sz, "sxz") == 0) m_data.push_back(SXZ);
		else if (strcmp(sz, "s1" ) == 0) m_data.push_back(S1 );
		else if (strcmp(sz, "s2" ) == 0) m_data.push_back(S2 );
		else if (strcmp(sz, "s3" ) == 0) m_data.push_back(S3 );
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
		else if (strcmp(sz, "crc") == 0) m_data.push_back(CRC);
		else if (strcmp(sz, "c1"  ) == 0) m_data.push_back(C1 );
		else if (strcmp(sz, "c2"  ) == 0) m_data.push_back(C2 );
		else if (strcmp(sz, "c3"  ) == 0) m_data.push_back(C3 );
		else if (strcmp(sz, "c4"  ) == 0) m_data.push_back(C4 );
		else if (strcmp(sz, "c5"  ) == 0) m_data.push_back(C5 );
		else if (strcmp(sz, "c6"  ) == 0) m_data.push_back(C6 );
		else if (strcmp(sz, "j1x" ) == 0) m_data.push_back(J1X);
		else if (strcmp(sz, "j1y" ) == 0) m_data.push_back(J1Y);
		else if (strcmp(sz, "j1z" ) == 0) m_data.push_back(J1Z);
		else if (strcmp(sz, "j2x" ) == 0) m_data.push_back(J2X);
		else if (strcmp(sz, "j2y" ) == 0) m_data.push_back(J2Y);
		else if (strcmp(sz, "j2z" ) == 0) m_data.push_back(J2Z);
		else if (strcmp(sz, "j3x" ) == 0) m_data.push_back(J3X);
		else if (strcmp(sz, "j3y" ) == 0) m_data.push_back(J3Y);
		else if (strcmp(sz, "j3z" ) == 0) m_data.push_back(J3Z);
		else if (strcmp(sz, "j4x" ) == 0) m_data.push_back(J4X);
		else if (strcmp(sz, "j4y" ) == 0) m_data.push_back(J4Y);
		else if (strcmp(sz, "j4z" ) == 0) m_data.push_back(J4Z);
		else if (strcmp(sz, "j5x" ) == 0) m_data.push_back(J5X);
		else if (strcmp(sz, "j5y" ) == 0) m_data.push_back(J5Y);
		else if (strcmp(sz, "j5z" ) == 0) m_data.push_back(J5Z);
		else if (strcmp(sz, "j6x" ) == 0) m_data.push_back(J6X);
		else if (strcmp(sz, "j6y" ) == 0) m_data.push_back(J6Y);
		else if (strcmp(sz, "j6z" ) == 0) m_data.push_back(J6Z);
		else if (strcmp(sz, "psi" ) == 0) m_data.push_back(PSI);
		else if (strcmp(sz, "Iex" ) == 0) m_data.push_back(IEX);
		else if (strcmp(sz, "Iey" ) == 0) m_data.push_back(IEY);
		else if (strcmp(sz, "Iez" ) == 0) m_data.push_back(IEZ);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double ElementDataRecord::Evaluate(int item, int ndata)
{
	FEMesh& mesh = m_pfem->GetMesh();

	// make sure we have an ELT
	if (m_ELT.empty()) BuildELT();

	// find the element
	assert((item >= 1) && (item <= mesh.Elements()));
	ELEMREF e = m_ELT[item-1];
	assert((e.ndom != -1) && (e.nid != -1));
	FEElement* pe = &mesh.Domain(e.ndom).ElementRef(e.nid);

	// see if we need to calculate strain
	bool bE = false;
	if ((ndata>=EX)&&(ndata<=E3)) bE = true;

	// see if ndata requires eigen values
	bool blE = false; if ((ndata>=E1)&&(ndata<=E3)) blE = true;
	bool bls = false; if ((ndata>=S1)&&(ndata<=S3)) bls = true;

	// calculate the return val
	double val = 0;
	mat3ds E;
	double lE[3], ls[3];
	if (dynamic_cast<FESolidElement*>(pe)) 
	{
		// this is a solid element
		FESolidElement& el = dynamic_cast<FESolidElement&>(*pe);

		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
			if (bE) E = pt.Strain();
			if (blE) E.exact_eigen(lE);
			if (bls) pt.s.exact_eigen(ls);

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
			case E1: val += lE[0]; break;
			case E2: val += lE[1]; break;
			case E3: val += lE[2]; break;
			case SX: val += pt.s.xx(); break;
			case SY: val += pt.s.yy(); break;
			case SZ: val += pt.s.zz(); break;
			case SXY: val += pt.s.xy(); break;
			case SYZ: val += pt.s.yz(); break;
			case SXZ: val += pt.s.xz(); break;
			case S1: val += ls[0]; break;
			case S2: val += ls[1]; break;
			case S3: val += ls[2]; break;
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

			FEBiphasicMaterialPoint* ppt = el.m_State[i]->ExtractData<FEBiphasicMaterialPoint>();
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

			FESoluteMaterialPoint* spt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
			if (spt)
			{
				switch (ndata)
				{
					case C: val += spt->m_ca; break;
					case JX: val += spt->m_j.x; break;
					case JY: val += spt->m_j.y; break;
					case JZ: val += spt->m_j.z; break;
				}
			}

			FESaltMaterialPoint* stt = el.m_State[i]->ExtractData<FESaltMaterialPoint>();
			if (stt)
			{
				switch (ndata)
				{
					case C1: val += stt->m_ca[0]; break;
					case J1X: val += stt->m_j[0].x; break;
					case J1Y: val += stt->m_j[0].y; break;
					case J1Z: val += stt->m_j[0].z; break;
					case C2: val += stt->m_ca[1]; break;
					case J2X: val += stt->m_j[1].x; break;
					case J2Y: val += stt->m_j[1].y; break;
					case J2Z: val += stt->m_j[1].z; break;
					case PSI: val += stt->m_psi; break;
					case IEX: val += stt->m_Ie.x; break;
					case IEY: val += stt->m_Ie.y; break;
					case IEZ: val += stt->m_Ie.z; break;
				}
			}
			FESolutesMaterialPoint* sst = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
			if (sst)
			{
				switch (ndata)
				{
					case C1: val += sst->m_ca[0]; break;
					case J1X: val += sst->m_j[0].x; break;
					case J1Y: val += sst->m_j[0].y; break;
					case J1Z: val += sst->m_j[0].z; break;
					case C2: val += sst->m_ca[1]; break;
					case J2X: val += sst->m_j[1].x; break;
					case J2Y: val += sst->m_j[1].y; break;
					case J2Z: val += sst->m_j[1].z; break;
					case C3: val += sst->m_ca[2]; break;
					case J3X: val += sst->m_j[2].x; break;
					case J3Y: val += sst->m_j[2].y; break;
					case J3Z: val += sst->m_j[2].z; break;
					case C4: val += sst->m_ca[3]; break;
					case J4X: val += sst->m_j[3].x; break;
					case J4Y: val += sst->m_j[3].y; break;
					case J4Z: val += sst->m_j[3].z; break;
					case C5: val += sst->m_ca[4]; break;
					case J5X: val += sst->m_j[4].x; break;
					case J5Y: val += sst->m_j[4].y; break;
					case J5Z: val += sst->m_j[4].z; break;
					case C6: val += sst->m_ca[5]; break;
					case J6X: val += sst->m_j[5].x; break;
					case J6Y: val += sst->m_j[5].y; break;
					case J6Z: val += sst->m_j[5].z; break;
					case PSI: val += sst->m_psi; break;
					case IEX: val += sst->m_Ie.x; break;
					case IEY: val += sst->m_Ie.y; break;
					case IEZ: val += sst->m_Ie.z; break;
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
	FEMesh& m = m_pfem->GetMesh();
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
	FEMesh& m = m_pfem->GetMesh();
	int n = m.Elements();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}
