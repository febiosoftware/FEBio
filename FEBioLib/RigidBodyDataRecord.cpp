#include "stdafx.h"
#include "RigidBodyDataRecord.h"
#include "FERigid.h"
#include "FERigidBody.h"

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
	FEMesh& mesh = m_pfem->m_mesh;
	int nrb = item - 1;
	if ((nrb < 0) || (nrb >= m_pfem->Materials())) return 0;

	double val = 0;
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(nrb));
	assert(pm);
	if (pm == 0) return 0;

	// find the rigid body that has this material
	int NRB = m_pfem->m_Obj.size();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_pfem->m_Obj[i]);
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
