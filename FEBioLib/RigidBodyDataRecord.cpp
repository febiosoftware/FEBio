#include "stdafx.h"
#include "RigidBodyDataRecord.h"
#include <FEBioMech/FERigid.h>
#include <FEBioMech/FERigidBody.h>

//-----------------------------------------------------------------------------
double FELogRigidBodyPosX::value(FERigidBody& rb) { return rb.m_rt.x; }
double FELogRigidBodyPosY::value(FERigidBody& rb) { return rb.m_rt.y; }
double FELogRigidBodyPosZ::value(FERigidBody& rb) { return rb.m_rt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyQuatX::value(FERigidBody& rb) { return rb.m_qt.x; }
double FELogRigidBodyQuatY::value(FERigidBody& rb) { return rb.m_qt.y; }
double FELogRigidBodyQuatZ::value(FERigidBody& rb) { return rb.m_qt.z; }
double FELogRigidBodyQuatW::value(FERigidBody& rb) { return rb.m_qt.w; }

//-----------------------------------------------------------------------------
double FELogRigidBodyForceX::value(FERigidBody& rb) { return rb.m_Fr.x; }
double FELogRigidBodyForceY::value(FERigidBody& rb) { return rb.m_Fr.y; }
double FELogRigidBodyForceZ::value(FERigidBody& rb) { return rb.m_Fr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyTorqueX::value(FERigidBody& rb) { return rb.m_Mr.x; }
double FELogRigidBodyTorqueY::value(FERigidBody& rb) { return rb.m_Mr.y; }
double FELogRigidBodyTorqueZ::value(FERigidBody& rb) { return rb.m_Mr.z; }

//-----------------------------------------------------------------------------
void RigidBodyDataRecord::Parse(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x" ) == 0) m_Data.push_back(new FELogRigidBodyPosX(m_pfem));
		else if (strcmp(sz, "y" ) == 0) m_Data.push_back(new FELogRigidBodyPosY(m_pfem));
		else if (strcmp(sz, "z" ) == 0) m_Data.push_back(new FELogRigidBodyPosZ(m_pfem));
		else if (strcmp(sz, "qx") == 0) m_Data.push_back(new FELogRigidBodyQuatX(m_pfem));
		else if (strcmp(sz, "qy") == 0) m_Data.push_back(new FELogRigidBodyQuatY(m_pfem));
		else if (strcmp(sz, "qz") == 0) m_Data.push_back(new FELogRigidBodyQuatZ(m_pfem));
		else if (strcmp(sz, "qw") == 0) m_Data.push_back(new FELogRigidBodyQuatW(m_pfem));
		else if (strcmp(sz, "Fx") == 0) m_Data.push_back(new FELogRigidBodyForceX(m_pfem));
		else if (strcmp(sz, "Fy") == 0) m_Data.push_back(new FELogRigidBodyForceY(m_pfem));
		else if (strcmp(sz, "Fz") == 0) m_Data.push_back(new FELogRigidBodyForceZ(m_pfem));
		else if (strcmp(sz, "Mx") == 0) m_Data.push_back(new FELogRigidBodyTorqueX(m_pfem));
		else if (strcmp(sz, "My") == 0) m_Data.push_back(new FELogRigidBodyTorqueY(m_pfem));
		else if (strcmp(sz, "Mz") == 0) m_Data.push_back(new FELogRigidBodyTorqueZ(m_pfem));
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}
//-----------------------------------------------------------------------------
double RigidBodyDataRecord::Evaluate(int item, int ndata)
{
	FEMesh& mesh = m_pfem->GetMesh();
	int nrb = item - 1;
	if ((nrb < 0) || (nrb >= m_pfem->Materials())) return 0;

	double val = 0;
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(nrb));
	assert(pm);
	if (pm == 0) return 0;

	// find the rigid body that has this material
	int NRB = m_pfem->Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& RB = dynamic_cast<FERigidBody&>(*m_pfem->Object(i));
		if (RB.m_mat == nrb) return m_Data[ndata]->value(RB);
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
