#include "stdafx.h"
#include "ObjectDataRecord.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include <FECore/FEMaterial.h>
#include "FEMechModel.h"
#include "FERigidMaterial.h"

//-----------------------------------------------------------------------------
void ObjectDataRecord::Parse(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		FELogObjectData* pdata = fecore_new<FELogObjectData>(FEOBJLOGDATA_ID, sz, m_pfem);
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double ObjectDataRecord::Evaluate(int item, int ndata)
{
	FEMechModel* fem = dynamic_cast<FEMechModel*>(m_pfem);

	FEMesh& mesh = fem->GetMesh();
	int nrb = item - 1;
	if ((nrb < 0) || (nrb >= m_pfem->Materials())) return 0;

	double val = 0;
	FEMaterial* pm = m_pfem->GetMaterial(nrb);
	assert(pm);
	if (pm == 0) return 0;

	// find the rigid body that has this material
	FERigidSystem& rs = *fem->GetRigidSystem();
	int NRB = rs.Objects();
	for (int i=0; i<NRB; ++i)
	{
		FERigidBody& obj = *rs.Object(i);
		if (obj.GetMaterialID() == nrb) return m_Data[ndata]->value(obj);
	}

	return val;
}

//-----------------------------------------------------------------------------
void ObjectDataRecord::SelectAllItems()
{
	int n = 0, i;
	for (i=0; i<m_pfem->Materials(); ++i)
	{
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(i));
		if (pm == 0) ++n;
	}

	if (n > 0)
	{
		m_item.resize(n);
		n = 0;
		for (i=0; i<m_pfem->Materials(); ++i)
		{
			FERigidMaterial* pm  = dynamic_cast<FERigidMaterial*>(m_pfem->GetMaterial(i));
			if (pm == 0)
			{
				m_item[n++] = i+1;
			}
		}
	}
}
