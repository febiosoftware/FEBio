#include "stdafx.h"
#include "FERigidConnector.h"
#include "FERigidSystem.h"
#include "FERigidBody.h"
#include "FEMechModel.h"
#include "FERigidMaterial.h"
#include <FECore/FEMaterial.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidConnector, FENLConstraint);
	ADD_PARAMETER(m_nRBa, "body_a"        );
	ADD_PARAMETER(m_nRBb, "body_b"        );
END_FECORE_CLASS();

int FERigidConnector::m_ncount = 0;

//-----------------------------------------------------------------------------
FERigidConnector::FERigidConnector(FEModel* pfem) : FENLConstraint(pfem) 
{
    m_F = m_M = vec3d(0, 0, 0);
	m_rbA = m_rbB = 0;
};

//-----------------------------------------------------------------------------
FERigidConnector::~FERigidConnector() {}

//-----------------------------------------------------------------------------
bool FERigidConnector::Init()
{
	// do base class first
	if (FENLConstraint::Init() == false) return false;

	// When the rigid spring is read in, the ID's correspond to the rigid materials.
	// Now we want to make the ID's refer to the rigid body ID's
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_nRBa - 1));
	if (pm == nullptr)
	{
		felog.printbox("FATAL ERROR", "Rigid connector %d (spring) does not connect two rigid bodies\n", m_nID + 1);
		return false;
	}
	m_nRBa = pm->GetRigidBodyID();

	pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(m_nRBb - 1));
	if (pm == nullptr)
	{
		felog.printbox("FATAL ERROR", "Rigid connector %d (spring) does not connect two rigid bodies\n", m_nID);
		return false;
	}
	m_nRBb = pm->GetRigidBodyID();

	// get the actual rigid bodies
	FERigidSystem& rs = *fem.GetRigidSystem();
	m_rbA = rs.Object(m_nRBa);
	m_rbB = rs.Object(m_nRBb);

	return true;
}

//-----------------------------------------------------------------------------
void FERigidConnector::BuildMatrixProfile(FEGlobalMatrix& M)
{
    vector<int> lm(12);
                    
    int* lm1 = m_rbA->m_LM;
    int* lm2 = m_rbB->m_LM;
                    
    for (int j=0; j<6; ++j) lm[j  ] = lm1[j];
    for (int j=0; j<6; ++j) lm[j+6] = lm2[j];
    M.build_add(lm);
}

//-----------------------------------------------------------------------------
void FERigidConnector::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	if (ar.IsSaving())
	{
        ar << m_nID;
        ar << m_nRBa << m_nRBb;
        ar << m_F << m_M;
	}
	else
	{
        ar >> m_nID;
        ar >> m_nRBa >> m_nRBb;
        ar >> m_F >> m_M;

		// get the actual rigid bodies
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidSystem& rs = *fem.GetRigidSystem();
		m_rbA = rs.Object(m_nRBa);
		m_rbB = rs.Object(m_nRBb);
	}
}
