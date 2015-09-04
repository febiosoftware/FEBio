#include "stdafx.h"
#include "FESolidAnalysis.h"
#include "FERigidMaterial.h"
#include "FEContactInterface.h"
#include "FEUncoupledMaterial.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FESolidAnalysis::FESolidAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_SOLID) 
{

}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FESolidAnalysis::Activate()
{
	// initialize base class data
	FEAnalysis::Activate();

	if (m_nanalysis == FE_DYNAMIC)
	{
		FEMesh& mesh = m_fem.GetMesh();
		FERigidSystem& rigid = *m_fem.GetRigidSystem();

		// set the initial velocities of all rigid nodes
		for (int i=0; i<mesh.Nodes(); ++i)
		{
			FENode& n = mesh.Node(i);
			if (n.m_rid >= 0)
			{
				FERigidBody& rb = *rigid.Object(n.m_rid);
				vec3d V = rb.m_vt;
				vec3d W = rb.m_wt;
				vec3d r = n.m_rt - rb.m_rt;

				vec3d v = V + (W ^ r); 
				n.m_vp = n.m_vt = v;

				vec3d a = (W ^ V)*2.0 + (W ^ (W ^ r));
				n.m_ap = n.m_at = a;
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESolidAnalysis::Deactivate()
{
	FEAnalysis::Deactivate();

	// store the current rigid body reaction forces
	FERigidSystem& rigid = *m_fem.GetRigidSystem();
	for (int i=0; i<rigid.Objects(); ++i)
	{
		FERigidBody& RB = *rigid.Object(i);
		RB.m_Fp = RB.m_Fr;
		RB.m_Mp = RB.m_Mr;
	}
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FEExplicitSolidAnalysis::Activate()
{
	// initialize base class data
	return FEAnalysis::Activate();
}

//-----------------------------------------------------------------------------
//! This function is called before the analysis is solved and initializes all
//! analysis data, such as determine active boundary conditions, initializes
//! equation numbers (the latter is actually done by the FESolver class).
bool FELinearSolidAnalysis::Activate()
{
	// initialize base class data
	return FEAnalysis::Activate();
}
