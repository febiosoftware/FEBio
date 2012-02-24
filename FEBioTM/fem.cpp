#include "stdafx.h"
#include "fem.h"
#include "FEBioLib/FEBioImport.h"
#include "FEBioLib/FERigidBody.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "Document.h"
#include "FECore/FEException.h"

//-----------------------------------------------------------------------------
FEM::FEM()
{
	m_pTask = 0;
}

//-----------------------------------------------------------------------------
FEM::FEM(CTask* pt)
{
	m_pTask = pt;
}

//-----------------------------------------------------------------------------
void FEM::CheckInterruption()
{
	if (m_pTask && (m_pTask->GetStatus() == CTask::CANCELLED)) 
	{
		throw ExitRequest();
	}
	Fl::check();
}

//-----------------------------------------------------------------------------
static FEM* pfem_copy = 0;

//-----------------------------------------------------------------------------
void FEM::PushState()
{
	if (pfem_copy == 0) pfem_copy = new FEM;
	pfem_copy->ShallowCopy(*this);
}

//-----------------------------------------------------------------------------
void FEM::PopState()
{
	assert(pfem_copy);
	ShallowCopy(*pfem_copy);
	delete pfem_copy;
	pfem_copy = 0;
}

//-----------------------------------------------------------------------------
// This function is used when pushing the FEM state data. Since we don't need
// to copy all the data, this function only copies the data that needs to be 
// restored for a running restart.

// TODO: Shallow copy nonlinear constraints
void FEM::ShallowCopy(FEM& fem)
{
	int i;

	// copy time data
	m_ftime = fem.m_ftime;

	// copy the mesh
	m_mesh = fem.m_mesh;

	// copy rigid body data
	if (m_Obj.empty())
	{
		for (int i=0; i<(int) fem.m_Obj.size();	++i)
		{
			FERigidBody* prb = new FERigidBody(this);
			m_Obj.push_back(prb);
		}
	}
	assert(m_Obj.size() == fem.m_Obj.size());
	for (i=0; i<(int) m_Obj.size(); ++i) m_Obj[i]->ShallowCopy(fem.m_Obj[i]);

	// copy contact data
	if (m_CI.empty())
	{
		FEContactInterface* pci;
		for (int i=0; i<fem.ContactInterfaces(); ++i)
		{
			switch (fem.m_CI[i]->Type())
			{
			case FE_CONTACT_SLIDING    : pci = new FESlidingInterface  (this); break;
			case FE_FACET2FACET_SLIDING: pci = new FEFacet2FacetSliding(this); break;
			case FE_CONTACT_TIED       : pci = new FETiedInterface     (this); break;
			case FE_CONTACT_RIGIDWALL  : pci = new FERigidWallInterface(this); break;
			case FE_CONTACT_SLIDING2   : pci = new FESlidingInterface2 (this); break;
			case FE_PERIODIC_BOUNDARY  : pci = new FEPeriodicBoundary  (this); break;
			case FE_SURFACE_CONSTRAINT : pci = new FESurfaceConstraint (this); break;
			case FE_CONTACT_SLIDING3   : pci = new FESlidingInterface3 (this); break;
			default:
				assert(false);
			}
			m_CI.push_back(pci);
		}
	}
	assert(ContactInterfaces() == fem.ContactInterfaces());
	for (i=0; i<ContactInterfaces(); ++i) m_CI[i]->ShallowCopy(*fem.m_CI[i]);
}
