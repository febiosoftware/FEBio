#include "stdafx.h"
#include "fem.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioMech/FERigidBody.h"
#include "FEBioMech/FESlidingInterface.h"
#include "FEBioMech/FETiedInterface.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEFacet2FacetSliding.h"
#include "FEBioMech/FEFacet2FacetTied.h"
#include "FEBioMix/FESlidingInterface2.h"
#include "FEBioMix/FESlidingInterface3.h"
#include "FEBioMech/FEPeriodicBoundary.h"
#include "FEBioMech/FESurfaceConstraint.h"
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
	// get a lock
	Fl::lock();
	if (m_pTask && (m_pTask->GetStatus() == CTask::CANCELLED)) 
	{
		// release lock
		Fl::unlock();
		throw ExitRequest();
	}
	// release lock
	Fl::unlock();
}
