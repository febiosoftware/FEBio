#include "stdafx.h"
#include "FEBiphasicDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEBiphasicDomain::FEBiphasicDomain(FEModel* pfem) : FEElasticDomain(pfem)
{
    m_pMat = 0;
    m_dofP = pfem->GetDOFIndex("p");
    m_dofQ = pfem->GetDOFIndex("q");
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
}
