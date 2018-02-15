#include "stdafx.h"
#include "FEElasticDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEElasticDomain::FEElasticDomain(FEModel* pfem) : m_pfem(pfem) 
{
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");

    m_dofU = pfem->GetDOFIndex("u");
    m_dofV = pfem->GetDOFIndex("v");
    m_dofW = pfem->GetDOFIndex("w");
    
    m_dofUP = pfem->GetDOFIndex("up");
    m_dofVP = pfem->GetDOFIndex("vp");
    m_dofWP = pfem->GetDOFIndex("wp");
    
    m_dofRU = pfem->GetDOFIndex("Ru");
	m_dofRV = pfem->GetDOFIndex("Rv");
	m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    
    m_dofVU  = pfem->GetDOFIndex("vu");
    m_dofVV  = pfem->GetDOFIndex("vv");
    m_dofVW  = pfem->GetDOFIndex("vw");
    m_dofVUP = pfem->GetDOFIndex("vup");
    m_dofVVP = pfem->GetDOFIndex("vvp");
    m_dofVWP = pfem->GetDOFIndex("vwp");
    
    m_dofAU  = pfem->GetDOFIndex("au");
    m_dofAV  = pfem->GetDOFIndex("av");
    m_dofAW  = pfem->GetDOFIndex("aw");
    m_dofAUP = pfem->GetDOFIndex("aup");
    m_dofAVP = pfem->GetDOFIndex("avp");
    m_dofAWP = pfem->GetDOFIndex("awp");
}
