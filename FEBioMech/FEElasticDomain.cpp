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
    
    m_dofRU = pfem->GetDOFIndex("Ru");
	m_dofRV = pfem->GetDOFIndex("Rv");
	m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofAU  = pfem->GetDOFIndex("au");
    m_dofAV  = pfem->GetDOFIndex("av");
    m_dofAW  = pfem->GetDOFIndex("aw");
}
