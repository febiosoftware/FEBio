#include "stdafx.h"
#include "FEElasticDomain.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEElasticDomain::FEElasticDomain(FEModel* pfem) : m_pfem(pfem) 
{
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");

    m_dofRU = pfem->GetDOFIndex("Ru");
    m_dofRV = pfem->GetDOFIndex("Rv");
    m_dofRW = pfem->GetDOFIndex("Rw");
    
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
    
    m_dofSXP = pfem->GetDOFIndex("sxp");
    m_dofSYP = pfem->GetDOFIndex("syp");
    m_dofSZP = pfem->GetDOFIndex("szp");
    
    m_dofVX = pfem->GetDOFIndex("vx");
    m_dofVY = pfem->GetDOFIndex("vy");
    m_dofVZ = pfem->GetDOFIndex("vz");
    
    m_dofSVX  = pfem->GetDOFIndex("svx");
    m_dofSVY  = pfem->GetDOFIndex("svy");
    m_dofSVZ  = pfem->GetDOFIndex("svz");
    m_dofSVXP = pfem->GetDOFIndex("svxp");
    m_dofSVYP = pfem->GetDOFIndex("svyp");
    m_dofSVZP = pfem->GetDOFIndex("svzp");
    
    m_dofSAX  = pfem->GetDOFIndex("sax");
    m_dofSAY  = pfem->GetDOFIndex("say");
    m_dofSAZ  = pfem->GetDOFIndex("saz");
    m_dofSAXP = pfem->GetDOFIndex("saxp");
    m_dofSAYP = pfem->GetDOFIndex("sayp");
    m_dofSAZP = pfem->GetDOFIndex("sazp");
}
