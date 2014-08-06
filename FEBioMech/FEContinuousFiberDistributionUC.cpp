//
//  FEContinuousFiberDistributionUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEContinuousFiberDistributionUC.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEContinuousFiberDistributionUC, FEUncoupledMaterial)
    ADD_PARAMETER(m_a, FE_PARAM_VEC3D, "a");
    ADD_PARAMETER(m_d, FE_PARAM_VEC3D, "d");
END_PARAMETER_LIST();



//-----------------------------------------------------------------------------
int FEContinuousFiberDistributionUC::Properties()
{
	return 3;
}

//-----------------------------------------------------------------------------
//! get a specific material property
FECoreBase* FEContinuousFiberDistributionUC::GetProperty(int i)
{
	switch(i)
	{
        case 0: return m_pFmat; break;
        case 1: return m_pFDD; break;
        case 2: return m_pFint; break;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEContinuousFiberDistributionUC::FindPropertyIndex(const char* szname)
{
	if      (strcmp(szname, "fibers"      ) == 0) return 0;
	else if (strcmp(szname, "distribution") == 0) return 1;
	else if (strcmp(szname, "scheme"      ) == 0) return 2;
	else return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEContinuousFiberDistributionUC::SetProperty(int i, FECoreBase* pm)
{
	switch(i)
	{
        case 0:
		{
			m_pFmat = dynamic_cast<FEElasticFiberMaterialUC*>(pm);
			if ((m_pFmat == 0) || (m_pFmat->IsRigid())) return false;
			return true;
		}
            break;
        case 1:
		{
			m_pFDD = dynamic_cast<FEFiberDensityDistribution*>(pm);
			if (m_pFDD == 0) return false;
			return true;
		}
            break;
        case 2:
		{
			m_pFint = dynamic_cast<FEFiberIntegrationSchemeUC*>(pm);
			if (m_pFint == 0) return false;
			return true;
		}
            break;
	}
	return false;
}

//-----------------------------------------------------------------------------
void FEContinuousFiberDistributionUC::Init()
{
    // set parent materials
    m_pFmat->SetParent(this);
    m_pFDD->SetParent(this);
    m_pFint->SetParent(this);
    
    // propagate pointers to fiber material and density distribution
    // to fiber integration scheme
    m_pFint->m_pFmat = m_pFmat;
    m_pFint->m_pFDD = m_pFDD;
    
    // initialize fiber integration scheme
    m_pFint->Init();
    
    // normalize axes vectors and evaluate local orientation
    vec3d a = m_a; a.unit();
    vec3d d = m_d; d.unit();
    vec3d c = a ^ d; c.unit();
    vec3d b = c ^ a;
	m_Q[0][0] = a.x; m_Q[0][1] = b.x; m_Q[0][2] = c.x;
	m_Q[1][0] = a.y; m_Q[1][1] = b.y; m_Q[1][2] = c.y;
	m_Q[2][0] = a.z; m_Q[2][1] = b.z; m_Q[2][2] = c.z;
}

//-----------------------------------------------------------------------------
FEParam* FEContinuousFiberDistributionUC::GetParameter(const ParamString& s)
{
	if (s.count() == 1) return FEMaterial::GetParameter(s);
    
	if      (s == "fibers"      ) return m_pFmat->GetParameter(s.next());
	else if (s == "distribution") return m_pFDD ->GetParameter(s.next());
	else if (s == "scheme"      ) return m_pFint->GetParameter(s.next());
	return 0;
}
