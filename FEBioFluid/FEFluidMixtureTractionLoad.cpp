//
//  FEFluidMixtureTractionLoad.cpp
//  FEBioFluid
//
//  Created by Jay Shim on 3/2/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFluidMixtureTractionLoad.h"
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/FEMesh.h>
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidMixtureTractionLoad, FESurfaceLoad)
ADD_PARAMETER(m_scale, "scale"   );
ADD_PARAMETER(m_TC   , "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidMixtureTractionLoad::FEFluidMixtureTractionLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_scale = 1.0;
    m_TC = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidMixtureTractionLoad::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_TC.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
bool FEFluidMixtureTractionLoad::Init()
{
    m_dof.Clear();
    if (m_dof.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT)) == false) return false;
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidMixtureTractionLoad::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_psurf->LoadVector(R, m_dof, true, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
        
        // fluid traction
        vec3d t = m_TC(mp)*m_scale;
        vec3d f = t*((mp.dxr ^ mp.dxs).norm());
        
        double H = dof_a.shape;
        fa[0] = H * f.x;
        fa[1] = H * f.y;
        fa[2] = H * f.z;
    });
}

//-----------------------------------------------------------------------------
//! calculate traction stiffness (there is none)
void FEFluidMixtureTractionLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    
}
