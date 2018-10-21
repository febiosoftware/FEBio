#include "stdafx.h"
#include "FETractionLoad.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_FECORE_CLASS(FETractionLoad, FESurfaceTraction)
	ADD_PARAMETER(m_scale   , "scale");
	ADD_PARAMETER(m_traction, "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceTraction(pfem)
{
	m_scale = 1.0;
	m_traction = vec3d(0, 0, 0);

	// Since the traction is deformation independent, we need to set the linear flag
	SetLinear(true);
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_traction.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
vec3d FETractionLoad::Traction(const FESurfaceMaterialPoint& pt)
{
	vec3d t = m_traction(pt)*m_scale;
	return t;
}
