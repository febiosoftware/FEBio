#pragma once
#include "FESurfaceTraction.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! FETractionLoad is a surface that has a constant (deformation independant)
//! traction force on it.
//!
class FETractionLoad : public FESurfaceTraction
{
public:
	//! constructor
	FETractionLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

protected:
	//! calculate traction
	vec3d Traction(const FESurfaceMaterialPoint& el) override;

protected:
	double			m_scale;	//!< scale factor for traction
	FEParamVec3		m_traction;	//!< vector traction

	DECLARE_FECORE_CLASS();
};
