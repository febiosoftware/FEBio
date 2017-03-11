#pragma once
#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidTractionLoad is a fluid surface that has a prescribed
//! viscous traction vector on it.
//!
class FEFluidTractionLoad : public FESurfaceLoad
{
public:
	//! constructor
	FEFluidTractionLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! calculate traction stiffness (there is none)
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

private:
	double			m_scale;	//!< magnitude of traction load
	FESurfaceMap	m_TC;		//!< traction boundary cards

	int		m_dofVX;
	int		m_dofVY;
	int		m_dofVZ;

	DECLARE_PARAMETER_LIST();
};
