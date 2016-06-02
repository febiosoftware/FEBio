#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! Surface that sustains a heat flux boundary condition
//!
class FEHeatFlux : public FESurfaceLoad
{
public:
	//! constructor
	FEHeatFlux(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! stiffness matrix
	void StiffnessMatrix(const FETimePoint& tp, FESolver* psolver) {}
	
	//! residual
	void Residual(const FETimePoint& tp, FEGlobalVector& R);

protected:
	double	m_flux;	//!< heat flux
	FESurfaceMap	m_FC;

	DECLARE_PARAMETER_LIST();
};
