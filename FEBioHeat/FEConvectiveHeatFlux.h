#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! Surface that sustains a convective heat flux boundary condition
//!
class FEConvectiveHeatFlux : public FESurfaceLoad
{
public:
	//! constructor
	FEConvectiveHeatFlux(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! stiffness matrix (TODO: obsolete interface. Remove it.)
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) { assert(false); }

	//! stiffness matrix (new interface)
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp);
	
	//! residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R);

protected:
	void ElementStiffness(FESurfaceElement& el, matrix& ke, double hc);

private:
	double	m_hc;		//!< heat transfer coefficient
	double	m_Ta;		//!< ambient temperature
	FESurfaceMap	m_FC;

	DECLARE_PARAMETER_LIST();
};
