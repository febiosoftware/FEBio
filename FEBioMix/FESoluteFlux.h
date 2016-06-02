#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! The flux surface is a surface domain that sustains a solute flux boundary
//! condition
//!
class FESoluteFlux : public FESurfaceLoad
{
public:
	//! constructor
	FESoluteFlux(FEModel* pfem);
	
	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetSolute(int isol) { m_isol = isol; }
	
	//! calculate flux stiffness
	void StiffnessMatrix(const FETimePoint& tp, FESolver* psolver);
	
	//! calculate residual
	void Residual(const FETimePoint& tp, FEGlobalVector& R);
	
	void UnpackLM(FEElement& el, vector<int>& lm);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt);
	
	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt);
	
	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt);
	
protected:
	double	m_flux;		//!< flux scale factor magnitude
	bool	m_blinear;	//!< linear or not (true is non-follower, false is follower)
	int		m_isol;		//!< solute index
	FESurfaceMap	m_PC;		//!< solute flux boundary cards

protected:
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofC;

	DECLARE_PARAMETER_LIST();
};
