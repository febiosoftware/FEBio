#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! This boundary condition sustains a fluid flux on a surface
//!
class FEFluidFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		double	s[9];		// nodal scale factors
	};

public:
	//! constructor
	FEFluidFlux(FEModel* pfem);

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetMixture(bool bmix) { m_bmixture = bmix; }

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! get a flux BC
	LOAD& FluidFlux(int n) { return m_PC[n]; }

	//! calculate flux stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpStream& ar);

	//! unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! set an attribute of the surface load
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt, bool mixture);

	//! calculate stiffness for an element, for steady-state analysis
	void FluxStiffnessSS(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates volumetric flow rate due to flux, for steady-state analysis
	bool FlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower), for steady-state analysis
	bool LinearFlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt, bool mixture);

protected:
	double	m_flux;			//!< fluid flux
	bool	m_bmixture;		//!< mixture velocity or relative fluid flux
	bool	m_blinear;		//!< type (linear or nonlinear)

	// Fluid flux boundary data
	vector<LOAD>	m_PC;		//!< fluid flux boundary cards

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofP;
	int	m_dofVX;
	int	m_dofVY;
	int	m_dofVZ;

	DECLARE_PARAMETER_LIST();
};
