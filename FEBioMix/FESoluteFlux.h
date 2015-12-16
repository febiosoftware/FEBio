#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! The flux surface is a surface domain that sustains a solute flux boundary
//! condition
//!
class FESoluteFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		double	s[9];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FESoluteFlux(FEModel* pfem);
	
	//! allocate storage
	void Create(int n);

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetSolute(int isol) { m_isol = isol; }
	
	//! get a flux BC
	LOAD& SoluteFlux(int n) { return m_PC[n]; }
	
	//! calculate flux stiffness
	void StiffnessMatrix(FESolver* psolver);
	
	//! calculate residual
	void Residual(FEGlobalVector& R);
	
	//! serialize data
	void Serialize(DumpFile& ar);

	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! set an attribute of the surface load
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn, double dt);
	
	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt);
	
	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn, double dt);
	
protected:
	double	m_flux;		//!< flux magnitude
	bool	m_blinear;	//!< linear or not (true is non-follower, false is follower)
	int		m_isol;		//!< solute index

	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofC;

	// solute flux boundary data
	vector<LOAD>	m_PC;		//!< solute flux boundary cards

	DECLARE_PARAMETER_LIST();
};
