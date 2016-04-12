#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! Surface that sustains a heat flux boundary condition
//!
class FEHeatFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		double	s[9];		// nodal scale factors
	};

public:
	//! constructor
	FEHeatFlux(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! get a heat flux load BC
	LOAD& HeatFlux(int n) { return m_FC[n]; }

	//! stiffness matrix
	void StiffnessMatrix(FESolver* psolver) {}
	
	//! residual
	void Residual(FEGlobalVector& R);

	//! serialization
	void Serialize(DumpStream& ar);

public:
	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

public:
	double	m_flux;	//!< heat flux

protected:
	vector<LOAD>	m_FC;

	DECLARE_PARAMETER_LIST();
};
