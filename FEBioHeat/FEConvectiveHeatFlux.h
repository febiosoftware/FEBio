#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! Surface that sustains a convective heat flux boundary condition
//!
class FEConvectiveHeatFlux : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		double	hc;			// heat transfer coefficient
		double	s[9];		// nodal scale factors
	};

public:
	//! constructor
	FEConvectiveHeatFlux(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! get a heat flux load BC
	LOAD& HeatFlux(int n) { return m_FC[n]; }

	//! stiffness matrix
	void StiffnessMatrix(const FETimePoint& tp, FESolver* psolver);
	
	//! residual
	void Residual(const FETimePoint& tp, FEGlobalVector& R);

	//! serialization
	void Serialize(DumpStream& ar);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

private:
	double	m_hc;		//!< heat transfer coefficient
	double	m_Ta;		//!< ambient temperature

protected:
	void ElementStiffness(FESurfaceElement& el, matrix& ke, double hc);

protected:
	vector<LOAD>	m_FC;

	DECLARE_PARAMETER_LIST();
};
