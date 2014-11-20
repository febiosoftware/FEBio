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
		int		lc;			// load curve describing ambient temperature
	};

public:
	//! constructor
	FEConvectiveHeatFlux(FEModel* pfem);

	//! allocate storage
	void Create(int n);

	//! get a heat flux load BC
	LOAD& HeatFlux(int n) { return m_FC[n]; }

	//! stiffness matrix
	void StiffnessMatrix(FESolver* psolver);
	
	//! residual
	void Residual(FEGlobalVector& R);

	//! serialization
	void Serialize(DumpFile& ar);

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
