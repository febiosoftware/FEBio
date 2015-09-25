#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! FEFluidTractionLoad is a fluid surface that has a constant
//! traction force on it.
//!
class FEFluidTractionLoad : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		vec3d	s[9];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FEFluidTractionLoad(FEModel* pfem);

	//! allocate storage
	void Create(int n);

	//! get a traction load BC
	LOAD& TractionLoad(int n) { return m_TC[n]; }

	//! calculate traction stiffness (there is none)
	void StiffnessMatrix(FESolver* psolver) {}

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

public:
	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

private:
	double	m_scale;	//!< magnitude of traction load
	vec3d	m_traction;	//!< traction vector

protected:
	vector<LOAD>	m_TC;		//!< traction boundary cards

	DECLARE_PARAMETER_LIST();
};
