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
	};

public:
	//! constructor
	FEFluidTractionLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! get a traction load BC
	LOAD& TractionLoad(int n) { return m_TC[n]; }

	//! calculate traction stiffness (there is none)
	void StiffnessMatrix(const FETimePoint& tp, FESolver* psolver) {}

	//! calculate residual
	void Residual(const FETimePoint& tp, FEGlobalVector& R);

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

private:
	double	m_scale;	//!< magnitude of traction load
	vec3d	m_traction;	//!< traction vector

protected:
	vector<LOAD>	m_TC;		//!< traction boundary cards

	int		m_dofVX;
	int		m_dofVY;
	int		m_dofVZ;

	DECLARE_PARAMETER_LIST();
};
