#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! FETractionLoad is a surface that has a constant (deformation independant)
//! traction force on it.
//!
class FETractionLoad : public FESurfaceLoad
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
	FETractionLoad(FEModel* pfem);

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

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
