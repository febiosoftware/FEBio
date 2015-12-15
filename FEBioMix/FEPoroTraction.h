#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! This boundary condition applies a poro-elastic normal traction on a surface
//!
class FEPoroNormalTraction : public FESurfaceLoad
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
	FEPoroNormalTraction(FEModel* pfem);

	//! allocate storage
	void Create(int n);

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetEffective(bool beff) { m_beffective = beff; }

	//! get a pressure load BC
	LOAD& NormalTraction(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpFile& ar);

	//! unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! set an attribute of the surface load
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	//! calculate stiffness for an element
	void TractionStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn, bool effective, bool bsymm);

	//! Calculates external pressure forces
	bool TractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearTractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	double	m_traction;		//!< traction value
	bool	m_blinear;		//!< linear or not (true is non-follower, false is follower)
	bool	m_beffective;	//!< effective or total normal traction

	// pressure boundary data
	vector<LOAD>	m_PC;		//!< pressure boundary cards

	DECLARE_PARAMETER_LIST();
};
