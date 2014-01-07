#pragma once
#include "FECore/FESurfaceLoad.h"

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceLoad
{
public:
	struct LOAD
	{
		LOAD();
		double	s[8];		// nodal scale factors
		int		lc;			// load curve
	};

public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! allocate storage
	void Create(int n);

	//! get a pressure load BC
	LOAD& PressureLoad(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpFile& ar);

	//! set the linear flag
	void SetLinear(bool blinear) { m_blinear = blinear; }

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

public:
	//! set an attribute of the surface load
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	bool			m_blinear;	//!< pressure load type (linear or nonlinear)
	double			m_pressure;	//!< pressure value
	vector<LOAD>	m_PC;		//!< pressure load cards

	DECLARE_PARAMETER_LIST();
};
