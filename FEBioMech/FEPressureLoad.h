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
		double	s[9];		// nodal scale factors
	};

public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! get a pressure load BC
	LOAD& PressureLoad(int n) { return m_PC[n]; }

	//! calculate pressure stiffness
	void StiffnessMatrix(FESolver* psolver);

	//! calculate residual
	void Residual(FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpStream& ar);

	//! set the linear flag
	void SetLinear(bool blinear) { m_blinear = blinear; }

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	//! set an attribute of the surface load
	// NOTE: To be removed in FEBio 3.0
	bool SetAttribute(const char* szatt, const char* szval);

	//! set an attribute of a surface facet
	// NOTE: To be removed in FEBio 3.0
	bool SetFacetAttribute(int nface, const char* szatt, const char* szval);

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);

	void SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);
	void UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn);

	//! Calculates external pressure forces
	void PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	void LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn);

protected:
	bool			m_blinear;	//!< pressure load type (linear or nonlinear)
	double			m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation
	vector<LOAD>	m_PC;		//!< pressure load cards

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
