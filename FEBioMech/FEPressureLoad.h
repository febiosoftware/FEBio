#pragma once
#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceLoad
{
public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! calculate pressure stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver);

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R);

	//! serialize data
	void Serialize(DumpStream& ar);

	//! set the linear flag
	void SetLinear(bool blinear) { m_blinear = blinear; }

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

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
    bool            m_bshellb; //!< flag for prescribing pressure on shell bottom
	double			m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation
	FESurfaceMap	m_PC;		//!< pressure scale factors

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
    int	m_dofU;
    int	m_dofV;
    int	m_dofW;

	DECLARE_PARAMETER_LIST();
};
