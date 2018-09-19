#pragma once
#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>

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
	void SetSurface(FESurface* ps) override;

	//! calculate pressure stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;

	//! serialize data
	void Serialize(DumpStream& ar) override;

	//! set the linear flag
	void SetLinear(bool blinear) { m_blinear = blinear; }

	//! Check if this is a linear force or not
	bool IsLinear() { return m_blinear; }

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

protected:
	//! calculate stiffness for an element
	void PressureStiffness(FESurfaceElement& el, matrix& ke);

	void SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);
	void UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);

	//! Calculates external pressure forces
	void PressureForce(FESurfaceElement& el, vector<double>& fe);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	void LinearPressureForce(FESurfaceElement& el, vector<double>& fe);

protected:
	bool			m_blinear;	//!< pressure load type (linear or nonlinear)
    bool            m_bshellb; //!< flag for prescribing pressure on shell bottom
	FEParamDouble	m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation
	bool			m_bstiff;	//!< use stiffness or not

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
    int	m_dofSX;
    int	m_dofSY;
    int	m_dofSZ;

	DECLARE_PARAMETER_LIST();
};
