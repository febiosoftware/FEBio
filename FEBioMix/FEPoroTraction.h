#pragma once
#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! This boundary condition applies a poro-elastic normal traction on a surface
//!
class FEPoroNormalTraction : public FESurfaceLoad
{
public:
	//! constructor
	FEPoroNormalTraction(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetEffective(bool beff) { m_beffective = beff; }

	//! calculate pressure stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;

	//! unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm);

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
    bool    m_bshellb;      //!< flag for prescribing traction on shell bottom
	bool	m_beffective;	//!< effective or total normal traction

	// pressure boundary data
	FESurfaceMap	m_PC;		//!< pressure boundary cards

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofP;
    int	m_dofSX;
    int	m_dofSY;
    int	m_dofSZ;
    int	m_dofQ;

	DECLARE_FECORE_CLASS();
};
