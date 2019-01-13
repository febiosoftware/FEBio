#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! The flux surface is a surface domain that sustains a solute flux boundary
//! condition
//!
class FECORE_API FESoluteFlux : public FESurfaceLoad
{
public:
	//! constructor
	FESoluteFlux(FEModel* pfem);
	
	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	void SetLinear(bool blinear) { m_blinear = blinear; }

	void SetSolute(int isol) { m_isol = isol; }
	
	//! calculate flux stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;
	
	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
	
	void UnpackLM(FEElement& el, vector<int>& lm);

protected:
	//! calculate stiffness for an element
	void FluxStiffness(FESurfaceElement& el, matrix& ke, double dt);
	
	//! Calculates volumetric flow rate due to flux
	bool FlowRate(FESurfaceElement& el, vector<double>& fe, double dt);
	
	//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
	bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, double dt);
	
protected:
	FEParamDouble	m_flux;		//!< flux scale factor magnitude
	bool	m_blinear;	//!< linear or not (true is non-follower, false is follower)
    bool    m_bshellb;  //!< flag for prescribing flux on shell bottom
	int		m_isol;		//!< solute index

protected:
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofC;
    int	m_dofSX;
    int	m_dofSY;
    int	m_dofSZ;
    int	m_dofD;

	DECLARE_FECORE_CLASS();
};
