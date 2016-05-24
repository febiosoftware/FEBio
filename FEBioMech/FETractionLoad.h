#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FETractionLoad is a surface that has a constant (deformation independant)
//! traction force on it.
//!
class FETractionLoad : public FESurfaceLoad
{
public:
	//! constructor
	FETractionLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps);

	//! calculate traction stiffness (there is none)
	void StiffnessMatrix(const FETimePoint& tp, FESolver* psolver) {}

	//! calculate residual
	void Residual(const FETimePoint& tp, FEGlobalVector& R);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

private:
	double	m_scale;	//!< magnitude of traction load

protected:
	FESurfaceMap	m_TC;		//!< traction boundary cards

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
