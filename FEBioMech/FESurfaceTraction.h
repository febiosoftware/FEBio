#pragma once
#include <FECore/FESurfaceLoad.h>

class FESurfaceTraction : public FESurfaceLoad
{
public:
	// Constructor
	FESurfaceTraction(FEModel* fem);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! helper function for getting nodal coordinates
	void GetNodalCoordinates(FESurfaceElement& el, vec3d* re);

	//! helper function for getting nodal coordinates
	void GetReferenceNodalCoordinates(FESurfaceElement& el, vec3d* re);

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;

	//! calculate stiffness
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! Is the traction linear?
	bool IsLinear() const { return m_blinear; }

	//! Set the traction load linear or not
	void SetLinear(bool b) { m_blinear = b; }

protected:
	//! element contribution to residual
	virtual void ElementResidual(FESurfaceElement& el, std::vector<double>& fe);

	//! element contribution to stiffness
	virtual void ElementStiffness(FESurfaceElement& el, matrix& ke) {}

	//! traction at material point
	virtual vec3d Traction(const FESurfaceMaterialPoint& mp) { return vec3d(0, 0, 0);  }

protected:
	// degrees of freedom
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	int	m_dofSX;
	int	m_dofSY;
	int	m_dofSZ;

	bool	m_blinear;	//!< is the load linear (i.e. it will be calculated in the reference frame and assummed deformation independent)
	bool	m_bshellb;	//!< flag for prescribing pressure on shell bottom

	DECLARE_FECORE_CLASS();
};
