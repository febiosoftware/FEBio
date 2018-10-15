#pragma once
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FESurface.h>
#include <FECore/tens3d.h>
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEMicroFlucSurface : public FESurface
{
public:
	//! constructor
	FEMicroFlucSurface(FEModel* fem);

	//! Initialization
	bool Init();

	//! copy data
	void CopyFrom(FEMicroFlucSurface& s);

public:
	vec3d SurfMicrofluc();

public:
	vec3d		m_Lm;	// Lagrange multipler microfluctuation
	vec3d		m_pv;	// "Pressure" vector
	vec3d		m_c;	// Microfluction across surface

	mat3d		m_Fm;	// Macroscopic deformation gradient
	tens3drs	m_Gm;	// Macroscopic deformation Hessian
};

//-----------------------------------------------------------------------------
// This class implements a constraint that tries to maintain the volume of the 
// enclosed space using an isochoric pressure.
class FE2OMicroConstraint : public FESurfaceConstraint
{
public:
	//! constructor
	FE2OMicroConstraint(FEModel* pfem);

	void Activate() override;
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void Serialize(DumpStream& ar) override;
	void CopyFrom(FENLConstraint* plc) override;

	// update state
	void Reset() override;
	void Update(const FETimeInfo& tp);

	FESurface* GetSurface() override;

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

public:
	FEMicroFlucSurface m_s;	//!< the bounding surface

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag

private:
	bool	m_binit;	//!< flag indicating whether the constraint is initialized

	int		m_dofX;
	int		m_dofY;
	int		m_dofZ;

	DECLARE_FECORE_CLASS();
};
