#pragma once
#include <FECore/FENLConstraint.h>
#include <FECore/FESurface.h>
#include <FECore/tens3d.h>
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEMicroFlucSurface : public FESurface
{
public:
	//! constructor
	FEMicroFlucSurface(FEMesh* pm);

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
class FE2OMicroConstraint : public FENLConstraint
{
public:
	//! constructor
	FE2OMicroConstraint(FEModel* pfem);

	void Activate();
	void Residual(FEGlobalVector& R, const FETimePoint& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
	bool Augment(int naug, const FETimePoint& tp);
	void Serialize(DumpStream& ar);
	void CopyFrom(FENLConstraint* plc);

	// update state
	void Reset();
	void Update(const FETimePoint& tp);

	FESurface* GetSurface(const char* sz);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M);

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

	DECLARE_PARAMETER_LIST();
};
