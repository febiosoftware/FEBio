#pragma once
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
class FEVolumeSurface : public FESurface
{
public:
	//! constructor
	FEVolumeSurface(FEMesh* pm);

	//! Initialization
	bool Init();

	//! copy data
	void CopyFrom(FEVolumeSurface& s);

	//! serialization
	void Serialize(DumpStream& ar);

public:
	double Volume();

public:
	double	m_Lp;	//!< Lagrange multipler pressure
	double	m_p;	//!< applied pressure (= Lp + eps*DV)
	double	m_V0;	//!< Initial volume
	double	m_Vt;	//!< current volume
};

//-----------------------------------------------------------------------------
// This class implements a constraint that tries to maintain the volume of the 
// enclosed space using an isochoric pressure.
class FEVolumeConstraint : public FESurfaceConstraint
{
public:
	//! constructor
	FEVolumeConstraint(FEModel* pfem);

	void Activate();
	void Residual(FEGlobalVector& R, const FETimeInfo& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);
	bool Augment(int naug, const FETimeInfo& tp);
	void Serialize(DumpStream& ar);
	void CopyFrom(FENLConstraint* plc);

	// update state
	void Reset();
	void Update(int niter, const FETimeInfo& tp);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// get the surface
	FESurface* GetSurface() override;

public:
	FEVolumeSurface m_s;	//!< the bounding surface

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag

private:
	bool	m_binit;	//!< flag indicating whether the constraint is initialized

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
