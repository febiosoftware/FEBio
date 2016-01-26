#pragma once
#include <FECore/FENLConstraint.h>
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
	void Serialize(DumpFile& ar);

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
class FEVolumeConstraint : public FENLConstraint
{
public:
	//! constructor
	FEVolumeConstraint(FEModel* pfem);

	void Activate();
	void Residual(FEGlobalVector& R, const FETimePoint& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
	bool Augment(int naug, const FETimePoint& tp);
	void Serialize(DumpFile& ar);
	void ShallowCopy(DumpStream& dmp, bool bsave);
	void CopyFrom(FENLConstraint* plc);

	// update state
	void Reset();
	void Update(const FETimePoint& tp);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M);

	FESurface* GetSurface(const char* sz);

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
