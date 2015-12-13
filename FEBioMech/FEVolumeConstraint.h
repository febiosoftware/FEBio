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

	FESurface* GetSurface(const char* sz);

public:
	FEVolumeSurface m_s;	//!< the bounding surface

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag

private:
	bool	m_binit;	//!< flag indicating whether the constraint is initialized

	DECLARE_PARAMETER_LIST();
};
