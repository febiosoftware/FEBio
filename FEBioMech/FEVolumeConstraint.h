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

	bool Init();
	void Residual(FEGlobalVector& R);
	void StiffnessMatrix(FESolver* psolver);
	bool Augment(int naug);
	void Serialize(DumpFile& ar);

	// update state
	void Reset();
	void Update();

	FESurface* GetSurface(const char* sz);

public:
	FEVolumeSurface m_s;	//!< the bounding surface

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag

	DECLARE_PARAMETER_LIST();
};
