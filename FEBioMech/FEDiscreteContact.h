#pragma once
#include <FECore/FENLConstraint.h>
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FEDiscreteSet;

//-----------------------------------------------------------------------------
class FEDiscreteContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEDiscreteContactSurface(FEModel* fem);

	//! Initialization
	bool Init();
};

//-----------------------------------------------------------------------------
class FEDiscreteContact : public FENLConstraint
{
	struct NODE
	{
		int					nid;	//!< (local) node ID
		FESurfaceElement*	pe;		//!< master surface element
		double				gap;	//!< gap distance
		double				Lm;		//!< Lagrange multiplier
		vec3d				nu;		//!< normal at master projection
		vec3d				q;		//!< projection point
		double				proj[2];	//!< iso-parametric coordinates of projection point
	};

public:
	FEDiscreteContact(FEModel* pfem);

public:
	bool Init();
	void Activate();

	void Residual(FEGlobalVector& R, const FETimePoint& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
	bool Augment(int naug, const FETimePoint& tp);
	void BuildMatrixProfile(FEGlobalMatrix& M);
	void Update(const FETimePoint& tp);

	void SetDiscreteSet(FEDiscreteSet* pset);

	FESurface* GetSurface(const char* sz) { return &m_surf; }

protected:
	void ProjectSurface(bool bupseg);
	void ContactNodalForce    (NODE& nodeData, FESurfaceElement& mel, vector<double>& fe);
	void ContactNodalStiffness(NODE& nodeData, FESurfaceElement& mel, matrix& ke);

protected:
	FEDiscreteContactSurface	m_surf;
	vector<NODE>	m_Node;
	double	m_normg0;
	bool	m_bfirst;

protected:
	bool	m_blaugon;	//!< augmentation flag
	double	m_altol;	//!< augmentation tolerance
	double	m_penalty;	//!< penalty parameter
	double	m_gaptol;	//!< gap tolerance
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations
	int		m_nsegup;	//!< number of segment updates (or zero)

	DECLARE_PARAMETER_LIST();
};
