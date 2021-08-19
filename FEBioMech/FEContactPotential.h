#pragma once
#include <FEBioMech/FEContactInterface.h>
#include <FEBioMech/FEContactSurface.h>
#include <set>
using namespace std;

class FEContactPotentialSurface : public FEContactSurface
{
public:
	class Data : public FEContactMaterialPoint
	{
	public:
		vec3d	m_tc;
	};

public:
	FEContactPotentialSurface(FEModel* fem);

	void GetContactTraction(int nelem, vec3d& tc) override;

	FEMaterialPoint* CreateMaterialPoint() override;
};

typedef FEContactPotentialSurface::Data FECPContactPoint;

class FEContactPotential : public FEContactInterface
{
public:
	FEContactPotential(FEModel* fem);

	// -- From FESurfacePairConstraint
public:
	//! return the primary surface
	FESurface* GetPrimarySurface() override;

	//! return the secondary surface
	FESurface* GetSecondarySurface() override;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	bool UseNodalIntegration() override;

	// Build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// update
	void Update() override;

	// init
	bool Init() override;

	// -- from FEContactInterface
public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	// Evaluates the contriubtion to the stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

protected:
	void ElementForce(FESurfaceElement& el1, FESurfaceElement& el2, vector<double>& fe);
	void ElementStiffness(FESurfaceElement& el1, FESurfaceElement& el2, matrix& ke);

	double PotentialDerive(double r);
	double PotentialDerive2(double r);

protected:
	FEContactPotentialSurface	m_surf1;
	FEContactPotentialSurface	m_surf2;

protected:
	double	m_kc;
	double	m_p;
	double	m_Rin;
	double	m_Rout;
	double	m_wtol;

	double	m_c1, m_c2;

	vector<	set<FESurfaceElement*> >			m_activeElements;
	vector< set<FESurfaceElement*> >	m_elemNeighbors;

	DECLARE_FECORE_CLASS();
};

