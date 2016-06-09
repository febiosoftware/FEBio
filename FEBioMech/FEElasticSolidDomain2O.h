#pragma once
#include "FEElasticSolidDomain.h"
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FESurface;

//-----------------------------------------------------------------------------
// This class implements a discontinuous-Galerkin formulation for gradient elasticity
class FEElasticSolidDomain2O : public FEElasticSolidDomain
{
protected:
	// Helper class for evaluating the discrete-Galerkin contribution
	// It stores the data needed for evaluating the integrals over the
	// internal surface.
	class FEInternalSurface2O 
	{
	public:
		struct Data
		{
			FEMaterialPoint*	m_pt[2];	//!< material point data for evaluating stresses
			vec3d				ksi[2];		//!< local element coordinates
			tens3drs			Qavg;		//!< average stress across interface
			tens5d				H[2];		//!< H stiffness
			tens6d				J[2];		//!< J stiffness
			tens6d				J0[2];		//!< initial J stiffness
			tens5d				H0[2];		//!< initial H stiffness
			tens6d 				J0avg;		//!< average initial higher order stiffess across interface
			mat3d				DgradU;		//!< displacement gradient jump across interface
		};

	public:
		FEInternalSurface2O ();

		// initialize the data structure
		bool Initialize(FEElasticSolidDomain2O* dom);

		int Elements() const { return m_ps->Elements(); }

		FESurfaceElement& Element(int i) { return m_ps->Element(i); }

		Data& GetData(int i) { return m_data[i]; }

		FESurface* GetSurface() { return m_ps; }

		double GetElementSize() const { return m_h; }

	private:
		FESurface*		m_ps;
		vector<Data>	m_data;
		double			m_h;	//!< element size
	};

public:
	//! constructor
	FEElasticSolidDomain2O(FEModel* pfem);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! build the matrix profile
	//! (overridden from FEDomain)
	void BuildMatrixProfile(FEGlobalMatrix& M);

	//! overridden from FEElasticSolidDomain
	void Update(const FETimeInfo& tp);

public:
	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! evaluate internal element forces
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! calculates the global stiffness matrix for this domain
	//! (overridden from FEElasticSolidDomain)
	void StiffnessMatrix(FESolver* psolver);

protected:
	// discontinuous-Galerkin contribution to residual
	void InternalForcesDG1(FEGlobalVector& R);
	void InternalForcesDG2(FEGlobalVector& R);
	void InternalForcesDG3(FEGlobalVector& R);

	void ElementInternalForce_PF(FESolidElement& el, vector<double>& fe);
	void ElementInternalForce_QG(FESolidElement& el, vector<double>& fe);

	// --- S T I F F N E S S ---
	//! calculates the solid element stiffness matrix
	void ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke);

	//! contributions from discontinuous Galerkin formulation
	void StiffnessMatrixDG(FESolver* psolver);
	void ElementStiffnessMatrixDG1(FESurfaceElement& el, FEInternalSurface2O::Data* pdata, matrix& ke);
	void ElementStiffnessMatrixDG2(FESurfaceElement& el, FEInternalSurface2O::Data* pdata, matrix& ke);
	void ElementStiffnessMatrixDG3(FESurfaceElement& el, FEInternalSurface2O::Data* pdata, matrix& ke);

private:
	void UpdateElementStress(int iel, double dt);
	void UpdateInternalSurfaceStresses();
	void UpdateKinematics();

public:
	// calculate gradient of deformation gradient
	void defhess(FESolidElement &el, int n, tens3drs &G);
	void defhess(FESolidElement &el, double r, double s, double t, tens3drs &G);

	// calculates derivatives of shape functions
	void shape_gradient(const FESolidElement& el, int n, vec3d* G);
	void shape_gradient(const FESolidElement& el, double r, double s, double t, vec3d* G);

	// Calculates second derivative of shape function N[node]
	void shape_gradient2(const FESolidElement& el, vec3d* X, int n, mat3d* H);
	void shape_gradient2(const FESolidElement& el, vec3d* X, double r, double s, double t, mat3d* H);

protected:
	FEInternalSurface2O	m_surf;
	bool	m_binitJ0;	//!< flag indicating J0 has been initialized
};
