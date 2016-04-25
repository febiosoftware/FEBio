#pragma once
#include "FEElasticSolidDomain.h"
#include "FECore/tens3d.h"
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
//! This class implements a domain used in an elastic remodeling problem.
//! It differs from the FEElasticSolidDomain in that it adds a stiffness matrix
//! due to the deformation dependent density.
class FEElasticMultiscaleDomain2O : public FEElasticSolidDomain
{
	// Helper class for evaluating the discrete-Galerkin contribution
	// It stores the data needed for evaluating the integrals over the
	// internal surface.
	class FEInternalSurface2O 
	{
	public:
		struct Data
		{
			FEMaterialPoint*	m_pt[2];
			vec3d	ksi[2];
		};

	public:
		FEInternalSurface2O ();

		// initialize the data structure
		bool Initialize(FEElasticMultiscaleDomain2O* dom);

		int Elements() const { return m_ps->Elements(); }

		FESurfaceElement& Element(int i) { return m_ps->Element(i); }

		Data& GetData(int i) { return m_data[i]; }

	private:
		FESurface*		m_ps;
		vector<Data>	m_data;
	};

public:
	//! constructor
	FEElasticMultiscaleDomain2O(FEModel* pfem);
	
	//! initialize class
	bool Initialize(FEModel& fem);

	//! initialize elements
	void InitElements();

	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	void UpdateElementStress(int iel, double dt);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! overridden from FEElasticSolidDomain
	void Update();

protected:
	void InternalWorkFlux(FEGlobalVector& R);
	void InternalElementWorkFlux(FESurfaceElement& el, vector<double>& fe, int& nd);

	void ElementInternalForce_PF(FESolidElement& el, vector<double>& fe);
	void ElementInternalForce_QG(FESolidElement& el, vector<double>& fe);

private:
	void UpdateInternalSurfaceStresses();

public:
	// --- S T I F F N E S S ---
	//! calculates the solid element stiffness matrix
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);
	void ElementMaterialStiffness(FESolidElement &el, matrix &ke);

	void defhess(FESolidElement &el, int n, tens3drs &G);
	void defhess(FESolidElement &el, double r, double s, double t, tens3drs &G);

private:
	FEInternalSurface2O	m_surf;
};
