#pragma once
#include "FESurface.h"
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
//! This class describes a solute flux on a surface element
//! bc = 0 for solvent flux

class FESoluteFlux : public FEBoundaryCondition
	{
	public:
		FESoluteFlux() { s[0] = s[1] = s[2] = s[3] = 1.0; bc = 0; blinear = false; }
		
	public:
		double	s[4];		// nodal scale factors
		int		face;		// face number
		int		lc;			// load curve
		int		bc;			// degree of freedom
		bool	blinear;	// linear or not (true is non-follower, false is follower)
	};

//-----------------------------------------------------------------------------
//! The flux surface is a surface domain that sustains a solute flux boundary
//! condition
//!
class FESoluteFluxSurface : public FESurface
	{
	public:
		//! constructor
		FESoluteFluxSurface(FEMesh* pm) : FESurface(pm) {}
		
		//! allocate storage
		void create(int n)
		{
			FESurface::create(n);
			m_PC.resize(n);
		}
		
		//! clone
		FEDomain* Clone()
		{
			FESoluteFluxSurface* ps = new FESoluteFluxSurface(m_pMesh);
			ps->m_PC = m_PC;
			return ps;
		}
		
		//! get a flux BC
		FESoluteFlux& SoluteFlux(int n) { return m_PC[n]; }
		
		//! calculate flux stiffness
		void StiffnessMatrix(FESolidSolver* psolver);
		
		//! calculate residual
		void Residual(FESolidSolver* psolver, vector<double>& R);
		
		//! serialize data
		void Serialize(FEM& fem, DumpFile& ar);
		
	protected:
		//! calculate stiffness for an element
		void FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& vn);
		
		//! Calculates volumetric flow rate due to flux
		bool FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn);
		
		//! Calculates the linear volumetric flow rate due to flux (ie. non-follower)
		bool LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& vn);
		
	protected:
		// solute flux boundary data
		vector<FESoluteFlux>	m_PC;		//!< solute flux boundary cards
	};
