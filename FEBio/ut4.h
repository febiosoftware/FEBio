#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Domain for nodally integrated tet elements 

//! This class implements the uniform nodal strain tetrahedron with
//! isochoric stabilization as described by Gee, Dohrmann, Key and Wall (2009)
//!
class FEUT4Domain : public FEElasticSolidDomain
{
	struct UT4NODE
	{
		int		inode;	// index of FE node
		double	Vi;		// reference nodal volume
		double	vi;		// current nodal volume
	};

public:
	//! constructor
	FEUT4Domain(FEMesh* pm);

	//! clone function
	FEDomain* Clone()
	{
		FEUT4Domain* pd = new FEUT4Domain(m_pMesh);
		pd->m_Elem = m_Elem; pd->m_pMesh = m_pMesh;
		return pd;
	}

	//! initialization function
	bool Initialize(FEM& fem);

	//! calculates the total residual
	void Residual(FESolidSolver* psolver, vector<double>& R);

	//! Update stresses
	void UpdateStresses(FEM& fem);

protected:
	//! calculates the nodal residual
	void NodalResidual(FESolidSolver* psolver, vector<double>& R);

	//! calculates the element residual
	void ElementResidual(FESolidSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FESolidElement& el, vector<double>& fe);

protected:
	//! calculate the volume of a tet element
	double TetVolume(vec3d* r);

private:
	double	m_alpha;	//!< stabilization factor alpha

	vector<int>		m_tag;	//!< nodal tags
	vector<UT4NODE>	m_Node;	//!< Nodal data
};
