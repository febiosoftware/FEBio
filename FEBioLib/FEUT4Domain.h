#pragma once
#include "FECore/FENodeElemList.h"
#include "FECore/tens4d.h"
#include "FECore/FEElasticMaterial.h"
#include "FEElasticSolidDomain.h"

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
		mat3d	Fi;		// average deformation gradient
		mat3ds	si;		// nodal stress
	};

public:
	//! constructor
	FEUT4Domain(FEMesh* pm, FEMaterial* pmat);

	//! destructor
	~FEUT4Domain();

	//! clone function
	FEDomain* Clone();

	//! data serialization
	void Serialize(DumpFile& ar);

	//! return the node-element list for this domain
	FENodeElemList& GetNodeElemList() { return m_NEL; }

	//! initialization function
	bool Initialize(FEModel& fem);

public: // overrides from FEElasticDomain

	//! Update stresses
	void UpdateStresses(FEModel& fem);

	//! calculates the internal force vector
	void InternalForces(FENLSolver* psolver, vector<double>& R);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FENLSolver* psolver);

protected:
	//! calculates the nodal internal forces
	void NodalInternalForces(FENLSolver* psolver, vector<double>& R);

	//! calculates the element internal forces
	void ElementInternalForces(FENLSolver* psolver, vector<double>& R);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FESolidElement& el, vector<double>& fe);

protected:
	//! Calculates the element stiffness matrix
	void ElementalStiffnessMatrix(FENLSolver* psolver);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(FEModel& fem, int iel, matrix& ke);

	//! Calculates the nodal stiffness matrix
	void NodalStiffnessMatrix(FENLSolver* psolver);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void ElementMaterialStiffness(FEModel& fem, FESolidElement& el, matrix& ke);

	//! nodal geometry stiffness contribution
	void NodalGeometryStiffness(UT4NODE& node, matrix& ke);

	//! nodal material stiffness contribution
	void NodalMaterialStiffness(UT4NODE& node, matrix& ke, FESolidMaterial* pme);

protected:
	//! calculate the volume of a tet element
	double TetVolume(vec3d* r);

	tens4ds Cvol(const tens4ds& C, const mat3ds& S);

public:
	static double	m_alpha;	//!< stabilization factor alpha
	static bool		m_bdev;		//!< use deviatoric components only for the element contribution

private:
	vector<int>		m_tag;	//!< nodal tags
	vector<UT4NODE>	m_Node;	//!< Nodal data
	vector<double>	m_Ve0;	//!< initial element volumes

	double	(*m_Be)[6][3];
	double	(*m_DB)[6][3];
	double	(*m_Ge)[4][3];

	FENodeElemList	m_NEL;
};
