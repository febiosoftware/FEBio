#pragma once
#include <FECore/FEDiscreteDomain.h>
#include "FEElasticDomain.h"
#include "FESpringMaterial.h"

//-----------------------------------------------------------------------------
//! domain for deformable springs
class FEDeformableSpringDomain : public FEDiscreteDomain, public FEElasticDomain
{
public:
	//! constructor
	FEDeformableSpringDomain(FEModel* pfem);

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	void Activate() override;

public: // overridden from FEElasticDomain
	//! build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	//! calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver) override;
	void MassMatrix(FESolver* psolver, double scale) override {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) override {}

	//! Calculates inertial forces for dynamic problems | todo implement (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) override { }

	//! update domain data
	void Update(const FETimeInfo& tp) override {}

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}

protected:
	double InitialLength();
	double CurrentLength();

protected:
	FESpringMaterial*	m_pMat;
	double				m_kbend;	// bending stiffness
	double				m_kstab;	// stabilization penalty
	double				m_L0;	//!< initial spring length

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! domain for deformable springs
//! This approach assumes that the nodes are distributed evenly between anchor 
//! points. An anchor is a point that is constrained (e.g. prescribed, or in contact).
class FEDeformableSpringDomain2 : public FEDiscreteDomain, public FEElasticDomain
{
	struct NodeData
	{
		bool	banchor;
	};

public:
	//! constructor
	FEDeformableSpringDomain2(FEModel* pfem);

	//! Unpack LM data
	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() override { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat) override;

	//! initialization
	bool Init() override;

	//! activation
	void Activate() override;

public:
	//! Set the position of a node
	void SetNodePosition(int node, const vec3d& r);

	//! Anchor (or release) a node
	void AnchorNode(int node, bool banchor);

	//! see if a node is anchored
	bool IsAnchored(int node) { return m_nodeData[node].banchor; }

	//! Update the position of all the nodes
	void UpdateNodes();

	//! Get the net force on this node
	vec3d NodalForce(int node);

	//! get net spring force
	double SpringForce();

	//! tangent
	vec3d Tangent(int node);

public: // overridden from FEElasticDomain

	//! calculate stiffness matrix
	void StiffnessMatrix(FESolver* psolver) override;
	void MassMatrix(FESolver* psolver, double scale) override {}
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf) override {}

	//! Calculates inertial forces for dynamic problems | todo implement (removed assert DSR)
	void InertialForces(FEGlobalVector& R, vector<double>& F) override { }

	//! update domain data
	void Update(const FETimeInfo& tp) override;

	//! internal stress forces
	void InternalForces(FEGlobalVector& R) override;

	//! calculate bodyforces (not used since springs are considered mass-less)
	void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override {}

public:
	double InitialLength();
	double CurrentLength();

protected:
	FESpringMaterial*	m_pMat;
	double				m_L0;	//!< initial wire length
	double				m_Lt;	//!< current wire length
	vector<NodeData>	m_nodeData;
};
