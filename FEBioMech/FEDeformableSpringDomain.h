/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FECore/FEDiscreteDomain.h>
#include <FECore/FEDofList.h>
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

	//! get the total dofs
	const FEDofList& GetDOFList() const override;

public: // overridden from FEElasticDomain
	//! build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	//! calculate stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;
	void MassMatrix(FELinearSystem& LS, double scale) override {}
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override {}

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

protected:
	FEDofList	m_dofU, m_dofR, m_dof;

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

	//! get the total dofs
	const FEDofList& GetDOFList() const override;

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
	void StiffnessMatrix(FELinearSystem& LS) override;
	void MassMatrix(FELinearSystem& LS, double scale) override {}
	void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override {}

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

protected:
	FEDofList	m_dofU, m_dofR, m_dof;
};
