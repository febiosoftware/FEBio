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
#include "FECore/FENodeElemList.h"
#include "FECore/tens4d.h"
#include "FEElasticMaterial.h"
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! Domain for nodally integrated tet elements 

//! This class implements the uniform nodal strain tetrahedron with
//! isochoric stabilization as described by Gee, Dohrmann, Key and Wall (2009)
//!
class FEBIOMECH_API FEUT4Domain : public FEElasticSolidDomain
{
public:
	struct UT4NODE
	{
		int		inode;	// index of FE node
		double	Vi;		// reference nodal volume
		double	vi;		// current nodal volume
		mat3d	Fi;		// average deformation gradient
		mat3ds	si;		// nodal stress

		void Serialize(DumpStream& ar);
	};

public:
	//! constructor
	FEUT4Domain(FEModel* pfem);

	//! destructor
	~FEUT4Domain();

	//! data serialization
	void Serialize(DumpStream& ar) override;

	//! get nodal data
	int UT4Nodes() { return (int) m_NODE.size(); }
	UT4NODE& UT4Node(int i) { return m_NODE[i]; }

	//! return the node-element list for this domain
	FENodeElemList& GetNodeElemList() { return m_NEL; }

	//! initialization function
	bool Init() override;

	//! build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	//! Set UT4 parameters
	void SetUT4Parameters(double alpha, bool bdev);

	//! override Create so we can grab the ut4 parameters
	bool Create(int nelems, FE_Element_Spec espec) override;

public: // overrides from FEElasticDomain

	//! Update domain data
	void Update(const FETimeInfo& tp) override;

	//! calculates the internal force vector
	void InternalForces(FEGlobalVector& R) override;

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FELinearSystem& LS) override;

protected:
	//! calculates the nodal internal forces
	void NodalInternalForces(FEGlobalVector& R);

	//! calculates the element internal forces
	void ElementInternalForces(FEGlobalVector& R);

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForces(FESolidElement& el, vector<double>& fe);

protected:
	//! Calculates the element stiffness matrix
	void ElementalStiffnessMatrix(FELinearSystem& LS);

	//! calculates the solid element stiffness matrix
	void ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke) override;

	//! Calculates the nodal stiffness matrix
	void NodalStiffnessMatrix(FELinearSystem& LS);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(FESolidElement& el, matrix& ke) override;

	//! material stiffness component
	void ElementMaterialStiffness(FESolidElement& el, matrix& ke) override;

	//! nodal geometry stiffness contribution
	void NodalGeometryStiffness(UT4NODE& node, matrix& ke);

	//! nodal material stiffness contribution
	void NodalMaterialStiffness(UT4NODE& node, matrix& ke, FESolidMaterial* pme);

protected:
	//! calculate the volume of a tet element
	double TetVolume(vec3d* r);

	tens4ds Cvol(const tens4ds& C, const mat3ds& S);

	double	m_alpha;	//!< stabilization factor alpha
	bool	m_bdev;		//!< use deviatoric components only for the element contribution

private:
	vector<int>		m_tag;	//!< nodal tags
	vector<UT4NODE>	m_NODE;	//!< Nodal data
	vector<double>	m_Ve0;	//!< initial element volumes

	double	(*m_Be)[6][3];
	double	(*m_DB)[6][3];
	double	(*m_Ge)[4][3];

	FENodeElemList	m_NEL;

	DECLARE_FECORE_CLASS();
};
