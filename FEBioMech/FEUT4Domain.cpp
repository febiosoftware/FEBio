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



#include "stdafx.h"
#include "FEUT4Domain.h"
#include "FECore/FEMesh.h"
#include "FECore/FEModel.h"
#include "FECore/FEGlobalMatrix.h"

//-----------------------------------------------------------------------------
// This function converts the Cauchy stress to a 2nd-PK stress
mat3ds cauchy_to_pk2(mat3ds& s, mat3d& F)
{
	mat3d Fi = F.inverse();
	mat3d Fit = Fi.transpose();
	double J = F.det();

	mat3d S = (Fi*s*Fit)*J;
	return S.sym();
}

//-----------------------------------------------------------------------------
// This function converts the 2nd-PK stress to a Cauchy stress
mat3ds pk2_to_cauchy(mat3ds& S, mat3d& F)
{
	mat3d Ft = F.transpose();
	double Ji = 1/F.det();

	mat3d s = (F*S*Ft)*Ji;
	return s.sym();
}

//-----------------------------------------------------------------------------
// This function converts the spatial tangent to the material tangent
tens4ds spatial_to_material(tens4ds& c, mat3d& F)
{
	mat3d Fi = F.inverse();
	double J = F.det();

	double C[6][6] = {0};

	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
				{
					C[0][0] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[0][l]*c(i,j,k,l);
					C[0][1] += J*Fi[0][i]*Fi[0][j]*Fi[1][k]*Fi[1][l]*c(i,j,k,l);
					C[0][2] += J*Fi[0][i]*Fi[0][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[0][3] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[0][4] += J*Fi[0][i]*Fi[0][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[0][5] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[1][1] += J*Fi[1][i]*Fi[1][j]*Fi[1][k]*Fi[1][l]*c(i,j,k,l);
					C[1][2] += J*Fi[1][i]*Fi[1][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[1][3] += J*Fi[1][i]*Fi[1][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[1][4] += J*Fi[1][i]*Fi[1][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[1][5] += J*Fi[1][i]*Fi[1][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[2][2] += J*Fi[2][i]*Fi[2][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[2][3] += J*Fi[2][i]*Fi[2][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[2][4] += J*Fi[2][i]*Fi[2][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[2][5] += J*Fi[2][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[3][3] += J*Fi[0][i]*Fi[1][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[3][4] += J*Fi[0][i]*Fi[1][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[3][5] += J*Fi[0][i]*Fi[1][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[4][4] += J*Fi[1][i]*Fi[2][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[4][5] += J*Fi[1][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[5][5] += J*Fi[0][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);
				}

	return tens4ds(C);
}

//-----------------------------------------------------------------------------
// This function converts the material tangent to the spatial tangent
tens4ds material_to_spatial(tens4ds& C, mat3d& F)
{
	double Ji = 1/F.det();

	double c[6][6] = {0};

	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
				{
					c[0][0] += Ji*F[0][i]*F[0][j]*F[0][k]*F[0][l]*C(i,j,k,l);
					c[0][1] += Ji*F[0][i]*F[0][j]*F[1][k]*F[1][l]*C(i,j,k,l);
					c[0][2] += Ji*F[0][i]*F[0][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[0][3] += Ji*F[0][i]*F[0][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[0][4] += Ji*F[0][i]*F[0][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[0][5] += Ji*F[0][i]*F[0][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[1][1] += Ji*F[1][i]*F[1][j]*F[1][k]*F[1][l]*C(i,j,k,l);
					c[1][2] += Ji*F[1][i]*F[1][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[1][3] += Ji*F[1][i]*F[1][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[1][4] += Ji*F[1][i]*F[1][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[1][5] += Ji*F[1][i]*F[1][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[2][2] += Ji*F[2][i]*F[2][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[2][3] += Ji*F[2][i]*F[2][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[2][4] += Ji*F[2][i]*F[2][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[2][5] += Ji*F[2][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[3][3] += Ji*F[0][i]*F[1][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[3][4] += Ji*F[0][i]*F[1][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[3][5] += Ji*F[0][i]*F[1][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[4][4] += Ji*F[1][i]*F[2][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[4][5] += Ji*F[1][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[5][5] += Ji*F[0][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);
				}

	return tens4ds(c);
}

//-----------------------------------------------------------------------------
void FEUT4Domain::UT4NODE::Serialize(DumpStream& ar)
{
	ar & inode;
	ar & Vi & vi;
	ar & Fi & si;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEUT4Domain, FEElasticSolidDomain)
	ADD_PARAMETER(m_alpha, "alpha");
	ADD_PARAMETER(m_bdev, "iso_stab");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor for the UT4Domain
FEUT4Domain::FEUT4Domain(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
	m_alpha = 0.05;
	m_bdev = false;

	m_Be = 0;
	m_DB = 0;
	m_Ge = 0;
}

//-----------------------------------------------------------------------------
void FEUT4Domain::Serialize(DumpStream& ar)
{
	FEElasticSolidDomain::Serialize(ar);
	ar & m_alpha & m_bdev;
	ar & m_tag;
	ar & m_NODE;
	ar & m_Ve0;

	if (ar.IsLoading())
	{
		// create the node-element list
		m_NEL.Create(*this);
		
		// find the largest valence
		int Nmax = m_NEL.MaxValence();

		// allocate buffers
		m_Ge = new double[Nmax*4][4][3];
		m_Be = new double[Nmax*4][6][3];
		m_DB = new double[Nmax*4][6][3];
	}
}

//-----------------------------------------------------------------------------
FEUT4Domain::~FEUT4Domain()
{
	if (m_DB) delete [] m_DB;
	if (m_Be) delete [] m_Be;
	if (m_Ge) delete [] m_Ge;
}

//-----------------------------------------------------------------------------
void FEUT4Domain::BuildMatrixProfile(FEGlobalMatrix& M)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    DOFS& fedofs = fem.GetDOFS();
    int MAX_NDOFS = fedofs.GetTotalDOFS();

	// we'll need the node-element list
	FENodeElemList& NEL = GetNodeElemList();
	assert(NEL.Size() > 0);

	vector<int> LM, elm;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		int NE = NEL.Valence(i);
		if (NE > 0)
		{
			LM.assign(NE*4*MAX_NDOFS, -1);
			FEElement** ppe = NEL.ElementList(i);
			for (int n=0; n<NE; ++n)
			{
				UnpackLM(*ppe[n], elm);
				for (int j=0; j<(int)elm.size(); ++j) LM[n*4*MAX_NDOFS + j] = elm[j];
			}
			M.build_add(LM);
		}
	}
}

//-----------------------------------------------------------------------------
//! Set UT4 parameters
void FEUT4Domain::SetUT4Parameters(double alpha, bool bdev)
{
	m_alpha = alpha;
	m_bdev = bdev;
}

//-----------------------------------------------------------------------------
//! override Create so we can grab the ut4 parameters
bool FEUT4Domain::Create(int nelems, FE_Element_Spec espec)
{
	// NOTE: Commented this out, since during restart the espec is a dummy variable. 
	// some sanity checks first
//	if (espec.eclass != FE_Element_Class::FE_ELEM_SOLID) return false;
//	if (espec.eshape != FE_Element_Shape::ET_TET4) return false;

	m_alpha = espec.m_ut4_alpha;
	m_bdev = espec.m_ut4_bdev;
	espec.m_but4 = true;

	return FEElasticSolidDomain::Create(nelems, espec);
}

//-----------------------------------------------------------------------------
//! Initialization for the UT4Domain.
//! Note that we first initialize the base class before initializing the domain.
bool FEUT4Domain::Init()
{
	// first call the base class
	if (FEElasticSolidDomain::Init() == false) return false;

	// next, we need to identify all the nodes that belong to this domain
	// we do this by looping over all the elements and tagging the nodes
	int NN = m_pMesh->Nodes();
	m_tag.assign(NN, -1);

	int i, j;
	int N = (int)m_Node.size();
	for (i=0; i<N; ++i) m_tag[m_Node[i]] = i;

	// allocate node structure
	m_NODE.resize(N);

	// initialize nodes
	for (i=0; i<N; ++i)
	{
		UT4NODE& n = m_NODE[i];
		n.inode = m_Node[i];
		n.Vi = n.vi = 0.0;
		n.si.zero();
	}

	// calculate the reference nodal volume
	// we do this here since this volume never changes
	int NE = Elements();
	m_Ve0.resize(NE);
	double Ve;
	vec3d r0[4];
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];

		// get the current element coordinates
		for (j=0; j<4; ++j) r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;

		// calculate the initial volume
		m_Ve0[i] = Ve = TetVolume(r0);

		// now assign one-quart to each node
		for (int j=0; j<4; ++j) m_NODE[ m_tag[el.m_node[j]]].Vi += 0.25*Ve;
	}

	// create the node-element list
	m_NEL.Create(*this);

	// find the largest valence
	int Nmax = m_NEL.MaxValence();

	// allocate buffers
	m_Ge = new double[Nmax*4][4][3];
	m_Be = new double[Nmax*4][6][3];
	m_DB = new double[Nmax*4][6][3];

	return true;
}

//-----------------------------------------------------------------------------
//! Update the nodal and element stresses
void FEUT4Domain::Update(const FETimeInfo& tp)
{
	// updating the element stresses is easy, since we only
	// need to call the base class
	FEElasticSolidDomain::Update(tp);

	// next we update the nodal data
	int i, j, NE = Elements();
	for (i=0; i<(int) m_NODE.size(); ++i) { m_NODE[i].vi = 0; m_NODE[i].Fi.zero(); }
	double Ve, ve;
	mat3d Fe;
	vec3d rt[4];
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];

		// nodal coordinates
		for (j=0; j<4; ++j) rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;

		// calculate the volume
		Ve = m_Ve0[i];
		ve = TetVolume(rt);

		// calculate the deformation gradient
		defgrad(el, Fe, 0);

		// now assign one-quart to each node
		for (j=0; j<4; ++j) 
		{
			UT4NODE& n = m_NODE[ m_tag[el.m_node[j]]];
			n.vi += 0.25*ve;
			n.Fi += Fe*(0.25*Ve / n.Vi);
		}
	}

	// create a material point
	// TODO: this will set the Q variable to a unit-matrix
	//		 in other words, we loose the material axis orientation
	//		 For now, I solve this by copying the Q parameter
	//       from the first element that the node connects to
	FEElasticMaterialPoint ep;
	FEMaterialPoint mp(&ep);
	mp.Init();

	// loop over all the nodes
	for (i=0; i<(int) m_NODE.size(); ++i)
	{
		UT4NODE& node = m_NODE[i];

		// set the material point data
		mp.m_r0 = m_pMesh->Node(node.inode).m_r0;
		mp.m_rt = m_pMesh->Node(node.inode).m_rt;

		ep.m_F = node.Fi;
		ep.m_J = ep.m_F.det();

		// calculate the stress
		node.si = m_pMat->Stress(mp);
	}
}

//-----------------------------------------------------------------------------
//! Calculate volume for tets using the determinant rule
double FEUT4Domain::TetVolume(vec3d* r)
{
	double A[4][4];
	A[0][0] = r[0].x; A[0][1] = r[1].x; A[0][2] = r[2].x; A[0][3] = r[3].x;
	A[1][0] = r[0].y; A[1][1] = r[1].y; A[1][2] = r[2].y; A[1][3] = r[3].y;
	A[2][0] = r[0].z; A[2][1] = r[1].z; A[2][2] = r[2].z; A[2][3] = r[3].z;
	A[3][0] =    1.0; A[3][1] =    1.0; A[3][2] =    1.0; A[3][3] =    1.0;

	double V = 0;
	V += mat3d(A[0][1], A[0][2], A[0][3], A[1][1], A[1][2], A[1][3], A[2][1], A[2][2], A[2][3]).det();
	V -= mat3d(A[0][0], A[0][2], A[0][3], A[1][0], A[1][2], A[1][3], A[2][0], A[2][2], A[2][3]).det();
	V += mat3d(A[0][0], A[0][1], A[0][3], A[1][0], A[1][1], A[1][3], A[2][0], A[2][1], A[2][3]).det();
	V -= mat3d(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]).det();

	return V/6.0;
}

//-----------------------------------------------------------------------------
//! The residual is defined as the sum of the nodal residual and the
//! element residual
void FEUT4Domain::InternalForces(FEGlobalVector& R)
{
	// Calculate the nodal contribution
	NodalInternalForces(R);

	// Calculate the element contribution
	ElementInternalForces(R);
}

//-----------------------------------------------------------------------------
//! This function calculates the nodal contribution to the residual
void FEUT4Domain::NodalInternalForces(FEGlobalVector& R)
{
	// inverse jacobian matrix
	double Ji[3][3];

	// shape function derivatives
	double Gx[4], Gy[4], Gz[4];

	// B-matrix
	double Be[6][3] = {0};

	int i, j, n;
	double Ve;
	vector<int> elm;
	vector<int> LM(3);
	vector<int> en(1);
	vector<double> fe(3);

	// loop over all the nodes
	int NN = (int) m_NODE.size();
	for (i=0; i<NN; ++i)
	{
		UT4NODE& node = m_NODE[i];

		// get the elements that belong to this list
		int NE = m_NEL.Valence(node.inode);
		FEElement** ppel = m_NEL.ElementList(node.inode);
		int* peli = m_NEL.ElementIndexList(node.inode);

		// get the nodal deformation gradient
		mat3d FI = node.Fi;

		mat3ds S;
		if (m_bdev)
		{
			S = node.si - node.si.dev()*m_alpha;
		}
		else
		{
			S = node.si*(1 - m_alpha);
		}

		// convert the Cauchy stress to a PK2-stress
		S = cauchy_to_pk2(S, FI);

		// loop over all elements that belong to this node
		for (n=0; n<NE; ++n)
		{
			FESolidElement& el = dynamic_cast<FESolidElement&>(*ppel[n]);
			UnpackLM(el, elm);

			// calculate element volume
			// TODO: we should store this somewhere instead of recalculating it
			Ve = m_Ve0[peli[n]];

			// The '-' sign is so that the internal forces get subtracted from the residual (R = Fext - Fint)
			double w = -0.25* Ve;

			// calculate the jacobian
			invjac0(el, Ji, 0);

			// get the shape function derivatives
			const double* Gr, *Gs, *Gt;
			Gr = el.Gr(0);
			Gs = el.Gs(0);
			Gt = el.Gt(0);

			for (j=0; j<4; ++j)
			{
				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !
				Gx[j] = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
				Gy[j] = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
				Gz[j] = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];
			}

			// get the element deformation gradient
			mat3d Fe = FI;

			// loop over element nodes
			for (j=0; j<4; ++j)
			{
				// setup nonlinear element B-matrix
				Be[0][0] = Fe[0][0]*Gx[j]; Be[0][1] = Fe[1][0]*Gx[j]; Be[0][2] = Fe[2][0]*Gx[j];
				Be[1][0] = Fe[0][1]*Gy[j]; Be[1][1] = Fe[1][1]*Gy[j]; Be[1][2] = Fe[2][1]*Gy[j];
				Be[2][0] = Fe[0][2]*Gz[j]; Be[2][1] = Fe[1][2]*Gz[j]; Be[2][2] = Fe[2][2]*Gz[j];
				Be[3][0] = Fe[0][0]*Gy[j] + Fe[0][1]*Gx[j]; Be[3][1] = Fe[1][0]*Gy[j] + Fe[1][1]*Gx[j]; Be[3][2] = Fe[2][0]*Gy[j] + Fe[2][1]*Gx[j];
				Be[4][0] = Fe[0][1]*Gz[j] + Fe[0][2]*Gy[j]; Be[4][1] = Fe[1][1]*Gz[j] + Fe[1][2]*Gy[j]; Be[4][2] = Fe[2][1]*Gz[j] + Fe[2][2]*Gy[j];
				Be[5][0] = Fe[0][2]*Gx[j] + Fe[0][0]*Gz[j]; Be[5][1] = Fe[1][2]*Gx[j] + Fe[1][0]*Gz[j]; Be[5][2] = Fe[2][2]*Gx[j] + Fe[2][0]*Gz[j];

				LM[0] = elm[3*j  ];
				LM[1] = elm[3*j+1];
				LM[2] = elm[3*j+2];

				en[0] = el.m_node[j];
	
				// calculate nodal force
				fe[0] = w*(Be[0][0]*S.xx() + Be[1][0]*S.yy() + Be[2][0]*S.zz() + Be[3][0]*S.xy() + Be[4][0]*S.yz() + Be[5][0]*S.xz());
				fe[1] = w*(Be[0][1]*S.xx() + Be[1][1]*S.yy() + Be[2][1]*S.zz() + Be[3][1]*S.xy() + Be[4][1]*S.yz() + Be[5][1]*S.xz());
				fe[2] = w*(Be[0][2]*S.xx() + Be[1][2]*S.yy() + Be[2][2]*S.zz() + Be[3][2]*S.xy() + Be[4][2]*S.yz() + Be[5][2]*S.xz());

				// assemble in global residual
				R.Assemble(en, LM, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the element contribution to the residual
void FEUT4Domain::ElementInternalForces(FEGlobalVector& R)
{
	// element force vector
	vector<double> fe;

	vector<int> lm;

	int NE = (int)m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		ElementInternalForces(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEUT4Domain::ElementInternalForces(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;

	const double* Gr, *Gs, *Gt;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		detJt = invjact(el, Ji, n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds s = pt.m_s;

		// take the deviatoric component and multiply it
		// with the stabilization factor
		if (m_bdev)
		{
			s = s.dev()*m_alpha;
		}
		else
		{
			s = s*m_alpha;
		}

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix. The stiffness matrix is a sum of the
//! nodal and element stiffness matrices
void FEUT4Domain::StiffnessMatrix(FELinearSystem& LS)
{
	// calculate nodal stiffness matrix
	NodalStiffnessMatrix(LS);

	// calculate element stiffness matrix
	ElementalStiffnessMatrix(LS);
}

//-----------------------------------------------------------------------------
//! Calculates the nodal contribution to the global stiffness matrix
void FEUT4Domain::NodalStiffnessMatrix(FELinearSystem& LS)
{
	vector<int> elm;
	vector<int> LM;
	vector<int> en;

	// loop over all the nodes
	int NN = (int) m_NODE.size(), ni, nj;
	for (int i=0; i<NN; ++i)
	{
		// get the next node
		UT4NODE& node = m_NODE[i];
		FEElement** ppe = m_NEL.ElementList(node.inode);

		// get its valence
		int NE = m_NEL.Valence(node.inode);

		// allocate an element stiffness matrix
		FEElementMatrix ke(NE*4*3, NE*4*3); ke.zero();

		// calculate the geometry stiffness for this node
		NodalGeometryStiffness(node, ke);

		// calculate the material stiffness for this node
		NodalMaterialStiffness(node, ke, m_pMat);

		// it is assumed that the previous function only build the upper-triangular part
		// so now we build the lower-triangular by copying it from the upper-triangular part
		for (ni=0; ni<NE*12; ++ni)
		{
			for (nj=0; nj<ni; ++nj)
			{
				ke[ni][nj] = ke[nj][ni];
			}
		}

		// create the LM and the en array
		LM.resize(NE*4*3);
		en.resize(NE*4  );
		for (ni=0; ni<NE; ++ni)
		{
			FEElement& el = *ppe[ni];
			UnpackLM(el, elm);
			for (int i=0; i<4; ++i)
			{
				LM[ni*4*3+3*i  ] = elm[3*i  ];
				LM[ni*4*3+3*i+1] = elm[3*i+1];
				LM[ni*4*3+3*i+2] = elm[3*i+2];

				en[ni*4+i] = el.m_node[i];
			}
		}
	
		// assemble the stiffness
		ke.SetNodes(en);
		ke.SetIndices(LM);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates the nodal geometry stiffness contribution
void FEUT4Domain::NodalGeometryStiffness(UT4NODE& node, matrix& ke)
{
	int i, j, ni, nj;

	// get the element list 
	int NE = m_NEL.Valence(node.inode);
	FEElement** ppe = m_NEL.ElementList(node.inode);
	int* peli = m_NEL.ElementIndexList(node.inode);

	// get the nodal stress
	mat3ds S;
	if (m_bdev)
	{
		S = node.si - node.si.dev()*m_alpha;
	}
	else
	{
		S = node.si*(1 - m_alpha);
	}

	// convert Cauchy stress to PK2 stress
	S = cauchy_to_pk2(S, node.Fi);

	// calculate the gradN matrices
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& ei = dynamic_cast<FESolidElement&>(*ppe[ni]);

		// calculate the jacobian
		double Ji[3][3];
		invjac0(ei, Ji, 0);

		// get the shape function derivatives
		const double* Gr, *Gs, *Gt;
		Gr = ei.Gr(0);
		Gs = ei.Gs(0);
		Gt = ei.Gt(0);

		for (j=0; j<4; ++j)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			m_Ge[ni][j][0] = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
			m_Ge[ni][j][1] = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
			m_Ge[ni][j][2] = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];
		}
	}

	// stiffness matrices
	double kij;

	// loop over all the elements
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& ei = dynamic_cast<FESolidElement&>(*ppe[ni]);

		// calculate element volume
		double Vi = m_Ve0[peli[ni]];
		double wi = 0.25* Vi*node.Vi/ node.Vi;

		// loop over the elements again
		for (nj=ni; nj<NE; ++nj)
		{
			FESolidElement& ej = dynamic_cast<FESolidElement&>(*ppe[nj]);

			// calculate element volume
			double Vj = m_Ve0[peli[nj]];
			double wj = 0.25* Vj / node.Vi;

			// We're ready to rock and roll!
			double sg[3];
			for (i=0; i<4; ++i)
			{
				double (&Gi)[3] = *(m_Ge[ni] + i);
				int mi = ni*12+i*3;
				int j0 = (ni==nj?i:0);
				for (j=j0; j<4; ++j)
				{
					double (&Gj)[3] = *(m_Ge[nj] + j);
					int mj = nj*12+j*3;

					sg[0] = S.xx()*Gj[0] + S.xy()*Gj[1] + S.xz()*Gj[2];
					sg[1] = S.xy()*Gj[0] + S.yy()*Gj[1] + S.yz()*Gj[2];
					sg[2] = S.xz()*Gj[0] + S.yz()*Gj[1] + S.zz()*Gj[2];

					kij = wi*wj*(Gi[0]*sg[0]+Gi[1]*sg[1]+Gi[2]*sg[2]);

					// copy to element stiffness matrix
					ke[mi  ][mj  ] += kij;
					ke[mi+1][mj+1] += kij;
					ke[mi+2][mj+2] += kij;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the volumetric contribution of the spatial tangent stiffness
//!
tens4ds FEUT4Domain::Cvol(const tens4ds& C, const mat3ds& S)
{
	double p = -S.tr()/3.0;

	mat3dd I(1);	// Identity
	tens4ds I4  = dyad4s(I);

	// Note the slightly different form than in the paper.
	// This is because Cvol needs to have the proper symmetries
	return I4*(2*p) + dyad1s(I, S)/3.0 + dyad1s(I, C.dot(I))/6.0;
}

//-----------------------------------------------------------------------------
//! Calculates the nodal material stiffness contribution
void FEUT4Domain::NodalMaterialStiffness(UT4NODE& node, matrix& ke, FESolidMaterial* pme)
{
	// get the number of elements this nodes connects
	int NE = m_NEL.Valence(node.inode);
	FEElement** ppe = m_NEL.ElementList(node.inode);
	int* peli = m_NEL.ElementIndexList(node.inode);

	// create a material point
	FEElasticMaterialPoint ep;
	FEMaterialPoint mp(&ep);
	mp.Init();

	// set the material point data
	mp.m_r0 = m_pMesh->Node(node.inode).m_r0;
	mp.m_rt = m_pMesh->Node(node.inode).m_rt;

	ep.m_F = node.Fi;
	ep.m_J = ep.m_F.det();

	// set the Cauchy-stress
	ep.m_s = node.si;

	// Calculate the spatial tangent
	tens4ds C = pme->Tangent(mp);

	// Next, we need to subtract the volumetric contribution Cvol
	if (m_bdev)
	{
		// subtract the isochoric component from C;
		// C = C - a*Ciso = C - (a*(C - Cvol)) = (1-a)*C + a*Cvol
		C = C*(1 - m_alpha) + Cvol(C, ep.m_s)*m_alpha;
	}
	else
	{
		C = C*(1 - m_alpha);
	}

	// convert spatial matrix to material
	C = spatial_to_material(C, ep.m_F);

	// extract the 'D' matrix
	double D[6][6] = {0};
	C.extract(D);

	// loop over all the elements
	int i, j, ni, nj;

	// calculate B-matrices
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& el = dynamic_cast<FESolidElement&>(*ppe[ni]);

		// calculate the jacobian
		double Ji[3][3];
		invjac0(el, Ji, 0);

		// get the shape function derivatives
		const double* Gr, *Gs, *Gt;
		Gr = el.Gr(0);
		Gs = el.Gs(0);
		Gt = el.Gt(0);

		// get element deformation gradient
		mat3d Fe = node.Fi;

		double Gx, Gy, Gz;
		for (j=0; j<4; ++j)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
			Gy = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
			Gz = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];

			double (&Bi)[6][3] = *(m_Be+(4*ni+j));
			Bi[0][0] = Fe[0][0]*Gx; Bi[0][1] = Fe[1][0]*Gx; Bi[0][2] = Fe[2][0]*Gx;
			Bi[1][0] = Fe[0][1]*Gy; Bi[1][1] = Fe[1][1]*Gy; Bi[1][2] = Fe[2][1]*Gy;
			Bi[2][0] = Fe[0][2]*Gz; Bi[2][1] = Fe[1][2]*Gz; Bi[2][2] = Fe[2][2]*Gz;
			Bi[3][0] = Fe[0][0]*Gy + Fe[0][1]*Gx; Bi[3][1] = Fe[1][0]*Gy + Fe[1][1]*Gx; Bi[3][2] = Fe[2][0]*Gy + Fe[2][1]*Gx;
			Bi[4][0] = Fe[0][1]*Gz + Fe[0][2]*Gy; Bi[4][1] = Fe[1][1]*Gz + Fe[1][2]*Gy; Bi[4][2] = Fe[2][1]*Gz + Fe[2][2]*Gy;
			Bi[5][0] = Fe[0][2]*Gx + Fe[0][0]*Gz; Bi[5][1] = Fe[1][2]*Gx + Fe[1][0]*Gz; Bi[5][2] = Fe[2][2]*Gx + Fe[2][0]*Gz;

			double (&DBi)[6][3] = *(m_DB+(4*ni+j));
			DBi[0][0] = (D[0][0]*Bi[0][0]+D[0][1]*Bi[1][0]+D[0][2]*Bi[2][0]+D[0][3]*Bi[3][0]+D[0][4]*Bi[4][0]+D[0][5]*Bi[5][0]);
			DBi[0][1] = (D[0][0]*Bi[0][1]+D[0][1]*Bi[1][1]+D[0][2]*Bi[2][1]+D[0][3]*Bi[3][1]+D[0][4]*Bi[4][1]+D[0][5]*Bi[5][1]);
			DBi[0][2] = (D[0][0]*Bi[0][2]+D[0][1]*Bi[1][2]+D[0][2]*Bi[2][2]+D[0][3]*Bi[3][2]+D[0][4]*Bi[4][2]+D[0][5]*Bi[5][2]);

			DBi[1][0] = (D[1][0]*Bi[0][0]+D[1][1]*Bi[1][0]+D[1][2]*Bi[2][0]+D[1][3]*Bi[3][0]+D[1][4]*Bi[4][0]+D[1][5]*Bi[5][0]);
			DBi[1][1] = (D[1][0]*Bi[0][1]+D[1][1]*Bi[1][1]+D[1][2]*Bi[2][1]+D[1][3]*Bi[3][1]+D[1][4]*Bi[4][1]+D[1][5]*Bi[5][1]);
			DBi[1][2] = (D[1][0]*Bi[0][2]+D[1][1]*Bi[1][2]+D[1][2]*Bi[2][2]+D[1][3]*Bi[3][2]+D[1][4]*Bi[4][2]+D[1][5]*Bi[5][2]);

			DBi[2][0] = (D[2][0]*Bi[0][0]+D[2][1]*Bi[1][0]+D[2][2]*Bi[2][0]+D[2][3]*Bi[3][0]+D[2][4]*Bi[4][0]+D[2][5]*Bi[5][0]);
			DBi[2][1] = (D[2][0]*Bi[0][1]+D[2][1]*Bi[1][1]+D[2][2]*Bi[2][1]+D[2][3]*Bi[3][1]+D[2][4]*Bi[4][1]+D[2][5]*Bi[5][1]);
			DBi[2][2] = (D[2][0]*Bi[0][2]+D[2][1]*Bi[1][2]+D[2][2]*Bi[2][2]+D[2][3]*Bi[3][2]+D[2][4]*Bi[4][2]+D[2][5]*Bi[5][2]);

			DBi[3][0] = (D[3][0]*Bi[0][0]+D[3][1]*Bi[1][0]+D[3][2]*Bi[2][0]+D[3][3]*Bi[3][0]+D[3][4]*Bi[4][0]+D[3][5]*Bi[5][0]);
			DBi[3][1] = (D[3][0]*Bi[0][1]+D[3][1]*Bi[1][1]+D[3][2]*Bi[2][1]+D[3][3]*Bi[3][1]+D[3][4]*Bi[4][1]+D[3][5]*Bi[5][1]);
			DBi[3][2] = (D[3][0]*Bi[0][2]+D[3][1]*Bi[1][2]+D[3][2]*Bi[2][2]+D[3][3]*Bi[3][2]+D[3][4]*Bi[4][2]+D[3][5]*Bi[5][2]);

			DBi[4][0] = (D[4][0]*Bi[0][0]+D[4][1]*Bi[1][0]+D[4][2]*Bi[2][0]+D[4][3]*Bi[3][0]+D[4][4]*Bi[4][0]+D[4][5]*Bi[5][0]);
			DBi[4][1] = (D[4][0]*Bi[0][1]+D[4][1]*Bi[1][1]+D[4][2]*Bi[2][1]+D[4][3]*Bi[3][1]+D[4][4]*Bi[4][1]+D[4][5]*Bi[5][1]);
			DBi[4][2] = (D[4][0]*Bi[0][2]+D[4][1]*Bi[1][2]+D[4][2]*Bi[2][2]+D[4][3]*Bi[3][2]+D[4][4]*Bi[4][2]+D[4][5]*Bi[5][2]);

			DBi[5][0] = (D[5][0]*Bi[0][0]+D[5][1]*Bi[1][0]+D[5][2]*Bi[2][0]+D[5][3]*Bi[3][0]+D[5][4]*Bi[4][0]+D[5][5]*Bi[5][0]);
			DBi[5][1] = (D[5][0]*Bi[0][1]+D[5][1]*Bi[1][1]+D[5][2]*Bi[2][1]+D[5][3]*Bi[3][1]+D[5][4]*Bi[4][1]+D[5][5]*Bi[5][1]);
			DBi[5][2] = (D[5][0]*Bi[0][2]+D[5][1]*Bi[1][2]+D[5][2]*Bi[2][2]+D[5][3]*Bi[3][2]+D[5][4]*Bi[4][2]+D[5][5]*Bi[5][2]);
		}
	}

	// loop over the elements again
	// (we only build the upper-triangular (block) matrix
	for (ni=0; ni<NE; ++ni)
	{
		// calculate element volume
		double Vi = m_Ve0[peli[ni]];
		double wi = 0.25* Vi*node.Vi/ node.Vi;

		for (nj=ni; nj<NE; ++nj)
		{
			// calculate element volume
			double Vj = m_Ve0[peli[nj]];
			double wj = 0.25* Vj / node.Vi;

			double wij = wi*wj;

			// We're ready to rock and roll!
			for (i=0; i<4; ++i)
			{
				double (&Bi)[6][3] = *(m_Be+(ni*4 + i));
				int mi = ni*12+i*3;
				int j0 = (nj==ni?i:0);

				for (j=j0; j<4; ++j)
				{
					// calculate the Bi*D*Bj term
					double (&DBj)[6][3] = *(m_DB+(nj*4 + j));
					int mj = nj*12+j*3;

					ke[mi  ][mj  ] += wij*(Bi[0][0]*DBj[0][0]+Bi[1][0]*DBj[1][0]+Bi[2][0]*DBj[2][0]+Bi[3][0]*DBj[3][0]+Bi[4][0]*DBj[4][0]+Bi[5][0]*DBj[5][0]);
					ke[mi  ][mj+1] += wij*(Bi[0][0]*DBj[0][1]+Bi[1][0]*DBj[1][1]+Bi[2][0]*DBj[2][1]+Bi[3][0]*DBj[3][1]+Bi[4][0]*DBj[4][1]+Bi[5][0]*DBj[5][1]);
					ke[mi  ][mj+2] += wij*(Bi[0][0]*DBj[0][2]+Bi[1][0]*DBj[1][2]+Bi[2][0]*DBj[2][2]+Bi[3][0]*DBj[3][2]+Bi[4][0]*DBj[4][2]+Bi[5][0]*DBj[5][2]);

					ke[mi+1][mj  ] += wij*(Bi[0][1]*DBj[0][0]+Bi[1][1]*DBj[1][0]+Bi[2][1]*DBj[2][0]+Bi[3][1]*DBj[3][0]+Bi[4][1]*DBj[4][0]+Bi[5][1]*DBj[5][0]);
					ke[mi+1][mj+1] += wij*(Bi[0][1]*DBj[0][1]+Bi[1][1]*DBj[1][1]+Bi[2][1]*DBj[2][1]+Bi[3][1]*DBj[3][1]+Bi[4][1]*DBj[4][1]+Bi[5][1]*DBj[5][1]);
					ke[mi+1][mj+2] += wij*(Bi[0][1]*DBj[0][2]+Bi[1][1]*DBj[1][2]+Bi[2][1]*DBj[2][2]+Bi[3][1]*DBj[3][2]+Bi[4][1]*DBj[4][2]+Bi[5][1]*DBj[5][2]);

					ke[mi+2][mj  ] += wij*(Bi[0][2]*DBj[0][0]+Bi[1][2]*DBj[1][0]+Bi[2][2]*DBj[2][0]+Bi[3][2]*DBj[3][0]+Bi[4][2]*DBj[4][0]+Bi[5][2]*DBj[5][0]);
					ke[mi+2][mj+1] += wij*(Bi[0][2]*DBj[0][1]+Bi[1][2]*DBj[1][1]+Bi[2][2]*DBj[2][1]+Bi[3][2]*DBj[3][1]+Bi[4][2]*DBj[4][1]+Bi[5][2]*DBj[5][1]);
					ke[mi+2][mj+2] += wij*(Bi[0][2]*DBj[0][2]+Bi[1][2]*DBj[1][2]+Bi[2][2]*DBj[2][2]+Bi[3][2]*DBj[3][2]+Bi[4][2]*DBj[4][2]+Bi[5][2]*DBj[5][2]);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the element contribution to the global stiffness matrix
void FEUT4Domain::ElementalStiffnessMatrix(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// element stiffness matrix
	FEElementMatrix ke;

	vector<int> elm;

	// repeat over all solid elements
	int NE = (int)m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate the element stiffness matrix
		ElementStiffness(tp, iel, ke);

		// get the element equation numbers
		UnpackLM(el, elm);

		// assemble element matrix in global stiffness matrix
		ke.SetNodes(el.m_node);
		ke.SetIndices(elm);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FEUT4Domain::ElementStiffness(const FETimeInfo& tp, int iel, matrix &ke)
{
	FESolidElement& el = Element(iel);

	// TODO: I need to figure out how to deal with incompressible materials
	//       Incompressible materials require an additional dilatational 
	//       stiffness which is calculated differently from the geometric and
	//       material stiffnesses. 
	// TODO: I don't know if the comment above applies anymore

	// calculate material stiffness (i.e. constitutive component)
	ElementMaterialStiffness(el, ke);

	// calculate geometrical stiffness
	ElementGeometricalStiffness(el, ke);

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	int ndof = 3*el.Nodes();
	int i, j;
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n

void FEUT4Domain::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
{
	vec3d G[FEElement::MAX_NODES];

	// nr of nodes
	int neln = el.Nodes();

	// nr of integration points
	int nint = el.GaussPoints();

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// calculate geometrical element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian
		double detJt = ShapeGradient(el, n, G)*gw[n];

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds s = pt.m_s;

		// we work with the deviatoric component only
		if (m_bdev)
		{
			s = s.dev()*m_alpha;
		}
		else
		{
			s = s*m_alpha;
		}

		for (int i = 0; i<neln; ++i)
			for (int j = i; j<neln; ++j)
			{
				kab = (G[i]*( s*G[j]) )*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEUT4Domain::ElementMaterialStiffness(FESolidElement &el, matrix &ke)
{
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];


	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian and shape function gradients
		double detJt = ShapeGradient(el, n, G)*gw[n];

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// Calculate the tangent
		tens4ds C = m_pMat->Tangent(mp);

		// Next, we need to subtract the volumetric contribution Cvol
		if (m_bdev)
		{
			// subtract the volumetric tensor from C;
			C = (C - Cvol(C, pt.m_s))*m_alpha;
		}
		else
		{
			C = C*m_alpha;
		}

		// extract the 'D' matrix
		C.extract(D);

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (int i = 0, i3 = 0; i<neln; ++i, i3 += 3)
		{
			Gxi = G[i].x;
			Gyi = G[i].y;
			Gzi = G[i].z;

			for (int j = i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = G[j].x;
				Gyj = G[j].y;
				Gzj = G[j].z;

				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;

				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;

				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}
