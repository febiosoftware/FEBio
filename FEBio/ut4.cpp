#include "stdafx.h"
#include "ut4.h"
#include "FEMesh.h"
#include "FESolidSolver.h"
#include "FEMicroMaterial.h"

// set the default stabilization factor
double FEUT4Domain::m_alpha = 0.05;
bool FEUT4Domain::m_bdev = false;

//-----------------------------------------------------------------------------
//! Constructor for the UT4Domain
FEUT4Domain::FEUT4Domain(FEMesh *pm, FEMaterial* pmat) : FEElasticSolidDomain(pm, pmat)
{
	// overwrite the domain type
	m_ntype = FE_UT4_DOMAIN; 
}

//-----------------------------------------------------------------------------
//! Initialization for the UT4Domain.
//! Note that we first initialize the base class before initializing the domain.
bool FEUT4Domain::Initialize(FEM& fem)
{
	// first call the base class
	if (FEElasticSolidDomain::Initialize(fem) == false) return false;

	// next, we need to identify all the nodes that belong to this domain
	// we do this by looping over all the elements and tagging the nodes
	int NN = m_pMesh->Nodes();
	m_tag.assign(NN, 0);

	int i, NE = Elements();
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		assert((el.Type() == FE_TET) || (el.Type() == FE_TETG1));
		for (int j=0; j<4; ++j) m_tag[el.m_node[j]] = 1;
	}

	// count the tags
	int N = 0;
	for (i=0; i<NN; ++i) if (m_tag[i] == 1) ++N;
	if (N < 4) return false;

	// allocate node structure
	m_Node.resize(N);

	// initialize nodes
	N = 0;
	for (i=0; i<NN; ++i)
		if (m_tag[i] == 1)
		{
			UT4NODE& n = m_Node[N];
			m_tag[i] = N++;
			n.inode = i;
			n.Vi = n.vi = 0.0;
		}
		else m_tag[i] = -1;

	// calculate the reference nodal volume
	// we do this here since this volume never changes
	double Ve;
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		UnpackElement(el);

		// calculate the initial volume
		Ve = TetVolume(el.r0());

		// now assign one-quart to each node
		for (int j=0; j<4; ++j) m_Node[ m_tag[el.m_node[j]]].Vi += 0.25*Ve;
	}

	// create the node-element list
	m_NEL.Create(*this);

	return true;
}

//-----------------------------------------------------------------------------
//! Update the nodal and element stresses
void FEUT4Domain::UpdateStresses(FEM &fem)
{
	// updating the element stresses is easy, since we only
	// need to call the base class
	FEElasticSolidDomain::UpdateStresses(fem);

	// next we update the nodal data
	int i, NE = Elements();
	for (i=0; i<(int) m_Node.size(); ++i) { m_Node[i].vi = 0; m_Node[i].Fi.zero(); }
	double Ve, ve;
	mat3d Fe;
	for (i=0; i<NE; ++i)
	{
		FESolidElement& el = m_Elem[i];
		UnpackElement(el);

		// calculate the volume
		Ve = TetVolume(el.r0());
		ve = TetVolume(el.rt());

		// calculate the deformation gradient
		el.defgrad(Fe, 0);

		// now assign one-quart to each node
		for (int j=0; j<4; ++j) 
		{
			UT4NODE& n = m_Node[ m_tag[el.m_node[j]]];
			n.vi += 0.25*ve;
			n.Fi += Fe*(0.25*Ve / n.Vi);
		}
	}

	// let's calculate the stresses
	// get the elastic material
	FEElasticMaterial* pme = fem.GetElasticMaterial(m_pMat);

	// create a material point
	// TODO: this will set the Q variable to a unit-matrix
	//		 in other words, we loose the material axis orientation
	FEElasticMaterialPoint pt;
	pt.Init(true);

	// loop over all the nodes
	for (i=0; i<(int) m_Node.size(); ++i)
	{
		UT4NODE& node = m_Node[i];

		// set the material point data
		pt.r0 = m_pMesh->Node(node.inode).m_r0;
		pt.rt = m_pMesh->Node(node.inode).m_rt;

		pt.F = node.Fi;
		pt.J = node.vi / node.Vi;

		// if the material is incompressible we need to set some additional values
		// TODO: since for tets the 3-field formulation is pointless, I really need to
		//       change this.
		FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pme);
		if (pmi)
		{
			// calculate volume ratio
			pt.avgJ = pt.J;

			// Calculate pressure. 
			pt.avgp = pmi->Up(pt.J);
		}

		// calculate the stress
		node.si = pme->Stress(pt);
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
void FEUT4Domain::Residual(FESolidSolver *psolver, vector<double> &R)
{
	// Calculate the nodal contribution
	NodalResidual(psolver, R);

	// Calculate the element contribution
	ElementResidual(psolver, R);
}

//-----------------------------------------------------------------------------
//! This function calculates the nodal contribution to the residual
void FEUT4Domain::NodalResidual(FESolidSolver* psolver, vector<double>& R)
{
	// inverse jacobian matrix
	double Ji[3][3];

	// shape function derivatives
	double Gx[4], Gy[4], Gz[4];

	// B-matrix
	double Be[6][3] = {0};

	int i, j, n;
	double Ve;
	vector<int> LM; LM.resize(3);
	vector<int> en; en.resize(1);
	vector<double> fe; fe.resize(3);

	// loop over all the nodes
	for (i=0; i<m_pMesh->Nodes(); ++i)
	{
		// get the elements that belong to this list
		int NE = m_NEL.Valence(i);
		FEElement** ppel = m_NEL.ElementList(i);

		if (NE > 0)
		{
			assert(m_tag[i] >= 0);

			UT4NODE& node = m_Node[ m_tag[i] ];
			assert(node.inode == i);

			mat3ds S;
			if (m_bdev)
			{
				S = node.si - node.si.dev()*m_alpha;
			}
			else
			{
				S = node.si*(1 - m_alpha);
			}

			// loop over all elements that belong to this node
			for (n=0; n<NE; ++n)
			{
				FESolidElement& el = dynamic_cast<FESolidElement&>(*ppel[n]);
				UnpackElement(el);

				// calculate element volume
				// TODO: we should store this somewhere instead of recalculating it
				Ve = TetVolume(el.r0());

				// The '-' sign is so that the internal forces get subtracted from the residual (R = Fext - Fint)
				double w = -0.25* Ve*node.vi / node.Vi;

				// calculate the jacobian
				el.invjact(Ji, 0);

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

				// loop over element nodes
				for (j=0; j<4; ++j)
				{
					// setup the element B-matrix
					Be[0][0] = Gx[j]; Be[1][1] = Gy[j]; Be[2][2] = Gz[j];
					Be[3][0] = Gy[j]; Be[3][1] = Gx[j]; 
					Be[4][1] = Gz[j]; Be[4][2] = Gy[j]; 
					Be[5][0] = Gz[j]; Be[5][2] = Gx[j];

					LM[0] = el.LM()[3*j  ];
					LM[1] = el.LM()[3*j+1];
					LM[2] = el.LM()[3*j+2];

					en[0] = el.m_node[j];
		
					// calculate nodal force
					fe[0] = w*(Be[0][0]*S.xx() + Be[1][0]*S.yy() + Be[2][0]*S.zz() + Be[3][0]*S.xy() + Be[4][0]*S.yz() + Be[5][0]*S.xz());
					fe[1] = w*(Be[0][1]*S.xx() + Be[1][1]*S.yy() + Be[2][1]*S.zz() + Be[3][1]*S.xy() + Be[4][1]*S.yz() + Be[5][1]*S.xz());
					fe[2] = w*(Be[0][2]*S.xx() + Be[1][2]*S.yy() + Be[2][2]*S.zz() + Be[3][2]*S.xy() + Be[4][2]*S.yz() + Be[5][2]*S.xz());

					// assemble in global residual
					psolver->AssembleResidual(en, LM, fe, R);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the element contribution to the residual
void FEUT4Domain::ElementResidual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	// element force vector
	vector<double> fe;

	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

		// unpack the element
		UnpackElement(el);

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		ElementInternalForces(el, fe);

		// apply body forces
		if (fem.UseBodyForces())
		{
			// Note that this function is in the base class
			// I'm not sure, but I don't think I have to overwrite this one
			BodyForces(fem, el, fe);
		}

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
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
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds s = pt.s;

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
void FEUT4Domain::StiffnessMatrix(FESolidSolver *psolver)
{
	// calculate nodal stiffness matrix
	NodalStiffnessMatrix(psolver);

	// calculate element stiffness matrix
	ElementalStiffnessMatrix(psolver);
}

//-----------------------------------------------------------------------------
//! Calculates the nodal contribution to the global stiffness matrix
void FEUT4Domain::NodalStiffnessMatrix(FESolidSolver *psolver)
{
	// loop over all the nodes
	for (int i=0; i<m_pMesh->Nodes(); ++i)
	{
		// see if we need to process this node
		if (m_tag[i] >= 0)
		{
			// calculate the geometry stiffness for this node
			NodalGeometryStiffness(m_Node[m_tag[i]], psolver);

			// calculate the material stiffness for this node
			NodalMaterialStiffness(m_Node[m_tag[i]], psolver);
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the nodal geometry stiffness contribution
void FEUT4Domain::NodalGeometryStiffness(UT4NODE& node, FESolidSolver* psolver)
{
	int i, j, ni, nj;

	// get the element list 
	int NE = m_NEL.Valence(node.inode);
	FEElement** ppe = m_NEL.ElementList(node.inode);

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

	// create the LM and the en array
	vector<int> LM; LM.resize(NE*4*3);
	vector<int> en; en.resize(NE*4  );

	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& el = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(el);
		for (i=0; i<4; ++i)
		{
			LM[ni*4*3+3*i  ] = el.LM()[3*i  ];
			LM[ni*4*3+3*i+1] = el.LM()[3*i+1];
			LM[ni*4*3+3*i+2] = el.LM()[3*i+2];

			en[ni*4+i] = el.m_node[i];
		}
	}

	// calculate the gradN matrices
	double (*Ge)[4][3] = new double [NE][4][3];
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& ei = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(ei);

		// calculate the jacobian
		double Ji[3][3];
		ei.invjact(Ji, 0);

		// get the shape function derivatives
		const double* Gr, *Gs, *Gt;
		Gr = ei.Gr(0);
		Gs = ei.Gs(0);
		Gt = ei.Gt(0);

		for (j=0; j<4; ++j)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Ge[ni][j][0] = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
			Ge[ni][j][1] = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
			Ge[ni][j][2] = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];
		}
	}

	// stiffness matrices
	double kij;
	matrix ke; 
	ke.Create(NE*4*3, NE*4*3); 
	ke.zero();

	// loop over all the elements
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& ei = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(ei);

		// calculate element volume
		// TODO: we should store this somewhere instead of recalculating it
		double Vi = TetVolume(ei.r0());
		double wi = 0.25* Vi*node.vi/ node.Vi;

		// loop over the elements again
		for (nj=0; nj<NE; ++nj)
		{
			FESolidElement& ej = dynamic_cast<FESolidElement&>(*ppe[nj]);
			UnpackElement(ej);

			// calculate element volume
			// TODO: we should store this somewhere instead of recalculating it
			double Vj = TetVolume(ej.r0());
			double wj = 0.25* Vj / node.Vi;

			// We're ready to rock and roll!
			double sg[3];
			for (i=0; i<4; ++i)
			{
				double (&Gi)[3] = *(Ge[ni] + i);
				for (j=0; j<4; ++j)
				{
					double (&Gj)[3] = *(Ge[nj] + j);

					sg[0] = S.xx()*Gj[0] + S.xy()*Gj[1] + S.xz()*Gj[2];
					sg[1] = S.xy()*Gj[0] + S.yy()*Gj[1] + S.yz()*Gj[2];
					sg[2] = S.xz()*Gj[0] + S.yz()*Gj[1] + S.zz()*Gj[2];

					kij = wi*wj*(Gi[0]*sg[0]+Gi[1]*sg[1]+Gi[2]*sg[2]);

					// copy to element stiffness matrix
					ke[ni*12+i*3  ][nj*12+j*3  ] = kij;
					ke[ni*12+i*3+1][nj*12+j*3+1] = kij;
					ke[ni*12+i*3+2][nj*12+j*3+2] = kij;
				}
			}
		}
	}

	// assemble the stiffness
	psolver->AssembleStiffness(en, LM, ke);

	// clean up
	delete [] Ge;
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
void FEUT4Domain::NodalMaterialStiffness(UT4NODE& node, FESolidSolver* psolver)
{
	// get the FE data
	FEM& fem = psolver->m_fem;

	// Get the material for the domain
	FEElasticMaterial* pme = fem.GetElasticMaterial(m_pMat);

	// create a material point
	// TODO: this will set the Q variable to a unit-matrix
	//		 in other words, we loose the material axis orientation
	FEElasticMaterialPoint pt;
	pt.Init(true);

	// set the material point data
	pt.r0 = m_pMesh->Node(node.inode).m_r0;
	pt.rt = m_pMesh->Node(node.inode).m_rt;

	pt.F = node.Fi;
	pt.J = node.vi / node.Vi;

	pt.s = node.si;

	// if the material is incompressible we need to some additional values
	// since for tets the 3-field formulation is pointless, I really need to
	// change this.
	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pme);
	if (pmi)
	{
		// calculate volume ratio
		pt.avgJ = pt.J;

		// Calculate pressure. 
		pt.avgp = pmi->Up(pt.J);
	}

	// Calculate the tangent
	tens4ds C = pme->Tangent(pt);

	// if the material is incompressible we need to add
	// the dilatational part since that part is not calculated in the
	// Tangent function.
	if (pmi)
	{
		double dpdJ = pmi->Upp(pt.J);
		mat3dd I(1); // Identity
		C += dyad1s(I)*(pt.J*dpdJ);
	}

	// Next, we need to subtract the volumetric contribution Cvol
	if (m_bdev)
	{
		// subtract the isochoric component from C;
		// C = C - a*Ciso = C - (a*(C - Cvol)) = (1-a)*C + a*Cvol
		C = C*(1 - m_alpha) + Cvol(C, pt.s)*m_alpha;
	}
	else
	{
		C = C*(1 - m_alpha);
	}

	// extract the 'D' matrix
	double D[6][6] = {0};
	C.extract(D);

	// get the number of elements this nodes connects
	int NE = m_NEL.Valence(node.inode);
	FEElement** ppe = m_NEL.ElementList(node.inode);

	// loop over all the elements
	int i, j, ni, nj;
	double DB[6][3];
	double kij[3][3];

	matrix ke; ke.Create(NE*4*3, NE*4*3);

	// create the LM and the en array
	vector<int> LM; LM.resize(NE*4*3);
	vector<int> en; en.resize(NE*4  );
	vector<double> Ve; Ve.resize(NE);

	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& el = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(el);
		Ve[ni] = TetVolume(el.r0());

		for (i=0; i<4; ++i)
		{
			LM[ni*4*3+3*i  ] = el.LM()[3*i  ];
			LM[ni*4*3+3*i+1] = el.LM()[3*i+1];
			LM[ni*4*3+3*i+2] = el.LM()[3*i+2];

			en[ni*4+i] = el.m_node[i];
		}
	}

	// calculate B-matrices
	double (*Be)[6][3] = new double [NE*4][6][3];
	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& el = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(el);

		// calculate the jacobian
		double Ji[3][3];
		el.invjact(Ji, 0);

		// get the shape function derivatives
		const double* Gr, *Gs, *Gt;
		Gr = el.Gr(0);
		Gs = el.Gs(0);
		Gt = el.Gt(0);

		double Gx, Gy, Gz;
		for (j=0; j<4; ++j)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
			Gy = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
			Gz = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];

			Be[4*ni+j][0][0] = Gx; Be[4*ni+j][0][1] =  0; Be[4*ni+j][0][2] = 0;
			Be[4*ni+j][1][0] =  0; Be[4*ni+j][1][1] = Gy; Be[4*ni+j][1][2] = 0;
			Be[4*ni+j][2][0] =  0; Be[4*ni+j][2][1] =  0; Be[4*ni+j][2][2] = Gz;
			Be[4*ni+j][3][0] = Gy; Be[4*ni+j][3][1] = Gx; Be[4*ni+j][3][2] = 0; 
			Be[4*ni+j][4][0] =  0; Be[4*ni+j][4][1] = Gz; Be[4*ni+j][4][2] = Gy; 
			Be[4*ni+j][5][0] = Gz; Be[4*ni+j][5][1] =  0; Be[4*ni+j][5][2] = Gx;
		}
	}

	for (ni=0; ni<NE; ++ni)
	{
		FESolidElement& ei = dynamic_cast<FESolidElement&>(*ppe[ni]);
		UnpackElement(ei);

		// calculate element volume
		double Vi = Ve[ni];
		double wi = 0.25* Vi*node.vi/ node.Vi;

		// loop over the elements again
		for (nj=0; nj<NE; ++nj)
		{
			FESolidElement& ej = dynamic_cast<FESolidElement&>(*ppe[nj]);
			UnpackElement(ej);

			// calculate element volume
			double Vj = Ve[nj];
			double wj = 0.25* Vj / node.Vi;

			// We're ready to rock and roll!
			for (i=0; i<4; ++i)
			{
				double (&Bi)[6][3] = *(Be+(ni*4 + i));

				for (j=0; j<4; ++j)
				{
					double (&Bj)[6][3] = *(Be+(nj*4 + j));

					DB[0][0] = wj*(D[0][0]*Bj[0][0]+D[0][1]*Bj[1][0]+D[0][2]*Bj[2][0]+D[0][3]*Bj[3][0]+D[0][4]*Bj[4][0]+D[0][5]*Bj[5][0]);
					DB[0][1] = wj*(D[0][0]*Bj[0][1]+D[0][1]*Bj[1][1]+D[0][2]*Bj[2][1]+D[0][3]*Bj[3][1]+D[0][4]*Bj[4][1]+D[0][5]*Bj[5][1]);
					DB[0][2] = wj*(D[0][0]*Bj[0][2]+D[0][1]*Bj[1][2]+D[0][2]*Bj[2][2]+D[0][3]*Bj[3][2]+D[0][4]*Bj[4][2]+D[0][5]*Bj[5][2]);

					DB[1][0] = wj*(D[1][0]*Bj[0][0]+D[1][1]*Bj[1][0]+D[1][2]*Bj[2][0]+D[1][3]*Bj[3][0]+D[1][4]*Bj[4][0]+D[1][5]*Bj[5][0]);
					DB[1][1] = wj*(D[1][0]*Bj[0][1]+D[1][1]*Bj[1][1]+D[1][2]*Bj[2][1]+D[1][3]*Bj[3][1]+D[1][4]*Bj[4][1]+D[1][5]*Bj[5][1]);
					DB[1][2] = wj*(D[1][0]*Bj[0][2]+D[1][1]*Bj[1][2]+D[1][2]*Bj[2][2]+D[1][3]*Bj[3][2]+D[1][4]*Bj[4][2]+D[1][5]*Bj[5][2]);

					DB[2][0] = wj*(D[2][0]*Bj[0][0]+D[2][1]*Bj[1][0]+D[2][2]*Bj[2][0]+D[2][3]*Bj[3][0]+D[2][4]*Bj[4][0]+D[2][5]*Bj[5][0]);
					DB[2][1] = wj*(D[2][0]*Bj[0][1]+D[2][1]*Bj[1][1]+D[2][2]*Bj[2][1]+D[2][3]*Bj[3][1]+D[2][4]*Bj[4][1]+D[2][5]*Bj[5][1]);
					DB[2][2] = wj*(D[2][0]*Bj[0][2]+D[2][1]*Bj[1][2]+D[2][2]*Bj[2][2]+D[2][3]*Bj[3][2]+D[2][4]*Bj[4][2]+D[2][5]*Bj[5][2]);

					DB[3][0] = wj*(D[3][0]*Bj[0][0]+D[3][1]*Bj[1][0]+D[3][2]*Bj[2][0]+D[3][3]*Bj[3][0]+D[3][4]*Bj[4][0]+D[3][5]*Bj[5][0]);
					DB[3][1] = wj*(D[3][0]*Bj[0][1]+D[3][1]*Bj[1][1]+D[3][2]*Bj[2][1]+D[3][3]*Bj[3][1]+D[3][4]*Bj[4][1]+D[3][5]*Bj[5][1]);
					DB[3][2] = wj*(D[3][0]*Bj[0][2]+D[3][1]*Bj[1][2]+D[3][2]*Bj[2][2]+D[3][3]*Bj[3][2]+D[3][4]*Bj[4][2]+D[3][5]*Bj[5][2]);

					DB[4][0] = wj*(D[4][0]*Bj[0][0]+D[4][1]*Bj[1][0]+D[4][2]*Bj[2][0]+D[4][3]*Bj[3][0]+D[4][4]*Bj[4][0]+D[4][5]*Bj[5][0]);
					DB[4][1] = wj*(D[4][0]*Bj[0][1]+D[4][1]*Bj[1][1]+D[4][2]*Bj[2][1]+D[4][3]*Bj[3][1]+D[4][4]*Bj[4][1]+D[4][5]*Bj[5][1]);
					DB[4][2] = wj*(D[4][0]*Bj[0][2]+D[4][1]*Bj[1][2]+D[4][2]*Bj[2][2]+D[4][3]*Bj[3][2]+D[4][4]*Bj[4][2]+D[4][5]*Bj[5][2]);

					DB[5][0] = wj*(D[5][0]*Bj[0][0]+D[5][1]*Bj[1][0]+D[5][2]*Bj[2][0]+D[5][3]*Bj[3][0]+D[5][4]*Bj[4][0]+D[5][5]*Bj[5][0]);
					DB[5][1] = wj*(D[5][0]*Bj[0][1]+D[5][1]*Bj[1][1]+D[5][2]*Bj[2][1]+D[5][3]*Bj[3][1]+D[5][4]*Bj[4][1]+D[5][5]*Bj[5][1]);
					DB[5][2] = wj*(D[5][0]*Bj[0][2]+D[5][1]*Bj[1][2]+D[5][2]*Bj[2][2]+D[5][3]*Bj[3][2]+D[5][4]*Bj[4][2]+D[5][5]*Bj[5][2]);

					kij[0][0] = wi*(Bi[0][0]*DB[0][0]+Bi[1][0]*DB[1][0]+Bi[2][0]*DB[2][0]+Bi[3][0]*DB[3][0]+Bi[4][0]*DB[4][0]+Bi[5][0]*DB[5][0]);
					kij[0][1] = wi*(Bi[0][0]*DB[0][1]+Bi[1][0]*DB[1][1]+Bi[2][0]*DB[2][1]+Bi[3][0]*DB[3][1]+Bi[4][0]*DB[4][1]+Bi[5][0]*DB[5][1]);
					kij[0][2] = wi*(Bi[0][0]*DB[0][2]+Bi[1][0]*DB[1][2]+Bi[2][0]*DB[2][2]+Bi[3][0]*DB[3][2]+Bi[4][0]*DB[4][2]+Bi[5][0]*DB[5][2]);

					kij[1][0] = wi*(Bi[0][1]*DB[0][0]+Bi[1][1]*DB[1][0]+Bi[2][1]*DB[2][0]+Bi[3][1]*DB[3][0]+Bi[4][1]*DB[4][0]+Bi[5][1]*DB[5][0]);
					kij[1][1] = wi*(Bi[0][1]*DB[0][1]+Bi[1][1]*DB[1][1]+Bi[2][1]*DB[2][1]+Bi[3][1]*DB[3][1]+Bi[4][1]*DB[4][1]+Bi[5][1]*DB[5][1]);
					kij[1][2] = wi*(Bi[0][1]*DB[0][2]+Bi[1][1]*DB[1][2]+Bi[2][1]*DB[2][2]+Bi[3][1]*DB[3][2]+Bi[4][1]*DB[4][2]+Bi[5][1]*DB[5][2]);

					kij[2][0] = wi*(Bi[0][2]*DB[0][0]+Bi[1][2]*DB[1][0]+Bi[2][2]*DB[2][0]+Bi[3][2]*DB[3][0]+Bi[4][2]*DB[4][0]+Bi[5][2]*DB[5][0]);
					kij[2][1] = wi*(Bi[0][2]*DB[0][1]+Bi[1][2]*DB[1][1]+Bi[2][2]*DB[2][1]+Bi[3][2]*DB[3][1]+Bi[4][2]*DB[4][1]+Bi[5][2]*DB[5][1]);
					kij[2][2] = wi*(Bi[0][2]*DB[0][2]+Bi[1][2]*DB[1][2]+Bi[2][2]*DB[2][2]+Bi[3][2]*DB[3][2]+Bi[4][2]*DB[4][2]+Bi[5][2]*DB[5][2]);

					// copy to element stiffness matrix
					ke[ni*12+i*3  ][nj*12+j*3  ] = kij[0][0]; ke[ni*12+i*3  ][nj*12+j*3+1] = kij[0][1]; ke[ni*12+i*3  ][nj*12+j*3+2] = kij[0][2];
					ke[ni*12+i*3+1][nj*12+j*3  ] = kij[1][0]; ke[ni*12+i*3+1][nj*12+j*3+1] = kij[1][1]; ke[ni*12+i*3+1][nj*12+j*3+2] = kij[1][2];
					ke[ni*12+i*3+2][nj*12+j*3  ] = kij[2][0]; ke[ni*12+i*3+2][nj*12+j*3+1] = kij[2][1]; ke[ni*12+i*3+2][nj*12+j*3+2] = kij[2][2];
				}
			}
		}
	}

	// assemble the stiffness
	psolver->AssembleStiffness(en, LM, ke);

	// cleanup
	delete [] Be;
}

//-----------------------------------------------------------------------------
//! Calculates the element contribution to the global stiffness matrix
void FEUT4Domain::ElementalStiffnessMatrix(FESolidSolver *psolver)
{
	FEM& fem = psolver->m_fem;

	// element stiffness matrix
	matrix ke;

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		UnpackElement(el);

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.Create(ndof, ndof);
		ke.zero();

		// calculate the element stiffness matrix
		ElementStiffness(fem, el, ke);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, el.LM(), ke);
	}
}

//-----------------------------------------------------------------------------
void FEUT4Domain::ElementStiffness(FEM &fem, FESolidElement &el, matrix &ke)
{
	// TODO: I need to figure out how to deal with incompressible materials
	//       Incompressible materials require an additional dilatational 
	//       stiffness which is calculated differently from the geometric and
	//       material stiffnesses. 

	// calculate material stiffness (i.e. constitutive component)
	MaterialStiffness(fem, el, ke);

	// calculate geometrical stiffness
	GeometricalStiffness(el, ke);

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

void FEUT4Domain::GeometricalStiffness(FESolidElement &el, matrix &ke)
{
	int n, i, j;

	double Gx[8], Gy[8], Gz[8];
	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// nr of nodes
	int neln = el.Nodes();

	// nr of integration points
	int nint = el.GaussPoints();

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// calculate geometrical element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// get the material point data
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds s = pt.s;

		// we work with the deviatoric component only
		if (m_bdev)
		{
			s = s.dev()*m_alpha;
		}
		else
		{
			s = s*m_alpha;
		}

		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
					   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
					   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FEUT4Domain::MaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	assert(fem.m_pStep->m_nModule != FE_POROELASTIC);

	int i, i3, j, j3, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	// NOTE: hard-coding of hex elements!
	// Gx = dH/dx
	double Gx[8], Gy[8], Gz[8];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// jacobian
	double Ji[3][3], detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	FESolidMaterial* pmat = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
//		el.defgrad(pt.F, n);
//		pt.J = el.detF(n);
		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;

		// Calculate the tangent
		tens4ds C = pmat->Tangent(mp);

		// if the material is incompressible we need to add
		// the dilatational part since that part is not calculated in the
		// Tangent function.
		FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pmat);
		if (pmi)
		{
			double dpdJ = pmi->Upp(pt.J);
			mat3dd I(1); // Identity
			C += dyad1s(I)*(pt.J*dpdJ);
		}

		// Next, we need to subtract the volumetric contribution Cvol
		if (m_bdev)
		{
			// subtract the volumetric tensor from C;
			C = (C - Cvol(C, pt.s))*m_alpha;
		}
		else
		{
			C = C*m_alpha;
		}

		// extract the 'D' matrix
		C.extract(D);

		if (dynamic_cast<FEMicroMaterial*>(pmat))
		{
			// the micro-material screws up the currently unpacked elements
			// so I have to unpack the element data again
			UnpackElement(el);
		}

		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = Gx[i];
			Gyi = Gy[i];
			Gzi = Gz[i];

			for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = Gx[j];
				Gyj = Gy[j];
				Gzj = Gz[j];

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
