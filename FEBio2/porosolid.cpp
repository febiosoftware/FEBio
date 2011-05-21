#include "stdafx.h"
#include "FEElasticSolidDomain.h"
#include "FESolidSolver.h"
#include "FEPoroElastic.h"
#include "FEMicroMaterial.h"
#include "log.h"

//-----------------------------------------------------------------------------
void FEPoroSolidDomain::Residual(FESolidSolver* psolver, vector<double>& R)
{
	int i, j;

	FEM& fem = psolver->m_fem;

	// make sure we are in poro-mode
	assert(fem.m_pStep->m_nModule == FE_POROELASTIC);

	// element force vector
	vector<double> fe;

	vector<int> elm;

	int NE = m_Elem.size();
	for (i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

		// this element should not be rigid
		assert(!el.IsRigid());

		//! this element should not be UDG
		assert(el.Type() != FE_UDGHEX);

		// unpack the element
		UnpackElement(el);

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		// (This function is inherited from FEElasticSolidDomain)
		InternalForces(el, fe);

		// apply body forces
		// TODO: can we calculate body-forces with our formulation
		//       of biphasic theory
		/*
		if (fem.UseBodyForces())
		{
			BodyForces(fem, el, fe);
		}
		*/

		// get element equation numbers
		UnpackLM(el, elm);

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, elm, fe, R);

		// do poro-elastic forces
		FEMaterial* pm = fem.GetMaterial(el.GetMatID());
		assert(dynamic_cast<FEPoroElastic*>(pm) != 0);

		// calculate fluid internal work
		InternalFluidWork(fem, el, fe);

		// add fluid work to global residual
		int neln = el.Nodes();
		int J;
		for (j=0; j<neln; ++j)
		{
			J = elm[3*neln+j];
			if (J >= 0) R[J] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FEPoroSolidDomain::InternalFluidWork(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int i, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// jacobian
	double Ji[3][3], detJ, J0i[3][3];

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;

	vec3d vt[8];
	double pn[8];
	for (i=0; i<neln; ++i)
	{
		vt[i] = m_pMesh->Node(el.m_node[i]).m_vt;
		pn[i] = m_pMesh->Node(el.m_node[i]).m_pt;
	}

	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);

	// gauss-weights
	double* wg = el.GaussWeights();

	FEMesh& mesh = fem.m_mesh;
	
	vec3d rp[8];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
	}

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}

	zero(fe);

	// get the time step value
	double dt = fem.m_pStep->m_dt;

	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEPoroElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FEPoroElasticMaterialPoint>());

		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);

		// we need to calculate the divergence of v. To do this we use
		// the formula div(v) = 1/J*dJdt, where J = det(F)
		el.invjac0(J0i, n);
		
		// next we calculate the deformation gradient
		mat3d Fp;
		Fp.zero();
		
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		H = el.H(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			GX = J0i[0][0]*Gr[i]+J0i[1][0]*Gs[i]+J0i[2][0]*Gt[i];
			GY = J0i[0][1]*Gr[i]+J0i[1][1]*Gs[i]+J0i[2][1]*Gt[i];
			GZ = J0i[0][2]*Gr[i]+J0i[1][2]*Gs[i]+J0i[2][2]*Gt[i];
			
			Fp[0][0] += rp[i].x*GX; Fp[1][0] += rp[i].y*GX; Fp[2][0] += rp[i].z*GX;
			Fp[0][1] += rp[i].x*GY; Fp[1][1] += rp[i].y*GY; Fp[2][1] += rp[i].z*GY;
			Fp[0][2] += rp[i].x*GZ; Fp[1][2] += rp[i].y*GZ; Fp[2][2] += rp[i].z*GZ;
			
			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
		}

		// next we get the determinant
		double Jp = Fp.det();
		double J = ept.J;
		
		// and then finally
		double divv = ((J-Jp)/dt)/J;

		// get the flux
		vec3d& w = pt.m_w;

		// update force vector
		for (i=0; i<neln; ++i)
		{
			fe[i] -= dt*(B1[i]*w.x+B2[i]*w.y+B3[i]*w.z - divv*H[i])*detJ*wg[n];
		}
	}

	return true;
}


//-----------------------------------------------------------------------------

void FEPoroSolidDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;

	// element stiffness matrix
	matrix ke;

	vector<int> elm;

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// this element should not be rigid
		assert(!el.IsRigid());

		UnpackElement(el);

		// get the elements material
		FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
		assert(dynamic_cast<FEPoroElastic*>(pmat) != 0);

		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*4;
		ke.Create(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementPoroStiffness(fem, el, ke);

		// get the element's LM vector
		UnpackLM(el, elm);

		// TODO: the problem here is that the LM array that is returned by the UnpackElement
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackElement function so that it returns
		// the LM vector in the right order for poroelastic elements.
		vector<int> lm(ndof);
		for (int i=0; i<neln; ++i)
		{
			lm[4*i  ] = elm[3*i];
			lm[4*i+1] = elm[3*i+1];
			lm[4*i+2] = elm[3*i+2];
			lm[4*i+3] = elm[3*neln+i];
		}
				
		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEPoroSolidDomain::ElementPoroStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;

	// jacobian
	double Ji[3][3], detJ, J0i[3][3];

	// Bp-matrix
	vector<double> B1(neln), B2(neln), B3(neln);
	double tmp;

	// permeability tensor
	double k[3][3];

	// zero stiffness matrix
	ke.zero();

	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	SolidElementStiffness(fem, el, ks);

	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[4*i  ][4*j] = ks[3*i  ][3*j  ]; ke[4*i  ][4*j+1] = ks[3*i  ][3*j+1]; ke[4*i  ][4*j+2] = ks[3*i  ][3*j+2];
			ke[4*i+1][4*j] = ks[3*i+1][3*j  ]; ke[4*i+1][4*j+1] = ks[3*i+1][3*j+1]; ke[4*i+1][4*j+2] = ks[3*i+1][3*j+2];
			ke[4*i+2][4*j] = ks[3*i+2][3*j  ]; ke[4*i+2][4*j+1] = ks[3*i+2][3*j+1]; ke[4*i+2][4*j+2] = ks[3*i+2][3*j+2];
		}

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		clog.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}

	// if we use the symmetric version of the poro-implementation
	// we have to multiply some stiffness terms with the time step
	bool bsymm = fem.m_bsym_poro;
	double dt = fem.m_pStep->m_dt;

	FEMesh& mesh = fem.m_mesh;

	// gauss-weights
	double* gw = el.GaussWeights();

	vec3d rp[8], v[8];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
		v[i] = mesh.Node(el.m_node[i]).m_vt;
	}

	// we'll need this later
	const int T[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};

	// loop over gauss-points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate jacobian
		el.invjact(Ji, n);
		detJ = el.detJt(n);

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		H = el.H(n);

		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate Bp matrix
			B1[i] = Gx;
			B2[i] = Gy;
			B3[i] = Gz;
		}

		// get the permeability tensor
		pm->Permeability(k, mp);

		// calculate the Q = Bt*k*B matrix
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				tmp = dt*detJ*gw[n];
				ke[4*i+3][4*j+3] -= tmp*(B1[i]*k[0][0]+B2[i]*k[1][0]+B3[i]*k[2][0])*B1[j];
				ke[4*i+3][4*j+3] -= tmp*(B1[i]*k[0][1]+B2[i]*k[1][1]+B3[i]*k[2][1])*B2[j];
				ke[4*i+3][4*j+3] -= tmp*(B1[i]*k[0][2]+B2[i]*k[1][2]+B3[i]*k[2][2])*B3[j];
			}

		// calculate the G-matrix
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				tmp = detJ*gw[n]*H[j];
				ke[4*i  ][4*j+3] -= tmp*B1[i];
				ke[4*i+1][4*j+3] -= tmp*B2[i];
				ke[4*i+2][4*j+3] -= tmp*B3[i];
			}

		if (bsymm)
		{
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					tmp = detJ*gw[n]*H[j];
					ke[4*j+3][4*i  ] -= tmp*B1[i];
					ke[4*j+3][4*i+1] -= tmp*B2[i];
					ke[4*j+3][4*i+2] -= tmp*B3[i];
				}
		}
		else
		{
			// we need to calculate the divergence of v. To do this we use
			// the formula div(v) = 1/J*dJdt, where J = det(F)
			el.invjac0(J0i, n);

			// next we calculate the deformation gradient
			mat3d Fp, gradv;
			Fp.zero();
			gradv.zero();

			for (i=0; i<neln; ++i)
			{
				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !
				Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
				Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
				Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

				GX = J0i[0][0]*Gr[i]+J0i[1][0]*Gs[i]+J0i[2][0]*Gt[i];
				GY = J0i[0][1]*Gr[i]+J0i[1][1]*Gs[i]+J0i[2][1]*Gt[i];
				GZ = J0i[0][2]*Gr[i]+J0i[1][2]*Gs[i]+J0i[2][2]*Gt[i];

				Fp[0][0] += rp[i].x*GX; Fp[1][0] += rp[i].y*GX; Fp[2][0] += rp[i].z*GX;
				Fp[0][1] += rp[i].x*GY; Fp[1][1] += rp[i].y*GY; Fp[2][1] += rp[i].z*GY;
				Fp[0][2] += rp[i].x*GZ; Fp[1][2] += rp[i].y*GZ; Fp[2][2] += rp[i].z*GZ;

				gradv[0][0] += v[i].x*Gx; gradv[1][0] += v[i].y*Gx; gradv[2][0] += v[i].z*Gx;
				gradv[0][1] += v[i].x*Gy; gradv[1][1] += v[i].y*Gy; gradv[2][1] += v[i].z*Gy;
				gradv[0][2] += v[i].x*Gz; gradv[1][2] += v[i].y*Gz; gradv[2][2] += v[i].z*Gz;
			}

			// next we get the determinant
			double Jp = Fp.det();
			double J = pt.J;

			// and then finally
			double divv = ((J-Jp)/dt)/J;

			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					tmp = dt*detJ*gw[n]*H[i];
					ke[4*i+3][4*j  ] -= tmp*(B1[j]*(divv+1/dt) - (gradv[0][0]*B1[j] + gradv[0][1]*B2[j] + gradv[0][2]*B3[j]));
					ke[4*i+3][4*j+1] -= tmp*(B2[j]*(divv+1/dt) - (gradv[1][0]*B1[j] + gradv[1][1]*B2[j] + gradv[1][2]*B3[j]));
					ke[4*i+3][4*j+2] -= tmp*(B3[j]*(divv+1/dt) - (gradv[2][0]*B1[j] + gradv[2][1]*B2[j] + gradv[2][2]*B3[j]));
				}

			FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();
			vec3d Dp = pt.m_gradp;

			tens4ds K = pm->Tangent_Permeability(mp);
			double D[6][6]; K.extract(D);
			double BKB[3][3];
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					// The multiplication BKB has been expanded since it was observed
					// that this creates a significant boost in speed
					{
						BKB[0][0]  = B1[i]*(B1[j]*D[T[0][0]][T[0][0]] + B2[j]*D[T[0][0]][T[0][1]] + B3[j]*D[T[0][0]][T[0][2]]); 
						BKB[0][0] += B2[i]*(B1[j]*D[T[1][0]][T[0][0]] + B2[j]*D[T[1][0]][T[0][1]] + B3[j]*D[T[1][0]][T[0][2]]); 
						BKB[0][0] += B3[i]*(B1[j]*D[T[2][0]][T[0][0]] + B2[j]*D[T[2][0]][T[0][1]] + B3[j]*D[T[2][0]][T[0][2]]); 

						BKB[0][1]  = B1[i]*(B1[j]*D[T[0][0]][T[1][0]] + B2[j]*D[T[0][0]][T[1][1]] + B3[j]*D[T[0][0]][T[1][2]]);
						BKB[0][1] += B2[i]*(B1[j]*D[T[1][0]][T[1][0]] + B2[j]*D[T[1][0]][T[1][1]] + B3[j]*D[T[1][0]][T[1][2]]);
						BKB[0][1] += B3[i]*(B1[j]*D[T[2][0]][T[1][0]] + B2[j]*D[T[2][0]][T[1][1]] + B3[j]*D[T[2][0]][T[1][2]]);

						BKB[0][2]  = B1[i]*(B1[j]*D[T[0][0]][T[2][0]] + B2[j]*D[T[0][0]][T[2][1]] + B3[j]*D[T[0][0]][T[2][2]]);
						BKB[0][2] += B2[i]*(B1[j]*D[T[1][0]][T[2][0]] + B2[j]*D[T[1][0]][T[2][1]] + B3[j]*D[T[1][0]][T[2][2]]);
						BKB[0][2] += B3[i]*(B1[j]*D[T[2][0]][T[2][0]] + B2[j]*D[T[2][0]][T[2][1]] + B3[j]*D[T[2][0]][T[2][2]]);

						BKB[1][0]  = B1[i]*(B1[j]*D[T[0][1]][T[0][0]] + B2[j]*D[T[0][1]][T[0][1]] + B3[j]*D[T[0][1]][T[0][2]]);
						BKB[1][0] += B2[i]*(B1[j]*D[T[1][1]][T[0][0]] + B2[j]*D[T[1][1]][T[0][1]] + B3[j]*D[T[1][1]][T[0][2]]);
						BKB[1][0] += B3[i]*(B1[j]*D[T[2][1]][T[0][0]] + B2[j]*D[T[2][1]][T[0][1]] + B3[j]*D[T[2][1]][T[0][2]]);

						BKB[1][1]  = B1[i]*(B1[j]*D[T[0][1]][T[1][0]] + B2[j]*D[T[0][1]][T[1][1]] + B3[j]*D[T[0][1]][T[1][2]]);
						BKB[1][1] += B2[i]*(B1[j]*D[T[1][1]][T[1][0]] + B2[j]*D[T[1][1]][T[1][1]] + B3[j]*D[T[1][1]][T[1][2]]);
						BKB[1][1] += B3[i]*(B1[j]*D[T[2][1]][T[1][0]] + B2[j]*D[T[2][1]][T[1][1]] + B3[j]*D[T[2][1]][T[1][2]]);

						BKB[1][2]  = B1[i]*(B1[j]*D[T[0][1]][T[2][0]] + B2[j]*D[T[0][1]][T[2][1]] + B3[j]*D[T[0][1]][T[2][2]]);
						BKB[1][2] += B2[i]*(B1[j]*D[T[1][1]][T[2][0]] + B2[j]*D[T[1][1]][T[2][1]] + B3[j]*D[T[1][1]][T[2][2]]);
						BKB[1][2] += B3[i]*(B1[j]*D[T[2][1]][T[2][0]] + B2[j]*D[T[2][1]][T[2][1]] + B3[j]*D[T[2][1]][T[2][2]]);

						BKB[2][0]  = B1[i]*(B1[j]*D[T[0][2]][T[0][0]] + B2[j]*D[T[0][2]][T[0][1]] + B3[j]*D[T[0][2]][T[0][2]]);
						BKB[2][0] += B2[i]*(B1[j]*D[T[1][2]][T[0][0]] + B2[j]*D[T[1][2]][T[0][1]] + B3[j]*D[T[1][2]][T[0][2]]);
						BKB[2][0] += B3[i]*(B1[j]*D[T[2][2]][T[0][0]] + B2[j]*D[T[2][2]][T[0][1]] + B3[j]*D[T[2][2]][T[0][2]]);

						BKB[2][1]  = B1[i]*(B1[j]*D[T[0][2]][T[1][0]] + B2[j]*D[T[0][2]][T[1][1]] + B3[j]*D[T[0][2]][T[1][2]]); 
						BKB[2][1] += B2[i]*(B1[j]*D[T[1][2]][T[1][0]] + B2[j]*D[T[1][2]][T[1][1]] + B3[j]*D[T[1][2]][T[1][2]]); 
						BKB[2][1] += B3[i]*(B1[j]*D[T[2][2]][T[1][0]] + B2[j]*D[T[2][2]][T[1][1]] + B3[j]*D[T[2][2]][T[1][2]]); 

						BKB[2][2]  = B1[i]*(B1[j]*D[T[0][2]][T[2][0]] + B2[j]*D[T[0][2]][T[2][1]] + B3[j]*D[T[0][2]][T[2][2]]);
						BKB[2][2] += B2[i]*(B1[j]*D[T[1][2]][T[2][0]] + B2[j]*D[T[1][2]][T[2][1]] + B3[j]*D[T[1][2]][T[2][2]]);
						BKB[2][2] += B3[i]*(B1[j]*D[T[2][2]][T[2][0]] + B2[j]*D[T[2][2]][T[2][1]] + B3[j]*D[T[2][2]][T[2][2]]);
					}

					tmp = dt*detJ*gw[n];
					ke[4*i+3][4*j  ] -= tmp*( BKB[0][0]*Dp.x + BKB[0][1]*Dp.y + BKB[0][2]*Dp.z);
					ke[4*i+3][4*j+1] -= tmp*( BKB[1][0]*Dp.x + BKB[1][1]*Dp.y + BKB[1][2]*Dp.z);
					ke[4*i+3][4*j+2] -= tmp*( BKB[2][0]*Dp.x + BKB[2][1]*Dp.y + BKB[2][2]*Dp.z);
				}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FEPoroSolidDomain::SolidElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	PoroMaterialStiffness(fem, el, ke);

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
//! Calculates element material stiffness element matrix

void FEPoroSolidDomain::PoroMaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	assert(fem.m_pStep->m_nModule == FE_POROELASTIC);

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

	// nodal pressures
	double pn[8];
	for (i=0; i<neln; ++i) pn[i] = m_pMesh->Node(el.m_node[i]).m_pt;

	// see if this is a poroelastic material
	FESolidMaterial* pmat = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));
	assert(dynamic_cast<FEPoroElastic*>(pmat));

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

		// evaluate fluid pressure at gauss-point
		FEPoroElasticMaterialPoint& ppt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());
		ppt.m_p = el.Evaluate(pn, n);

		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
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


//-----------------------------------------------------------------------------
void FEPoroSolidDomain::UpdateStresses(FEM &fem)
{
	int i, n;
	int nint, neln;
	double* gw;
	vec3d r0[8], rt[8];
	double pn[8];

	assert(fem.m_pStep->m_nModule == FE_POROELASTIC);
	FEMesh& mesh = *m_pMesh;

	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		assert(!el.IsRigid());

		assert(el.Type() != FE_UDGHEX);

		// unpack the element data
		UnpackElement(el);

		// get the number of integration points
		nint = el.GaussPoints();

		// number of nodes
		neln = el.Nodes();

		// nodal coordinates
		for (int j=0; j<neln; ++j)
		{
			r0[j] = mesh.Node(el.m_node[j]).m_r0;
			rt[j] = mesh.Node(el.m_node[j]).m_rt;
			pn[j] = mesh.Node(el.m_node[j]).m_pt;
		}

		// get the integration weights
		gw = el.GaussWeights();

		// get the material
		FESolidMaterial* pm = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));

		// extract the elastic component
		FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());

		assert(dynamic_cast<FEPoroElastic*>(pm) != 0);

		// loop over the integration points and calculate
		// the stress at the integration point
		for (n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

			// material point coordinates
			// TODO: I'm not entirly happy with this solution
			//		 since the material point coordinates are used by most materials.
			pt.r0 = el.Evaluate(r0, n);
			pt.rt = el.Evaluate(rt, n);
			
			// get the deformation gradient and determinant
			pt.J = el.defgrad(pt.F, n);

			// three-field element variables
			pt.avgJ = el.m_eJ;
			pt.avgp = el.m_ep;

			// poroelasticity data
			FEPoroElasticMaterialPoint& ppt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());

			// evaluate fluid pressure at gauss-point
			ppt.m_p = el.Evaluate(pn, n);

			// calculate the gradient of p at gauss-point
			ppt.m_gradp = el.gradient(pn, n);

			// calculate the stress at this material point
			pt.s = pm->Stress(mp);

			if (dynamic_cast<FEMicroMaterial*>(pme))
			{
				// the micro-material screws up the currently unpacked elements
				// so I have to unpack the element data again
				UnpackElement(el);
			}

			// for poroelastic materials also update the fluid flux
			FEPoroElastic* pmat = dynamic_cast<FEPoroElastic*>(pm);
			ppt.m_w = pmat->Flux(mp);
			ppt.m_pa = ppt.m_p;		// for consistency with biphasic material
		}
	}
}
