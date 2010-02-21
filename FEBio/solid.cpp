#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEDomain.h"
#include "log.h"
#include "FEPoroElastic.h"
#include "FEMicroMaterial.h"

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FESolidDomain::InternalForces(FESolidElement& el, vector<double>& fe)
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
		mat3ds& s = pt.s;

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
//! calculates the body forces

void FESolidDomain::BodyForces(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int i, n;
	double *H;

	// jacobian
	double detJ, dens;

	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));

	dens = pme->Density();

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double* gw = el.GaussWeights();

	// loop over integration points
	vec3d g = fem.m_acc*dens;
	for (n=0; n<nint; ++n)
	{
		detJ = el.detJ0(n)*gw[n];

		H = el.H(n);

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += H[i]*g.x*detJ;
			fe[3*i+1] += H[i]*g.y*detJ;
			fe[3*i+2] += H[i]*g.z*detJ;
		}						
	}
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for enhanced strain
//! solid elements.

void FESolidDomain::UDGInternalForces(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	// make sure this element is of the correct type
	assert(el.Type() == FE_UDGHEX);

	// get the stress data
	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	mat3ds& s = pt.s;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	m_pMesh->AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	m_pMesh->AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = m_pMesh->HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx, Gy, Gz;

	// calculate the internal force
	for (int i=0; i<8; ++i)
	{
		Gx = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];

		fe[3*i  ] -= ve*(Gx*s.xx() + Gy*s.xy() + Gz*s.xz());
		fe[3*i+1] -= ve*(Gx*s.xy() + Gy*s.yy() + Gz*s.yz());
		fe[3*i+2] -= ve*(Gx*s.xz() + Gy*s.yz() + Gz*s.zz());
	}

	// add hourglass forces
	UDGHourglassForces(fem, el, fe);
}


//-----------------------------------------------------------------------------
//! calculates the hourglass forces

void FESolidDomain::UDGHourglassForces(FEM& fem, FESolidElement &el, vector<double> &fe)
{
	int i;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	vec3d* r0 = el.r0();
	vec3d* rt = el.rt();

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	double GX[8], GY[8], GZ[8];
	m_pMesh->AvgCartDerivs(el, GX, GY, GZ);

	mat3d F;
	m_pMesh->AvgDefGrad(el, F, GX, GY, GZ);
/*
	double Ji0[3][3], Jt[3][3];
	el.invjac0(Ji0, 0);
	el.jact(Jt, 0);

	F[0][0] = Ji0[0][0]*Jt[0][0]+Ji0[0][1]*Jt[1][0]+Ji0[0][2]*Jt[2][0];
	F[0][1] = Ji0[0][0]*Jt[0][1]+Ji0[0][1]*Jt[1][1]+Ji0[0][2]*Jt[2][1];
	F[0][2] = Ji0[0][0]*Jt[0][2]+Ji0[0][1]*Jt[1][2]+Ji0[0][2]*Jt[2][2];

	F[1][0] = Ji0[1][0]*Jt[0][0]+Ji0[1][1]*Jt[1][0]+Ji0[1][2]*Jt[2][0];
	F[1][1] = Ji0[1][0]*Jt[0][1]+Ji0[1][1]*Jt[1][1]+Ji0[1][2]*Jt[2][1];
	F[1][2] = Ji0[1][0]*Jt[0][2]+Ji0[1][1]*Jt[1][2]+Ji0[1][2]*Jt[2][2];

	F[2][0] = Ji0[2][0]*Jt[0][0]+Ji0[2][1]*Jt[1][0]+Ji0[2][2]*Jt[2][0];
	F[2][1] = Ji0[2][0]*Jt[0][1]+Ji0[2][1]*Jt[1][1]+Ji0[2][2]*Jt[2][1];
	F[2][2] = Ji0[2][0]*Jt[0][2]+Ji0[2][1]*Jt[1][2]+Ji0[2][2]*Jt[2][2];
*/
	double u4 = 0, u5 = 0, u6 = 0, u7 = 0;
	double v4 = 0, v5 = 0, v6 = 0, v7 = 0;
	double w4 = 0, w5 = 0, w6 = 0, w7 = 0;

	u4 = x4 - (F[0][0]*X4 + F[0][1]*Y4 + F[0][2]*Z4);
	v4 = y4 - (F[1][0]*X4 + F[1][1]*Y4 + F[1][2]*Z4);
	w4 = z4 - (F[2][0]*X4 + F[2][1]*Y4 + F[2][2]*Z4);

	u5 = x5 - (F[0][0]*X5 + F[0][1]*Y5 + F[0][2]*Z5);
	v5 = y5 - (F[1][0]*X5 + F[1][1]*Y5 + F[1][2]*Z5);
	w5 = z5 - (F[2][0]*X5 + F[2][1]*Y5 + F[2][2]*Z5);

	u6 = x6 - (F[0][0]*X6 + F[0][1]*Y6 + F[0][2]*Z6);
	v6 = y6 - (F[1][0]*X6 + F[1][1]*Y6 + F[1][2]*Z6);
	w6 = z6 - (F[2][0]*X6 + F[2][1]*Y6 + F[2][2]*Z6);

	u7 = x7 - (F[0][0]*X7 + F[0][1]*Y7 + F[0][2]*Z7);
	v7 = y7 - (F[1][0]*X7 + F[1][1]*Y7 + F[1][2]*Z7);
	w7 = z7 - (F[2][0]*X7 + F[2][1]*Y7 + F[2][2]*Z7);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass forces
	double hg = fem.m_pStep->m_hg;

	for (i=0; i<8; ++i)
	{
		fe[3*i  ] -= hg*(g4[i]*u4 + g5[i]*u5 + g6[i]*u6 + g7[i]*u7);
		fe[3*i+1] -= hg*(g4[i]*v4 + g5[i]*v5 + g6[i]*v6 + g7[i]*v7);
		fe[3*i+2] -= hg*(g4[i]*w4 + g5[i]*w5 + g6[i]*w6 + g7[i]*w7);
	}
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces due to the fluid work
//! Note that we only use the first n entries in fe, where n is the number
//! of nodes

bool FESolidDomain::InternalFluidWork(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int i, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// jacobian
	double Ji[3][3], detJ, J0i[3][3];

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;

	vec3d* vt = el.vt();
	double* pn = el.pt();

	// Bp-matrix
	vector<double[3]> B(neln);

	// gauss-weights
	double* wg = el.GaussWeights();

	FEMesh& mesh = fem.m_mesh;
	
	vec3d rp[8];
	for (i=0; i<neln; ++i) 
	{
		rp[i] = mesh.Node(el.m_node[i]).m_rp;
	}

	// get the logfile
	Logfile& log = GetLogfile();

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		log.printbox("FATAL ERROR", "Incorrect material type\n");
		return false;
	}

	fe.zero();

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
			B[i][0] = Gx;
			B[i][1] = Gy;
			B[i][2] = Gz;
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
			fe[i] -= dt*(B[i][0]*w.x+B[i][1]*w.y+B[i][2]*w.z - divv*H[i])*detJ*wg[n];
		}
	}

	return true;
}


//-----------------------------------------------------------------------------
//! calculates dilatational element stiffness component for element iel

void FESolidDomain::DilatationalStiffness(FEM& fem, FESolidElement& elem, matrix& ke)
{
	int i, j, n;

	const int nint = elem.GaussPoints();
	const int neln = elem.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEElasticMaterial* pm = fem.GetElasticMaterial(elem.GetMatID());

	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pm);
	assert(pmi);

	// average global derivatives
	vector<double> gradN(3*neln);
	gradN.zero();

	// initial element volume
	double Ve = 0;

	// global derivatives of shape functions
	double Gx, Gy, Gz;
	const double *gw = elem.GaussWeights();

	// jacobian
	double Ji[3][3], detJt, detJ0;

	double *Gr, *Gs, *Gt;

	// repeat over gauss-points
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJ0 = elem.detJ0(n);
		detJt = elem.detJt(n);
		elem.invjact(Ji, n);

		detJt *= gw[n];

		Ve += detJ0*gw[n];

		Gr = elem.Gr(n);
		Gs = elem.Gs(n);
		Gt = elem.Gt(n);

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		for (i=0; i<neln; ++i)
		{
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			gradN[3*i  ] += Gx*detJt;
			gradN[3*i+1] += Gy*detJt;
			gradN[3*i+2] += Gz*detJt;
		}
	}

	// get effective modulus
	double k = pmi->Upp(elem.m_eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += elem.m_Lk*pmi->hpp(elem.m_eJ);

	// divide by initial volume
	k /= Ve;

	// calculate dilatational stiffness component
	// we only calculate the upper triangular part
	// since ke is symmetric.
	for (i=0; i<ndof; ++i)
		for (j=i; j<ndof; ++j)
			ke[i][j] += k*gradN[i]*gradN[j];
}


//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n

void FESolidDomain::GeometricalStiffness(FESolidElement &el, matrix &ke)
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
		mat3ds& s = pt.s;

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

void FESolidDomain::MaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
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

	// see if this is a poroelastic material
	FESolidMaterial* pmat = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));
	bool bporo = false;
	if ((fem.m_pStep->m_nModule == FE_POROELASTIC) && (dynamic_cast<FEPoroElastic*>(pmat))) bporo = true;

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

		if (bporo)
		{
			FEPoroElasticMaterialPoint& pt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());

			// evaluate fluid pressure at gauss-point
			pt.m_p = el.Evaluate(el.pt(), n);
		}

		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
		C.extract(D);

		if (dynamic_cast<FEMicroMaterial*>(pmat))
		{
			// the micro-material screws up the currently unpacked elements
			// so I have to unpack the element data again
			fem.m_mesh.UnpackElement(el);
		}

/*		if (m_fem.GetDebugFlag())
		{
			tens4ds t(D);
			if (IsPositiveDefinite(t) == false)
			{
				m_fem.m_log.printbox("WARNING", "Elasticity tensor is not positive-definite for\nelement %d at integration point %d.", el.m_nID, n+1);
			}
		}
*/
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
//! calculates the hourglass stiffness for UDG hex elements

void FESolidDomain::UDGHourglassStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	vec3d* r0 = el.r0();
	vec3d* rt = el.rt();

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	FEMesh& mesh = fem.m_mesh;

	double GX[8], GY[8], GZ[8];
	mesh.AvgCartDerivs(el, GX, GY, GZ);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass stiffness
	double hg = fem.m_pStep->m_hg;

	double kab;

	for (i=0; i<8; ++i)
	{
		for (j=i; j<8; ++j)
		{
			kab = hg*(g4[i]*g4[j] + g5[i]*g5[j] + g6[i]*g6[j] + g7[i]*g7[j]);

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
	}
}


//-----------------------------------------------------------------------------

void FESolidDomain::UDGDilatationalStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j;

	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEElasticMaterial* pm = fem.GetElasticMaterial(el.GetMatID());

	FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(pm);
	assert(pmi);

	FEMesh& mesh = fem.m_mesh;

	// calculate the average cartesian derivatives
	double Gx[8], Gy[8], Gz[8];
	mesh.AvgCartDerivs(el, Gx, Gy, Gz, 1);

	// calculate element volume
	double ve = mesh.HexVolume(el, 1);
	double Ve = mesh.HexVolume(el, 0);

	// get effective modulus
	double k = pmi->Upp(el.m_eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += el.m_Lk*pmi->hpp(el.m_eJ);

	// multiply with volume
	k *= (ve*ve)/Ve;

	// calculate dilatational stiffness component
	// we only calculate the upper triangular part
	// since ke is symmetric.
	for (i=0; i<8; ++i)
		for (j=i; j<8; ++j)
		{
			ke[3*i  ][3*j  ] += k*Gx[i]*Gx[j];
			ke[3*i  ][3*j+1] += k*Gx[i]*Gy[j];
			ke[3*i  ][3*j+2] += k*Gx[i]*Gz[j];

			ke[3*i+1][3*j  ] += k*Gy[i]*Gx[j];
			ke[3*i+1][3*j+1] += k*Gy[i]*Gy[j];
			ke[3*i+1][3*j+2] += k*Gy[i]*Gz[j];

			ke[3*i+2][3*j  ] += k*Gz[i]*Gx[j];
			ke[3*i+2][3*j+1] += k*Gz[i]*Gy[j];
			ke[3*i+2][3*j+2] += k*Gz[i]*Gz[j];
		}
}


//-----------------------------------------------------------------------------

void FESolidDomain::UDGGeometricalStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j;

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// get the material point data
	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// element's Cauchy-stress tensor at gauss point n
	// s is the voight vector
	mat3ds& s = pt.s;

	FEMesh& mesh = fem.m_mesh;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	mesh.AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	mesh.AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = mesh.HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	for (i=0; i<8; ++i)
		for (j=i; j<8; ++j)
		{
			kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
				   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
				   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*ve;

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
}


//-----------------------------------------------------------------------------

void FESolidDomain::UDGMaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	// make sure we have the right element type
	assert(el.Type() == FE_UDGHEX);

	int i, i3, j, j3;

	// Get the current element's data
	const int neln = el.Nodes();

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// see if this is a poroelastic material
	FESolidMaterial* pmat = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));
	bool bporo = false;
	if ((fem.m_pStep->m_nModule == FE_POROELASTIC) && (dynamic_cast<FEPoroElastic*>(pmat))) bporo = true;
	
	// for now we do not allow this element to be used in a poroelastic simulation
	assert(bporo == false);

	FEMaterialPoint& mp = *el.m_State[0];
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	FEMesh& mesh = fem.m_mesh;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	mesh.AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	mesh.AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = mesh.HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	// setup the material point
	// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
//	el.defgrad(pt.F, n);
//	pt.J = el.detF(n);
	pt.avgJ = el.m_eJ;
	pt.avgp = el.m_ep;

	// get the 'D' matrix
	tens4ds C = pmat->Tangent(mp);
	C.extract(D);

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

			ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*ve;
			ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*ve;
			ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*ve;

			ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*ve;
			ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*ve;
			ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*ve;

			ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*ve;
			ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*ve;
			ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*ve;
		}
	}
}


//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FESolidDomain::ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	// see if the material is incompressible
	FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());
	bool bdilst = false;
	if (dynamic_cast<FEIncompressibleMaterial*>(pme)) bdilst = true;

	if (el.Type() == FE_UDGHEX)
	{
		// calculate material stiffness
		UDGMaterialStiffness(fem, el, ke);

		// calculate geometrical stiffness
		UDGGeometricalStiffness(fem, el, ke);

		// Calculate dilatational stiffness, if necessary
		if (bdilst) UDGDilatationalStiffness(fem, el, ke);

		// add hourglass stiffness
		UDGHourglassStiffness(fem, el, ke);
	}
	else
	{
		// calculate material stiffness (i.e. constitutive component)
		MaterialStiffness(fem, el, ke);

		// calculate geometrical stiffness
		GeometricalStiffness(el, ke);

		// Calculate dilatational stiffness, if necessary
		if (bdilst) DilatationalStiffness(fem, el, ke);
	}

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
//! calculates element inertial stiffness matrix

void FESolidDomain::ElementInertialStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;

	// shape functions
	double* H;

	// jacobian
	double detJ0;

	// get the material
	FESolidMaterial* pm = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));

	double a = 4.0 / (fem.m_pStep->m_dt*fem.m_pStep->m_dt);
	double d = pm->Density();
	double kab;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		H = el.H(n);
		detJ0 = el.detJ0(n)*gw[n];
		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				kab = a*H[i]*H[j]*detJ0*d;
				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}	
	}

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];

}


//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel

bool FESolidDomain::ElementPoroStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double *Gr, *Gs, *Gt, *H;
	double Gx, Gy, Gz, GX, GY, GZ;

	// jacobian
	double Ji[3][3], detJ, J0i[3][3];

	// Bp-matrix
	vector<double[3]> B(neln);
	double tmp;

	// permeability tensor
	double k[3][3];

	// zero stiffness matrix
	ke.zero();

	// calculate solid stiffness matrix
	int ndof = 3*el.Nodes();
	matrix ks(ndof, ndof); ks.zero();
	ElementStiffness(fem, el, ks);

	// copy solid stiffness matrix into ke
	for (i=0; i<neln; ++i)
		for (j=0; j<neln; ++j)
		{
			ke[4*i  ][4*j] = ks[3*i  ][3*j  ]; ke[4*i  ][4*j+1] = ks[3*i  ][3*j+1]; ke[4*i  ][4*j+2] = ks[3*i  ][3*j+2];
			ke[4*i+1][4*j] = ks[3*i+1][3*j  ]; ke[4*i+1][4*j+1] = ks[3*i+1][3*j+1]; ke[4*i+1][4*j+2] = ks[3*i+1][3*j+2];
			ke[4*i+2][4*j] = ks[3*i+2][3*j  ]; ke[4*i+2][4*j+1] = ks[3*i+2][3*j+1]; ke[4*i+2][4*j+2] = ks[3*i+2][3*j+2];
		}

	// get the logfile
	Logfile& log = GetLogfile();

	// get the element's material
	FEPoroElastic* pm = dynamic_cast<FEPoroElastic*> (fem.GetMaterial(el.GetMatID()));
	if (pm == 0)
	{
		log.printbox("FATAL ERROR", "Incorrect material type\n");
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
			B[i][0] = Gx;
			B[i][1] = Gy;
			B[i][2] = Gz;
		}

		// get the permeability tensor
		pm->Permeability(k, mp);

		// calculate the Q = Bt*k*B matrix
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				tmp = dt*detJ*gw[n];
				ke[4*i+3][4*j+3] -= tmp*(B[i][0]*k[0][0]+B[i][1]*k[1][0]+B[i][2]*k[2][0])*B[j][0];
				ke[4*i+3][4*j+3] -= tmp*(B[i][0]*k[0][1]+B[i][1]*k[1][1]+B[i][2]*k[2][1])*B[j][1];
				ke[4*i+3][4*j+3] -= tmp*(B[i][0]*k[0][2]+B[i][1]*k[1][2]+B[i][2]*k[2][2])*B[j][2];
			}

		// calculate the G-matrix
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				tmp = detJ*gw[n]*H[j];
				ke[4*i  ][4*j+3] -= tmp*B[i][0];
				ke[4*i+1][4*j+3] -= tmp*B[i][1];
				ke[4*i+2][4*j+3] -= tmp*B[i][2];
			}

		if (bsymm)
		{
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					tmp = detJ*gw[n]*H[j];
					ke[4*j+3][4*i  ] -= tmp*B[i][0];
					ke[4*j+3][4*i+1] -= tmp*B[i][1];
					ke[4*j+3][4*i+2] -= tmp*B[i][2];
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
					ke[4*i+3][4*j  ] -= tmp*(B[j][0]*(divv+1/dt) - (gradv[0][0]*B[j][0] + gradv[0][1]*B[j][1] + gradv[0][2]*B[j][2]));
					ke[4*i+3][4*j+1] -= tmp*(B[j][1]*(divv+1/dt) - (gradv[1][0]*B[j][0] + gradv[1][1]*B[j][1] + gradv[1][2]*B[j][2]));
					ke[4*i+3][4*j+2] -= tmp*(B[j][2]*(divv+1/dt) - (gradv[2][0]*B[j][0] + gradv[2][1]*B[j][1] + gradv[2][2]*B[j][2]));
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
						BKB[0][0]  = B[i][0]*(B[j][0]*D[T[0][0]][T[0][0]] + B[j][1]*D[T[0][0]][T[0][1]] + B[j][2]*D[T[0][0]][T[0][2]]); 
						BKB[0][0] += B[i][1]*(B[j][0]*D[T[1][0]][T[0][0]] + B[j][1]*D[T[1][0]][T[0][1]] + B[j][2]*D[T[1][0]][T[0][2]]); 
						BKB[0][0] += B[i][2]*(B[j][0]*D[T[2][0]][T[0][0]] + B[j][1]*D[T[2][0]][T[0][1]] + B[j][2]*D[T[2][0]][T[0][2]]); 

						BKB[0][1]  = B[i][0]*(B[j][0]*D[T[0][0]][T[1][0]] + B[j][1]*D[T[0][0]][T[1][1]] + B[j][2]*D[T[0][0]][T[1][2]]);
						BKB[0][1] += B[i][1]*(B[j][0]*D[T[1][0]][T[1][0]] + B[j][1]*D[T[1][0]][T[1][1]] + B[j][2]*D[T[1][0]][T[1][2]]);
						BKB[0][1] += B[i][2]*(B[j][0]*D[T[2][0]][T[1][0]] + B[j][1]*D[T[2][0]][T[1][1]] + B[j][2]*D[T[2][0]][T[1][2]]);

						BKB[0][2]  = B[i][0]*(B[j][0]*D[T[0][0]][T[2][0]] + B[j][1]*D[T[0][0]][T[2][1]] + B[j][2]*D[T[0][0]][T[2][2]]);
						BKB[0][2] += B[i][1]*(B[j][0]*D[T[1][0]][T[2][0]] + B[j][1]*D[T[1][0]][T[2][1]] + B[j][2]*D[T[1][0]][T[2][2]]);
						BKB[0][2] += B[i][2]*(B[j][0]*D[T[2][0]][T[2][0]] + B[j][1]*D[T[2][0]][T[2][1]] + B[j][2]*D[T[2][0]][T[2][2]]);

						BKB[1][0]  = B[i][0]*(B[j][0]*D[T[0][1]][T[0][0]] + B[j][1]*D[T[0][1]][T[0][1]] + B[j][2]*D[T[0][1]][T[0][2]]);
						BKB[1][0] += B[i][1]*(B[j][0]*D[T[1][1]][T[0][0]] + B[j][1]*D[T[1][1]][T[0][1]] + B[j][2]*D[T[1][1]][T[0][2]]);
						BKB[1][0] += B[i][2]*(B[j][0]*D[T[2][1]][T[0][0]] + B[j][1]*D[T[2][1]][T[0][1]] + B[j][2]*D[T[2][1]][T[0][2]]);

						BKB[1][1]  = B[i][0]*(B[j][0]*D[T[0][1]][T[1][0]] + B[j][1]*D[T[0][1]][T[1][1]] + B[j][2]*D[T[0][1]][T[1][2]]);
						BKB[1][1] += B[i][1]*(B[j][0]*D[T[1][1]][T[1][0]] + B[j][1]*D[T[1][1]][T[1][1]] + B[j][2]*D[T[1][1]][T[1][2]]);
						BKB[1][1] += B[i][2]*(B[j][0]*D[T[2][1]][T[1][0]] + B[j][1]*D[T[2][1]][T[1][1]] + B[j][2]*D[T[2][1]][T[1][2]]);

						BKB[1][2]  = B[i][0]*(B[j][0]*D[T[0][1]][T[2][0]] + B[j][1]*D[T[0][1]][T[2][1]] + B[j][2]*D[T[0][1]][T[2][2]]);
						BKB[1][2] += B[i][1]*(B[j][0]*D[T[1][1]][T[2][0]] + B[j][1]*D[T[1][1]][T[2][1]] + B[j][2]*D[T[1][1]][T[2][2]]);
						BKB[1][2] += B[i][2]*(B[j][0]*D[T[2][1]][T[2][0]] + B[j][1]*D[T[2][1]][T[2][1]] + B[j][2]*D[T[2][1]][T[2][2]]);

						BKB[2][0]  = B[i][0]*(B[j][0]*D[T[0][2]][T[0][0]] + B[j][1]*D[T[0][2]][T[0][1]] + B[j][2]*D[T[0][2]][T[0][2]]);
						BKB[2][0] += B[i][1]*(B[j][0]*D[T[1][2]][T[0][0]] + B[j][1]*D[T[1][2]][T[0][1]] + B[j][2]*D[T[1][2]][T[0][2]]);
						BKB[2][0] += B[i][2]*(B[j][0]*D[T[2][2]][T[0][0]] + B[j][1]*D[T[2][2]][T[0][1]] + B[j][2]*D[T[2][2]][T[0][2]]);

						BKB[2][1]  = B[i][0]*(B[j][0]*D[T[0][2]][T[1][0]] + B[j][1]*D[T[0][2]][T[1][1]] + B[j][2]*D[T[0][2]][T[1][2]]); 
						BKB[2][1] += B[i][1]*(B[j][0]*D[T[1][2]][T[1][0]] + B[j][1]*D[T[1][2]][T[1][1]] + B[j][2]*D[T[1][2]][T[1][2]]); 
						BKB[2][1] += B[i][2]*(B[j][0]*D[T[2][2]][T[1][0]] + B[j][1]*D[T[2][2]][T[1][1]] + B[j][2]*D[T[2][2]][T[1][2]]); 

						BKB[2][2]  = B[i][0]*(B[j][0]*D[T[0][2]][T[2][0]] + B[j][1]*D[T[0][2]][T[2][1]] + B[j][2]*D[T[0][2]][T[2][2]]);
						BKB[2][2] += B[i][1]*(B[j][0]*D[T[1][2]][T[2][0]] + B[j][1]*D[T[1][2]][T[2][1]] + B[j][2]*D[T[1][2]][T[2][2]]);
						BKB[2][2] += B[i][2]*(B[j][0]*D[T[2][2]][T[2][0]] + B[j][1]*D[T[2][2]][T[2][1]] + B[j][2]*D[T[2][2]][T[2][2]]);
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
