#include "stdafx.h"
#include "FESolver.h"
#include "FEPoroElastic.h"
#include "FEViscoElasticMaterial.h"
#include <math.h>

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements
//! Note that we use a one-point gauss integration rule for the thickness
//! integration. This will integrate linear functions exactly.

void FESolver::InternalForces(FEShellElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	const double* Gr, *Gs, *H;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();
	double gt;

	double* h0 = el.m_h0;
	double za;

	double Nx, Ny, Nz;
	double Mx, My, Mz;

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEElasticMaterialPoint& pt = *(el.m_State[n]->ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.s;

		gt = el.gt(n);

		Gr = el.Hr(n);
		Gs = el.Hs(n);
		H  = el.H(n);

		for (i=0; i<neln; ++i)
		{
			za = 0.5*gt*h0[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Nx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i];
			Ny = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i];
			Nz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i];

			Mx = za*Ji[0][0]*Gr[i] + za*Ji[1][0]*Gs[i] + Ji[2][0]*0.5*h0[i]*H[i];
			My = za*Ji[0][1]*Gr[i] + za*Ji[1][1]*Gs[i] + Ji[2][1]*0.5*h0[i]*H[i];
			Mz = za*Ji[0][2]*Gr[i] + za*Ji[1][2]*Gs[i] + Ji[2][2]*0.5*h0[i]*H[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[6*i  ] -= ( Nx*s.xx()  +
				           Ny*s.xy() +
					       Nz*s.xz() )*detJt;

			fe[6*i+1] -= ( Ny*s.yy()  +
				           Nx*s.xy() +
					       Nz*s.yz() )*detJt;

			fe[6*i+2] -= ( Nz*s.zz()  +
				           Ny*s.yz() +
					       Nx*s.xz() )*detJt;

			fe[6*i+3] -= ( Mx*s.xx()  +
				           My*s.xy() +
					       Mz*s.xz() )*detJt;

			fe[6*i+4] -= ( My*s.yy()  +
				           Mx*s.xy() +
					       Mz*s.yz() )*detJt;

			fe[6*i+5] -= ( Mz*s.zz()  +
				           My*s.yz() +
					       Mx*s.xz() )*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the shell element stiffness matrix

void FESolver::ElementStiffness(FEShellElement& el, matrix& ke)
{
	int i, i6, j, j6, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 6*neln;

	// stiffness components for the initial stress component of stiffness matrix
	double kab;

	// global derivatives of shape functions
	// NOTE: hard-coding of quad elements!
	// Gx = dH/dx
	double Nx[4], Ny[4], Nz[4];
	double Mx[4], My[4], Mz[4];

	double Nxi, Nyi, Nzi;
	double Nxj, Nyj, Nzj;
	double Mxi, Myi, Mzi;
	double Mxj, Myj, Mzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][6];

	// element stress
	mat3ds s;

	// get the element's material
	FEMaterial* pm = m_fem.GetMaterial(el.GetMatID());

	// extract the elastic component
	FEElasticMaterial* pme = m_fem.GetElasticMaterial(el.GetMatID());

	double *Grn, *Gsn, *Hn;
	double Gr, Gs, H;

	// jacobian
	double Ji[3][3], detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate the average thickness
	double* h0 = el.m_h0, gt, za;

	// see if this is a poroelastic material
	bool bporo = false;
	if ((m_fem.m_pStep->m_itype == FE_STATIC_PORO) && (dynamic_cast<FEPoroElastic*>(pm))) bporo = true;

	// calculate element stiffness matrix
	ke.zero();
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *(el.m_State[n]);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate jacobian
		el.invjact(Ji, n);
		detJt = el.detJt(n)*gw[n];

		Grn = el.Hr(n);
		Gsn = el.Hs(n);
		Hn  = el.H(n);

		gt = el.gt(n);

		// ------------ constitutive component --------------

		// setup the material point
//		el.defgrad(pt.F, n);
//		pt.J = el.detF(n);

		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;
//		pt.Q = el.m_Q[n];


		if (bporo)
		{
			FEPoroElasticMaterialPoint& pt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());

			pt.m_p = 0;

			// evaluate fluid pressure at gauss-point
//			pt.p = el.Evaluate(el.pt(), n);
		}

		// get the 'D' matrix
		pm->Tangent(D, mp);

		for (i=0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			H  = Hn[i];

			za = 0.5*gt*h0[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Nx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs;
			Ny[i] = Ji[0][1]*Gr+Ji[1][1]*Gs;
			Nz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs;

			Mx[i] = za*Ji[0][0]*Gr + za*Ji[1][0]*Gs + Ji[2][0]*0.5*h0[i]*H;
			My[i] = za*Ji[0][1]*Gr + za*Ji[1][1]*Gs + Ji[2][1]*0.5*h0[i]*H;
			Mz[i] = za*Ji[0][2]*Gr + za*Ji[1][2]*Gs + Ji[2][2]*0.5*h0[i]*H;
		}

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (i=0, i6=0; i<neln; ++i, i6 += 6)
		{
			Nxi = Nx[i];
			Nyi = Ny[i];
			Nzi = Nz[i];

			Mxi = Mx[i];
			Myi = My[i];
			Mzi = Mz[i];

			for (j=i, j6 = i6; j<neln; ++j, j6 += 6)
			{
				Nxj = Nx[j];
				Nyj = Ny[j];
				Nzj = Nz[j];

				Mxj = Mx[j];
				Myj = My[j];
				Mzj = Mz[j];


				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Nxj+D[0][3]*Nyj+D[0][5]*Nzj);
				DBL[0][1] = (D[0][1]*Nyj+D[0][3]*Nxj+D[0][4]*Nzj);
				DBL[0][2] = (D[0][2]*Nzj+D[0][4]*Nyj+D[0][5]*Nxj);
				DBL[0][3] = (D[0][0]*Mxj+D[0][3]*Myj+D[0][5]*Mzj);
				DBL[0][4] = (D[0][1]*Myj+D[0][3]*Mxj+D[0][4]*Mzj);
				DBL[0][5] = (D[0][2]*Mzj+D[0][4]*Myj+D[0][5]*Mxj);

				DBL[1][0] = (D[1][0]*Nxj+D[1][3]*Nyj+D[1][5]*Nzj);
				DBL[1][1] = (D[1][1]*Nyj+D[1][3]*Nxj+D[1][4]*Nzj);
				DBL[1][2] = (D[1][2]*Nzj+D[1][4]*Nyj+D[1][5]*Nxj);
				DBL[1][3] = (D[1][0]*Mxj+D[1][3]*Myj+D[1][5]*Mzj);
				DBL[1][4] = (D[1][1]*Myj+D[1][3]*Mxj+D[1][4]*Mzj);
				DBL[1][5] = (D[1][2]*Mzj+D[1][4]*Myj+D[1][5]*Mxj);

				DBL[2][0] = (D[2][0]*Nxj+D[2][3]*Nyj+D[2][5]*Nzj);
				DBL[2][1] = (D[2][1]*Nyj+D[2][3]*Nxj+D[2][4]*Nzj);
				DBL[2][2] = (D[2][2]*Nzj+D[2][4]*Nyj+D[2][5]*Nxj);
				DBL[2][3] = (D[2][0]*Mxj+D[2][3]*Myj+D[2][5]*Mzj);
				DBL[2][4] = (D[2][1]*Myj+D[2][3]*Mxj+D[2][4]*Mzj);
				DBL[2][5] = (D[2][2]*Mzj+D[2][4]*Myj+D[2][5]*Mxj);

				DBL[3][0] = (D[3][0]*Nxj+D[3][3]*Nyj+D[3][5]*Nzj);
				DBL[3][1] = (D[3][1]*Nyj+D[3][3]*Nxj+D[3][4]*Nzj);
				DBL[3][2] = (D[3][2]*Nzj+D[3][4]*Nyj+D[3][5]*Nxj);
				DBL[3][3] = (D[3][0]*Mxj+D[3][3]*Myj+D[3][5]*Mzj);
				DBL[3][4] = (D[3][1]*Myj+D[3][3]*Mxj+D[3][4]*Mzj);
				DBL[3][5] = (D[3][2]*Mzj+D[3][4]*Myj+D[3][5]*Mxj);

				DBL[4][0] = (D[4][0]*Nxj+D[4][3]*Nyj+D[4][5]*Nzj);
				DBL[4][1] = (D[4][1]*Nyj+D[4][3]*Nxj+D[4][4]*Nzj);
				DBL[4][2] = (D[4][2]*Nzj+D[4][4]*Nyj+D[4][5]*Nxj);
				DBL[4][3] = (D[4][0]*Mxj+D[4][3]*Myj+D[4][5]*Mzj);
				DBL[4][4] = (D[4][1]*Myj+D[4][3]*Mxj+D[4][4]*Mzj);
				DBL[4][5] = (D[4][2]*Mzj+D[4][4]*Myj+D[4][5]*Mxj);

				DBL[5][0] = (D[5][0]*Nxj+D[5][3]*Nyj+D[5][5]*Nzj);
				DBL[5][1] = (D[5][1]*Nyj+D[5][3]*Nxj+D[5][4]*Nzj);
				DBL[5][2] = (D[5][2]*Nzj+D[5][4]*Nyj+D[5][5]*Nxj);
				DBL[5][3] = (D[5][0]*Mxj+D[5][3]*Myj+D[5][5]*Mzj);
				DBL[5][4] = (D[5][1]*Myj+D[5][3]*Mxj+D[5][4]*Mzj);
				DBL[5][5] = (D[5][2]*Mzj+D[5][4]*Myj+D[5][5]*Mxj);

				ke[i6  ][j6  ] += (Nxi*DBL[0][0] + Nyi*DBL[3][0] + Nzi*DBL[5][0] )*detJt;
				ke[i6  ][j6+1] += (Nxi*DBL[0][1] + Nyi*DBL[3][1] + Nzi*DBL[5][1] )*detJt;
				ke[i6  ][j6+2] += (Nxi*DBL[0][2] + Nyi*DBL[3][2] + Nzi*DBL[5][2] )*detJt;
				ke[i6  ][j6+3] += (Nxi*DBL[0][3] + Nyi*DBL[3][3] + Nzi*DBL[5][3] )*detJt;
				ke[i6  ][j6+4] += (Nxi*DBL[0][4] + Nyi*DBL[3][4] + Nzi*DBL[5][4] )*detJt;
				ke[i6  ][j6+5] += (Nxi*DBL[0][5] + Nyi*DBL[3][5] + Nzi*DBL[5][5] )*detJt;

				ke[i6+1][j6  ] += (Nyi*DBL[1][0] + Nxi*DBL[3][0] + Nzi*DBL[4][0] )*detJt;
				ke[i6+1][j6+1] += (Nyi*DBL[1][1] + Nxi*DBL[3][1] + Nzi*DBL[4][1] )*detJt;
				ke[i6+1][j6+2] += (Nyi*DBL[1][2] + Nxi*DBL[3][2] + Nzi*DBL[4][2] )*detJt;
				ke[i6+1][j6+3] += (Nyi*DBL[1][3] + Nxi*DBL[3][3] + Nzi*DBL[4][3] )*detJt;
				ke[i6+1][j6+4] += (Nyi*DBL[1][4] + Nxi*DBL[3][4] + Nzi*DBL[4][4] )*detJt;
				ke[i6+1][j6+5] += (Nyi*DBL[1][5] + Nxi*DBL[3][5] + Nzi*DBL[4][5] )*detJt;

				ke[i6+2][j6  ] += (Nzi*DBL[2][0] + Nyi*DBL[4][0] + Nxi*DBL[5][0] )*detJt;
				ke[i6+2][j6+1] += (Nzi*DBL[2][1] + Nyi*DBL[4][1] + Nxi*DBL[5][1] )*detJt;
				ke[i6+2][j6+2] += (Nzi*DBL[2][2] + Nyi*DBL[4][2] + Nxi*DBL[5][2] )*detJt;
				ke[i6+2][j6+3] += (Nzi*DBL[2][3] + Nyi*DBL[4][3] + Nxi*DBL[5][3] )*detJt;
				ke[i6+2][j6+4] += (Nzi*DBL[2][4] + Nyi*DBL[4][4] + Nxi*DBL[5][4] )*detJt;
				ke[i6+2][j6+5] += (Nzi*DBL[2][5] + Nyi*DBL[4][5] + Nxi*DBL[5][5] )*detJt;

				ke[i6+3][j6  ] += (Mxi*DBL[0][0] + Myi*DBL[3][0] + Mzi*DBL[5][0] )*detJt;
				ke[i6+3][j6+1] += (Mxi*DBL[0][1] + Myi*DBL[3][1] + Mzi*DBL[5][1] )*detJt;
				ke[i6+3][j6+2] += (Mxi*DBL[0][2] + Myi*DBL[3][2] + Mzi*DBL[5][2] )*detJt;
				ke[i6+3][j6+3] += (Mxi*DBL[0][3] + Myi*DBL[3][3] + Mzi*DBL[5][3] )*detJt;
				ke[i6+3][j6+4] += (Mxi*DBL[0][4] + Myi*DBL[3][4] + Mzi*DBL[5][4] )*detJt;
				ke[i6+3][j6+5] += (Mxi*DBL[0][5] + Myi*DBL[3][5] + Mzi*DBL[5][5] )*detJt;

				ke[i6+4][j6  ] += (Myi*DBL[1][0] + Mxi*DBL[3][0] + Mzi*DBL[4][0] )*detJt;
				ke[i6+4][j6+1] += (Myi*DBL[1][1] + Mxi*DBL[3][1] + Mzi*DBL[4][1] )*detJt;
				ke[i6+4][j6+2] += (Myi*DBL[1][2] + Mxi*DBL[3][2] + Mzi*DBL[4][2] )*detJt;
				ke[i6+4][j6+3] += (Myi*DBL[1][3] + Mxi*DBL[3][3] + Mzi*DBL[4][3] )*detJt;
				ke[i6+4][j6+4] += (Myi*DBL[1][4] + Mxi*DBL[3][4] + Mzi*DBL[4][4] )*detJt;
				ke[i6+4][j6+5] += (Myi*DBL[1][5] + Mxi*DBL[3][5] + Mzi*DBL[4][5] )*detJt;

				ke[i6+5][j6  ] += (Mzi*DBL[2][0] + Myi*DBL[4][0] + Mxi*DBL[5][0] )*detJt;
				ke[i6+5][j6+1] += (Mzi*DBL[2][1] + Myi*DBL[4][1] + Mxi*DBL[5][1] )*detJt;
				ke[i6+5][j6+2] += (Mzi*DBL[2][2] + Myi*DBL[4][2] + Mxi*DBL[5][2] )*detJt;
				ke[i6+5][j6+3] += (Mzi*DBL[2][3] + Myi*DBL[4][3] + Mxi*DBL[5][3] )*detJt;
				ke[i6+5][j6+4] += (Mzi*DBL[2][4] + Myi*DBL[4][4] + Mxi*DBL[5][4] )*detJt;
				ke[i6+5][j6+5] += (Mzi*DBL[2][5] + Myi*DBL[4][5] + Mxi*DBL[5][5] )*detJt;
			}
		}

		// ------------ initial stress component --------------
	
		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		s = pt.s;

		for (i=0; i<neln; ++i)
			for (j=i; j<neln; ++j)
			{
				// the v-v component
				kab = (Nx[i]*(s.xx()*Nx[j]+s.xy()*Ny[j]+s.xz()*Nz[j]) +
					   Ny[i]*(s.xy()*Nx[j]+s.yy()*Ny[j]+s.yz()*Nz[j]) + 
					   Nz[i]*(s.xz()*Nx[j]+s.yz()*Ny[j]+s.zz()*Nz[j]))*detJt;

				ke[6*i  ][6*j  ] += kab;
				ke[6*i+1][6*j+1] += kab;
				ke[6*i+2][6*j+2] += kab;

				// the v-t component
				kab = (Nx[i]*(s.xx()*Mx[j]+s.xy()*My[j]+s.xz()*Mz[j]) +
					   Ny[i]*(s.xy()*Mx[j]+s.yy()*My[j]+s.yz()*Mz[j]) + 
					   Nz[i]*(s.xz()*Mx[j]+s.yz()*My[j]+s.zz()*Mz[j]))*detJt;

				ke[6*i  ][6*j+3] += kab;
				ke[6*i+1][6*j+4] += kab;
				ke[6*i+2][6*j+5] += kab;

				// the t-v component
				kab = (Mx[i]*(s.xx()*Nx[j]+s.xy()*Ny[j]+s.xz()*Nz[j]) +
					   My[i]*(s.xy()*Nx[j]+s.yy()*Ny[j]+s.yz()*Nz[j]) + 
					   Mz[i]*(s.xz()*Nx[j]+s.yz()*Ny[j]+s.zz()*Nz[j]))*detJt;

				ke[6*i+3][6*j  ] += kab;
				ke[6*i+4][6*j+1] += kab;
				ke[6*i+5][6*j+2] += kab;

				// the t-t component
				kab = (Mx[i]*(s.xx()*Mx[j]+s.xy()*My[j]+s.xz()*Mz[j]) +
					   My[i]*(s.xy()*Mx[j]+s.yy()*My[j]+s.yz()*Mz[j]) + 
					   Mz[i]*(s.xz()*Mx[j]+s.yz()*My[j]+s.zz()*Mz[j]))*detJt;

				ke[6*i+3][6*j+3] += kab;
				ke[6*i+4][6*j+4] += kab;
				ke[6*i+5][6*j+5] += kab;
			}

	} // end loop over gauss-points

	// Dilatational stiffness component
	// Only for (nearly) incompressible materials
	if (dynamic_cast<FEIncompressibleMaterial*>(pme)) DilatationalStiffness(el, ke);

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}

//-----------------------------------------------------------------------------
//! Calculates body forces for shells

void FESolver::BodyForces(FEShellElement& el, vector<double>& fe)
{
	// get the element's material
	FEMaterial* pme = m_fem.GetMaterial(el.GetMatID());

	// material density
	double dens = pme->Density();

	// "gravity" force
	vec3d g = m_fem.m_acc*dens;

	// integration weights
	double* gw = el.GaussWeights();

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double *Hn, detJ;

	// calculate the average thickness
	double* h0 = el.m_h0, gt, za;

	// loop over integration points
	for (int n=0; n<nint; ++n)
	{
		detJ = el.detJ0(n)*gw[n];
		Hn  = el.H(n);
		gt = el.gt(n);

		for (int i=0; i<neln; ++i)
		{
			za = 0.5*gt*h0[i];

			fe[6*i  ] += Hn[i]*g.x*detJ;
			fe[6*i+1] += Hn[i]*g.y*detJ;
			fe[6*i+2] += Hn[i]*g.z*detJ;

			fe[6*i+3] += za*Hn[i]*g.x*detJ;
			fe[6*i+4] += za*Hn[i]*g.y*detJ;
			fe[6*i+5] += za*Hn[i]*g.z*detJ;
		}
	}
}
