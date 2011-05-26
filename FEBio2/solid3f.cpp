#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FE3FieldElasticSolidDomain::ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
	// calculate material stiffness (i.e. constitutive component)
	MaterialStiffness(fem, el, ke);

	// calculate geometrical stiffness
	GeometricalStiffness(el, ke);

	// Calculate dilatational stiffness
	DilatationalStiffness(fem, el, ke);

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
//! calculates dilatational element stiffness component for element iel

void FE3FieldElasticSolidDomain::DilatationalStiffness(FEM& fem, FESolidElement& elem, matrix& ke)
{
	int i, j, n;

	const int nint = elem.GaussPoints();
	const int neln = elem.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEElasticMaterial* pm = fem.GetElasticMaterial(elem.GetMatID());

	FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(pm);
	assert(pmi);

	// average global derivatives
	vector<double> gradN(3*neln);
	zero(gradN);

	// initial element volume
	double Ve = 0;

	// global derivatives of shape functions
	double Gx, Gy, Gz;
	const double *gw = elem.GaussWeights();

	// jacobian
	double Ji[3][3], Jt, J0;

	double *Gr, *Gs, *Gt;

	// repeat over gauss-points
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		J0 = detJ0(elem, n);
		Jt = invjact(elem, Ji, n);

		Jt *= gw[n];

		Ve += J0*gw[n];

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

			gradN[3*i  ] += Gx*Jt;
			gradN[3*i+1] += Gy*Jt;
			gradN[3*i+2] += Gz*Jt;
		}
	}

	// get effective modulus
	double k = pmi->UJJ(elem.m_eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
//	k += elem.m_Lk*pmi->hpp(elem.m_eJ);

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
//! Calculates element material stiffness element matrix

void FE3FieldElasticSolidDomain::MaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
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

	// get the material
	FEUncoupledMaterial& mat = dynamic_cast<FEUncoupledMaterial&>(*m_pMat);

	// we need the following tensors for the dilational stiffness
	mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds Cp = IxI - I4*2;

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n)*gw[n];

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;

		// get the material's tangent
		// Note that we are only grabbing the deviatoric tangent. 
		// The other tangent terms depend on the pressure p
		// which we seperately
		tens4ds C = Cp*el.m_ep + mat.DevTangent(mp);

		// get the 'D' matrix
		C.extract(D);

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
//! calculates element's geometrical stiffness component for each integration point

void FE3FieldElasticSolidDomain::GeometricalStiffness(FESolidElement &el, matrix &ke)
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
		detJt = invjact(el, Ji, n)*gw[n];

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
//! This function updates the stresses for elements using the three-field formulation.
//! For such elements, the stress is a sum of a deviatoric stress, calculate by the
//! material and a dilatational term.
void FE3FieldElasticSolidDomain::UpdateStresses(FEM &fem)
{
	int i, n;
	int nint, neln;
	double* gw;

	// get the material
	FEUncoupledMaterial& mat = *(dynamic_cast<FEUncoupledMaterial*>(m_pMat));

	for (i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the number of integration points
		nint = el.GaussPoints();

		// get the integration weights
		gw = el.GaussWeights();

		// number of nodes
		neln = el.Nodes();

		// nodal coordinates
		vec3d r0[8], rt[8];
		for (int j=0; j<neln; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		// calculate the average dilatation and pressure
		double v = 0, V = 0;

		for (n=0; n<nint; ++n)
		{
			v += detJt(el, n)*gw[n];
			V += detJ0(el, n)*gw[n];
		}

		// calculate volume ratio
		el.m_eJ = v / V;

		// Calculate pressure. This is a sum of a Lagrangian term and a penalty term
		//        <----- Lag. mult. ----->   <------ penalty ----->
//		el.m_ep = el.m_Lk*pmi->hp(el.m_eJ) + pmi->Up(el.m_eJ);
		el.m_ep = mat.UJ(el.m_eJ);

		// loop over the integration points and calculate
		// the stress at the integration point
		for (n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

			// material point coordinates
			// TODO: I'm not entirly happy with this solution
			//		 since the material point coordinates are not used by most materials.
			pt.r0 = el.Evaluate(r0, n);
			pt.rt = el.Evaluate(rt, n);

			// get the deformation gradient and determinant
			pt.J = defgrad(el, pt.F, n);

			// three-field element variables
			pt.avgJ = el.m_eJ;
			pt.avgp = el.m_ep;

			// calculate the stress at this material point
			// Note that we don't call the material's Stress member function.
			// The reason is that we need to use the averaged pressure for the element
			// and the Stress function uses the pointwise pressure. 
			// Therefore we call the DevStress function and add the pressure term
			// seperately. 
			pt.s = mat3dd(el.m_ep) + mat.DevStress(mp);
		}
	}
}
