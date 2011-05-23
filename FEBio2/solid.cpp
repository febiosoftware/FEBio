#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "log.h"
#include "FEPoroElastic.h"
#include "FEMicroMaterial.h"
#include "FETransverselyIsotropic.h"
#include "FESolidSolver.h"
#include "FEPointBodyForce.h"

//-----------------------------------------------------------------------------
// FEElasticSolidDomain
//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
}

//-----------------------------------------------------------------------------
bool FEElasticSolidDomain::Initialize(FEModel &mdl)
{
	// initialize base class
	FESolidDomain::Initialize(mdl);

	FEM& fem = dynamic_cast<FEM&>(mdl);

	bool bmerr = false;

	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		// unpack element data
		FESolidElement& el = m_Elem[i];

		try
		{
			UnpackElement(el);
		}
		catch (NegativeJacobian e)
		{
			clog.printbox("F A T A L   E R R O R", "A negative jacobian was detected at\n integration point %d of element %d.\nDid you use the right node numbering?", e.m_iel, e.m_ng);
			return false;
		}

		if (dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver))
		{
			// get the elements material
			FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());

			// set the local element coordinates
			if (pme)
			{
				if (pme->m_pmap)
				{
					for (int n=0; n<el.GaussPoints(); ++n)
					{
						FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
						pt.Q = pme->m_pmap->LocalElementCoord(el, n);
					}
				}
				else
				{
					if (fem.GetDebugFlag())
					{
						// If we get here, then the element has a user-defined fiber axis
						// we should check to see if it has indeed been specified.
						// TODO: This assumes that pt.Q will not get intialized to
						//		 a valid value. I should find another way for checking since I
						//		 would like pt.Q always to be initialized to a decent value.
						if (dynamic_cast<FETransverselyIsotropic*>(pme))
						{
							FEElasticMaterialPoint& pt = *el.m_State[0]->ExtractData<FEElasticMaterialPoint>();
							mat3d& m = pt.Q;
							if (fabs(m.det() - 1) > 1e-7)
							{
								// this element did not get specified a user-defined fiber direction
								clog.printbox("ERROR", "Solid element %d was not assigned a fiber direction.", i+1);
								bmerr = true;
							}
						}
					}
				}
			}
		}
	}

	return (bmerr == false);
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.m_State[j]->Init(false);
	}
}


//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Serialize(DumpFile &ar)
{
	FEM& fem = dynamic_cast<FEM&>(*ar.GetFEM());
	if (ar.IsSaving())
	{
		ar << m_Node;

		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FESolidElement& el = m_Elem[i];
			int nmat = el.GetMatID();
			ar << el.Type();
			
			ar << nmat;
			ar << el.m_nrigid;
			ar << el.m_nID;
			ar << el.m_node;

			ar << el.m_eJ;
			ar << el.m_ep;
			ar << el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}
	}
	else
	{
		ar >> m_Node;

		int n, mat;
		for (size_t i=0; i<m_Elem.size(); ++i)
		{
			FESolidElement& el = m_Elem[i];
			ar >> n;

			el.SetType(n);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nrigid;
			ar >> el.m_nID;
			ar >> el.m_node;

			ar >> el.m_eJ;
			ar >> el.m_ep;
			ar >> el.m_Lk;

			for (int j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(fem.GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// FE3FieldElasticSolidDomain
//-----------------------------------------------------------------------------
bool FE3FieldElasticSolidDomain::Initialize(FEModel &fem)
{
	// make sure the domain material uses an uncoupled formulation
	if (dynamic_cast<FEUncoupledMaterial*>(m_pMat) == 0) return false;
	return FEElasticSolidDomain::Initialize(fem);
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain::Residual(FESolidSolver *psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	// element force vector
	vector<double> fe;

	vector<int> lm;

	int NE = m_Elem.size();
	for (int i=0; i<NE; ++i)
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
		InternalForces(el, fe);

		// apply body forces
		if (fem.HasBodyForces()) BodyForces(fem, el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		psolver->AssembleResidual(el.m_node, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEElasticSolidDomain::InternalForces(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;
	mat3ds s;

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
		s = pt.s;

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

void FEElasticSolidDomain::BodyForces(FEM& fem, FESolidElement& el, vector<double>& fe)
{
	int NF = fem.BodyForces();
	for (int nf = 0; nf < NF; ++nf)
	{
		FEBodyForce& BF = *fem.GetBodyForce(nf);

		// calculate force scale values
		// the "-" sign is to be consistent with NIKE3D's convention
		vec3d g;
		if (BF.lc[0] >= 0) g.x = -fem.GetLoadCurve(BF.lc[0])->Value()*BF.s[0];
		if (BF.lc[1] >= 0) g.y = -fem.GetLoadCurve(BF.lc[1])->Value()*BF.s[1];
		if (BF.lc[2] >= 0) g.z = -fem.GetLoadCurve(BF.lc[2])->Value()*BF.s[2];

		// don't forget to multiply with the density
		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));
		double dens = pme->Density();
		g *= dens;

		// TODO: I don't like this but for now I'll hard-code the modification of the
		//       force center position
		if (dynamic_cast<FEPointBodyForce*>(&BF))
		{
			FEPointBodyForce* pf = dynamic_cast<FEPointBodyForce*>(&BF);
			if (pf->m_rlc[0] >= 0) pf->m_rc.x = fem.GetLoadCurve(pf->m_rlc[0])->Value();
			if (pf->m_rlc[1] >= 0) pf->m_rc.y = fem.GetLoadCurve(pf->m_rlc[1])->Value();
			if (pf->m_rlc[2] >= 0) pf->m_rc.z = fem.GetLoadCurve(pf->m_rlc[2])->Value();
		}

		// jacobian
		double detJ;
		double *H;
		double* gw = el.GaussWeights();
		vec3d f;

		// number of nodes
		int neln = el.Nodes();

		// nodal coordinates
		vec3d r0[8], rt[8];
		for (int i=0; i<neln; ++i)
		{
			r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
			rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
		}

		// loop over integration points
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			pt.r0 = el.Evaluate(r0, n);
			pt.rt = el.Evaluate(rt, n);

			detJ = el.detJ0(n)*gw[n];

			// get the force
			f = BF.force(mp);

			H = el.H(n);

			for (int i=0; i<neln; ++i)
			{
				fe[3*i  ] += H[i]*g.x*f.x*detJ;
				fe[3*i+1] += H[i]*g.y*f.y*detJ;
				fe[3*i+2] += H[i]*g.z*f.z*detJ;
			}						
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEElasticSolidDomain::BodyForceStiffness(FEM& fem, FESolidElement &el, matrix &ke)
{
	int neln = el.Nodes();
	int ndof = ke.columns()/neln;

	int NF = fem.BodyForces();
	for (int nf =0; nf < NF; ++nf)
	{
		FEBodyForce& BF = *fem.GetBodyForce(nf);

		// calculate force scale values
		// the "-" sign is to be consistent with NIKE3D's convention
		vec3d g;
		if (BF.lc[0] >= 0) g.x = -fem.GetLoadCurve(BF.lc[0])->Value()*BF.s[0];
		if (BF.lc[1] >= 0) g.y = -fem.GetLoadCurve(BF.lc[1])->Value()*BF.s[1];
		if (BF.lc[2] >= 0) g.z = -fem.GetLoadCurve(BF.lc[2])->Value()*BF.s[2];

		// don't forget to multiply with the density
		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));
		double dens = pme->Density();
		g *= dens;

		// jacobian
		double detJ;
		double *H;
		double* gw = el.GaussWeights();
		mat3ds K;

		// loop over integration points
		int nint = el.GaussPoints();
		int neln = el.Nodes();
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.m_State[n];
			detJ = el.detJ0(n)*gw[n];

			// get the stiffness
			K = BF.stiffness(mp);

			H = el.H(n);

			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					ke[ndof*i  ][ndof*j  ] += H[i]*H[j]*g.x*K(0,0)*detJ;
					ke[ndof*i  ][ndof*j+1] += H[i]*H[j]*g.x*K(0,1)*detJ;
					ke[ndof*i  ][ndof*j+2] += H[i]*H[j]*g.x*K(0,2)*detJ;

					ke[ndof*i+1][ndof*j  ] += H[i]*H[j]*g.y*K(1,0)*detJ;
					ke[ndof*i+1][ndof*j+1] += H[i]*H[j]*g.y*K(1,1)*detJ;
					ke[ndof*i+1][ndof*j+2] += H[i]*H[j]*g.y*K(1,2)*detJ;

					ke[ndof*i+2][ndof*j  ] += H[i]*H[j]*g.z*K(2,0)*detJ;
					ke[ndof*i+2][ndof*j+1] += H[i]*H[j]*g.z*K(2,1)*detJ;
					ke[ndof*i+2][ndof*j+2] += H[i]*H[j]*g.z*K(2,2)*detJ;
				}
		}	
	}
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n

void FEElasticSolidDomain::GeometricalStiffness(FESolidElement &el, matrix &ke)
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

void FEElasticSolidDomain::MaterialStiffness(FEM& fem, FESolidElement &el, matrix &ke)
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
		pt.avgJ = el.m_eJ;
		pt.avgp = el.m_ep;

		// get the 'D' matrix
		tens4ds C = pmat->Tangent(mp);
		C.extract(D);

		if (dynamic_cast<FEMicroMaterial*>(pmat))
		{
			// the micro-material screws up the currently unpacked elements
			// so I have to unpack the element data again
			UnpackElement(el);
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

void FEElasticSolidDomain::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;

	// element stiffness matrix
	matrix ke;

	vector<int> lm;

	// repeat over all solid elements
	int NE = m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// this element should not be rigid
		assert(!el.IsRigid());

		UnpackElement(el);

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.Create(ndof, ndof);
		ke.zero();

		// calculate the element stiffness matrix
		ElementStiffness(fem, el, ke);

		// add the inertial stiffness for dynamics
		if (fem.m_pStep->m_nanalysis == FE_DYNAMIC) ElementInertialStiffness(fem, el, ke);

		// add body force stiffness
		if (fem.HasBodyForces()) BodyForceStiffness(fem, el, ke);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the element stiffness matrix. It calls the material
//! stiffness function, the geometrical stiffness function and, if necessary, the
//! dilatational stiffness function. Note that these three functions only calculate
//! the upper diagonal matrix due to the symmetry of the element stiffness matrix
//! The last section of this function fills the rest of the element stiffness matrix.

void FEElasticSolidDomain::ElementStiffness(FEM& fem, FESolidElement& el, matrix& ke)
{
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
//! calculates element inertial stiffness matrix
void FEElasticSolidDomain::ElementInertialStiffness(FEM& fem, FESolidElement& el, matrix& ke)
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
void FEElasticSolidDomain::UpdateStresses(FEM &fem)
{
	int i, n;
	int nint, neln;
	double* gw;

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

		// nodall coordinates
		vec3d r0[8], rt[8];
		for (int j=0; j<neln; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		// get the integration weights
		gw = el.GaussWeights();

		// get the material
		FESolidMaterial* pm = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(el.GetMatID()));

		// extract the elastic component
		FEElasticMaterial* pme = fem.GetElasticMaterial(el.GetMatID());

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
			pt.J = defgrad(el, pt.F, n);

			// three-field element variables
			pt.avgJ = el.m_eJ;
			pt.avgp = el.m_ep;

			// calculate the stress at this material point
			pt.s = pm->Stress(mp);

			if (dynamic_cast<FEMicroMaterial*>(pme))
			{
				// the micro-material screws up the currently unpacked elements
				// so I have to unpack the element data again
				UnpackElement(el);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEElasticSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*MAX_NDOFS);
	
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		int* id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[0];
		lm[3*i+1] = id[1];
		lm[3*i+2] = id[2];

		// now the pressure dofs
		lm[3*N+i] = id[6];

		// rigid rotational dofs
		lm[4*N + 3*i  ] = id[7];
		lm[4*N + 3*i+1] = id[8];
		lm[4*N + 3*i+2] = id[9];

		// fill the rest with -1
		lm[7*N + 3*i  ] = -1;
		lm[7*N + 3*i+1] = -1;
		lm[7*N + 3*i+2] = -1;

		lm[10*N + i] = id[10];
		
		// concentration dofs
		lm[11*N + i] = id[11];
	}
}

//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure

void FEElasticSolidDomain::UnpackElement(FEElement& el, unsigned int nflag)
{
	// copy initial nodal coordinates
	// (also needed for traits data)
	if ((nflag & FE_UNPACK_R0) || (nflag & FE_UNPACK_TRAITS))
	{
		int N = el.Nodes();
		vec3d* r0 = el.r0();
		for (int i=0; i<N; ++i)
		{
			int n = el.m_node[i];

			FENode& node = m_pMesh->Node(n);

			r0[i] = node.m_r0;
		}
	}

	// copy current nodal coordinates
	// (also needed for traits data)
	if ((nflag & FE_UNPACK_RT) || (nflag & FE_UNPACK_TRAITS))
	{
		int N = el.Nodes();
		vec3d* rt = el.rt();
		for (int i=0; i<N; ++i)
		{
			int n = el.m_node[i];

			FENode& node = m_pMesh->Node(n);

			rt[i] = node.m_rt;
		}
	}

	// unpack the traits data
	if (nflag & FE_UNPACK_TRAITS) el.UnpackTraitsData(nflag);
}
