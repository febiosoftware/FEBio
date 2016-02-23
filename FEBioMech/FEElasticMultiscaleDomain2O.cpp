#include "stdafx.h"
#include "FEElasticMultiscaleDomain2O.h"
#include "FEMicroMaterial2O.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"

//-----------------------------------------------------------------------------
//! constructor
FEElasticMultiscaleDomain2O::FEElasticMultiscaleDomain2O(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialize element data
bool FEElasticMultiscaleDomain2O::Initialize(FEModel& fem)
{
	if (FEElasticSolidDomain::Initialize(fem) == false) return false;

	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
	FEMesh& m = *GetMesh();
		
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);
	FEModel& rve = pmat->m_mrve;
			
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int i=0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
			xt[i] = m.Node(el.m_node[i]).m_rt;
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
			pt.m_r0 = r0;
			pt.m_rt = rt;

			pt.m_J = defgrad(el, pt.m_F, j);
			defhess(el, mmpt2O.m_G, j);

			if (mmpt2O.m_rve_init == false)
			{
				mmpt2O.m_rve.CopyFrom(rve);
				mmpt2O.m_rve.Init();
				mmpt2O.m_rve_prev.CopyFrom(rve);
				mmpt2O.m_rve_prev.Init();
				mmpt2O.m_rve_init = true;
			}
		}
	}

	// create the probes
	int NP = pmat->Probes();
	for (int i=0; i<NP; ++i)
	{
		FEMicroProbe& p = pmat->Probe(i);
		FEElement* pel = FindElementFromID(p.m_neid);
		if (pel)
		{
			int nint = pel->GaussPoints();
			int ngp = p.m_ngp - 1;
			if ((ngp>=0)&&(ngp<nint))
			{
				FEMaterialPoint& mp = *pel->GetMaterialPoint(ngp);
				FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
				FERVEProbe* prve = new FERVEProbe(fem, mmpt.m_rve, p.m_szfile);
			}
			else return fecore_error("Invalid gausspt number for micro-probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
		}
		else return fecore_error("Invalid Element ID for micro probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::InternalForces(FEGlobalVector& R)
{
	// call base class first
	FEElasticSolidDomain::InternalForces(R);

	// add the flux parts
	InternalWorkFlux(R);
}

//-----------------------------------------------------------------------------
//! Evaluate contribution of discrete-Galerkin enforcement of stress flux.
void FEElasticMultiscaleDomain2O::InternalWorkFlux(FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	// loop over all elements
	int NE = Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = Element(i);
		int neln = el.Nodes();

		// allocate force vector
		int ndof = neln*3;
		fe.resize(ndof, 0.0);

		// evaluate force vector
		InternalElementWorkFlux(el, fe);

		// unpack the LM values
		UnpackLM(el, lm);

		// assemble 
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::InternalElementWorkFlux(FESolidElement& el, vector<double>& fe)
{
	FEMesh& mesh = *GetMesh();

	double Ji[3][3];
	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];

	double Hr[FEElement::MAX_NODES];
	double Hs[FEElement::MAX_NODES];

	int neln = el.Nodes();

	int face[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];

	// get the material
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(GetMaterial());
	assert(pmat);

	// get the integration rule
	FEElementLibrary& ELib = *FEElementLibrary::GetInstance();
	FESurfaceElementTraits& rule = dynamic_cast<FESurfaceElementTraits&>(*ELib.GetElementTraits(FE_TRI6G3));

	// loop over all the element's faces
	int nf = mesh.Faces(el);
	for (int j=0; j<nf; ++j)
	{
		// determine the sign of the contribution based on the neigbor ID
		double sgn = 1.0; // TODO: determine sign

		// get the nodal coordinates of the face
		int nn = mesh.GetFace(el, j, face);
		for (int k=0; k<nn; ++k) rt[k] = mesh.Node(face[k]).m_rt;

		// loop over integration rule
		int nint = rule.nint;
		vector<double>& gw = rule.gw;
		for (int n=0; n<nint; ++n)
		{
			// iso-parametric coordinates in surface element
			double h1 = rule.gr[n];
			double h2 = rule.gs[n];

			// calculate jacobian on surface
			rule.shape_deriv(Hr, Hs, h1, h2);
			vec3d r1(0,0,0), r2(0,0,0);
			for (int k=0; k<nn; ++k)
			{
				r1 += rt[k]*Hr[k];
				r2 += rt[k]*Hs[k];
			}
			vec3d nu = r1^r2;
			double J = nu.norm();

			// calculate iso-parametric coordinates in parent element
			double g1, g2, g3;
			switch (j)
			{
			case 0: g1 = h1; g2 = h2; g3 = 0.; break;
			case 1: g1 = h1; g2 = 0.; g3 = h2; break;
			case 2: g1 = 0.; g2 = h1; g3 = h2; break;
			case 3: g1 = 1.0-h1-h2; g2 = h1; g3 = h2; break;
			}

			// evaluate the tensor at this integration point
			// setup a material point
			FEElasticMaterialPoint pt;
			defgrad(el, pt.m_F, g1, g2, g3);
			pt.m_J = pt.m_F.det();
			mat3ds s = pmat->Stress(pt);

			// evaluate the spatial gradient of shape functions
			invjact(el, Ji, g1, g2, g3);
			el.shape_deriv(Gr, Gs, Gt, g1, g2, g3);
			for (int i=0; i<neln; ++i)
			{
				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !
				vec3d G;
				G.x = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
				G.y = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
				G.z = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

				// put it all together
				vec3d SgradN = s*(G*J);

				fe[3*i  ] += SgradN.x*gw[n]*sgn;
				fe[3*i+1] += SgradN.y*gw[n]*sgn;
				fe[3*i+2] += SgradN.z*gw[n]*sgn;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements
void FEElasticMultiscaleDomain2O::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

	double Gx, Gy, Gz;
	
	mat3ds s; s.zero();
	tens3ds tau; tau.zero();

	const double* Gr, *Gs, *Gt;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());

		// calculate the jacobian
		detJt = invjact(el, Ji, n);
		
		detJt *= gw[n];

		// get the stress vector for this integration point
		s = pt.m_s;
		tau = mmpt2O.m_tau;

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

		//// LTE 
		// Calculate G (Deformation Hessian)
		double *Grrn = el.Grr(n);
		double *Grsn = el.Grs(n);
		double *Grtn = el.Grt(n);

		double *Gsrn = el.Gsr(n);
		double *Gssn = el.Gss(n);
		double *Gstn = el.Gst(n);

		double *Gtrn = el.Gtr(n);
		double *Gtsn = el.Gts(n);
		double *Gttn = el.Gtt(n);


		for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			double Grri = Grrn[i];
			double Grsi = Grsn[i];
			double Grti = Grtn[i];

			double Gsri = Gsrn[i];
			double Gssi = Gssn[i];
			double Gsti = Gstn[i];

			double Gtri = Gtrn[i];
			double Gtsi = Gtsn[i];
			double Gtti = Gttn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !	
			// Calculate based on Gjk = Ji[l][j]*Glm*Ji[m][k]
			double GXX = Ji[0][0]*Grri*Ji[0][0] + Ji[0][0]*Grsi*Ji[1][0] + Ji[0][0]*Grti*Ji[2][0] + Ji[1][0]*Gsri*Ji[0][0] + Ji[1][0]*Gssi*Ji[1][0] + Ji[1][0]*Gsti*Ji[2][0] + Ji[2][0]*Gtri*Ji[0][0] + Ji[2][0]*Gtsi*Ji[1][0] + Ji[2][0]*Gtti*Ji[2][0];
			double GXY = Ji[0][0]*Grri*Ji[0][1] + Ji[0][0]*Grsi*Ji[1][1] + Ji[0][0]*Grti*Ji[2][1] + Ji[1][0]*Gsri*Ji[0][1] + Ji[1][0]*Gssi*Ji[1][1] + Ji[1][0]*Gsti*Ji[2][1] + Ji[2][0]*Gtri*Ji[0][1] + Ji[2][0]*Gtsi*Ji[1][1] + Ji[2][0]*Gtti*Ji[2][1];
			double GXZ = Ji[0][0]*Grri*Ji[0][2] + Ji[0][0]*Grsi*Ji[1][2] + Ji[0][0]*Grti*Ji[2][2] + Ji[1][0]*Gsri*Ji[0][2] + Ji[1][0]*Gssi*Ji[1][2] + Ji[1][0]*Gsti*Ji[2][2] + Ji[2][0]*Gtri*Ji[0][2] + Ji[2][0]*Gtsi*Ji[1][2] + Ji[2][0]*Gtti*Ji[2][2];
		
			double GYX = Ji[0][1]*Grri*Ji[0][0] + Ji[0][1]*Grsi*Ji[1][0] + Ji[0][1]*Grti*Ji[2][0] + Ji[1][1]*Gsri*Ji[0][0] + Ji[1][1]*Gssi*Ji[1][0] + Ji[1][1]*Gsti*Ji[2][0] + Ji[2][1]*Gtri*Ji[0][0] + Ji[2][1]*Gtsi*Ji[1][0] + Ji[2][1]*Gtti*Ji[2][0];
			double GYY = Ji[0][1]*Grri*Ji[0][1] + Ji[0][1]*Grsi*Ji[1][1] + Ji[0][1]*Grti*Ji[2][1] + Ji[1][1]*Gsri*Ji[0][1] + Ji[1][1]*Gssi*Ji[1][1] + Ji[1][1]*Gsti*Ji[2][1] + Ji[2][1]*Gtri*Ji[0][1] + Ji[2][1]*Gtsi*Ji[1][1] + Ji[2][1]*Gtti*Ji[2][1];
			double GYZ = Ji[0][1]*Grri*Ji[0][2] + Ji[0][1]*Grsi*Ji[1][2] + Ji[0][1]*Grti*Ji[2][2] + Ji[1][1]*Gsri*Ji[0][2] + Ji[1][1]*Gssi*Ji[1][2] + Ji[1][1]*Gsti*Ji[2][2] + Ji[2][1]*Gtri*Ji[0][2] + Ji[2][1]*Gtsi*Ji[1][2] + Ji[2][1]*Gtti*Ji[2][2];

			double GZX = Ji[0][2]*Grri*Ji[0][0] + Ji[0][2]*Grsi*Ji[1][0] + Ji[0][2]*Grti*Ji[2][0] + Ji[1][2]*Gsri*Ji[0][0] + Ji[1][2]*Gssi*Ji[1][0] + Ji[1][2]*Gsti*Ji[2][0] + Ji[2][2]*Gtri*Ji[0][0] + Ji[2][2]*Gtsi*Ji[1][0] + Ji[2][2]*Gtti*Ji[2][0];
			double GZY = Ji[0][2]*Grri*Ji[0][1] + Ji[0][2]*Grsi*Ji[1][1] + Ji[0][2]*Grti*Ji[2][1] + Ji[1][2]*Gsri*Ji[0][1] + Ji[1][2]*Gssi*Ji[1][1] + Ji[1][2]*Gsti*Ji[2][1] + Ji[2][2]*Gtri*Ji[0][1] + Ji[2][2]*Gtsi*Ji[1][1] + Ji[2][2]*Gtti*Ji[2][1];
			double GZZ = Ji[0][2]*Grri*Ji[0][2] + Ji[0][2]*Grsi*Ji[1][2] + Ji[0][2]*Grti*Ji[2][2] + Ji[1][2]*Gsri*Ji[0][2] + Ji[1][2]*Gssi*Ji[1][2] + Ji[1][2]*Gsti*Ji[2][2] + Ji[2][2]*Gtri*Ji[0][2] + Ji[2][2]*Gtsi*Ji[1][2] + Ji[2][2]*Gtti*Ji[2][2];
			
			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -=  (Gx*s.xx() + Gy*s.xy() + Gz*s.xz())*detJt +
						   (tau.d[0]*GXX + tau.d[1]*(GXY + GYX) + tau.d[2]*(GXZ + GZX) + tau.d[3]*GYY + tau.d[4]*(GYZ + GZY) + tau.d[5]*GZZ)*detJt*detJt;

			fe[3*i+1] -=  (Gx*s.xy() + Gy*s.yy() + Gz*s.yz())*detJt +
		   				   (tau.d[1]*GXX + tau.d[3]*(GXY + GYX) + tau.d[4]*(GXZ + GZX) + tau.d[6]*GYY + tau.d[7]*(GYZ + GZY) + tau.d[8]*GZZ)*detJt*detJt;

			fe[3*i+2] -=  (Gx*s.xz() + Gy*s.yz() + Gz*s.zz())*detJt +
		   				   (tau.d[2]*GXX + tau.d[4]*(GXY + GYX) + tau.d[5]*(GXZ + GZX) + tau.d[7]*GYY + tau.d[8]*(GYZ + GZY) + tau.d[9]*GZZ)*detJt*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FEElasticMultiscaleDomain2O::UpdateElementStress(int iel, double dt)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];
	
	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
		rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
	}

	// get the integration weights
	double* gw = el.GaussWeights();

	// calculate the stress at this material point
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());
		
		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);

		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
		defhess(el, mmpt2O.m_G, n);

		pmat->Stress2O(mp);
	}
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for integration point n
void FEElasticMultiscaleDomain2O::ElementGeometricalStiffness(FESolidElement &el, matrix &ke)
{
	int n, i, j;

	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];
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

		//// LTE 
		double *Grrn = el.Grr(n);
		double *Grsn = el.Grs(n);
		double *Grtn = el.Grt(n);

		double *Gsrn = el.Gsr(n);
		double *Gssn = el.Gss(n);
		double *Gstn = el.Gst(n);

		double *Gtrn = el.Gtr(n);
		double *Gtsn = el.Gts(n);
		double *Gttn = el.Gtt(n);


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
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());
		
		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds s = pt.m_s;
		tens3ds tau = mmpt2O.m_tau;

		double Grrj, Grsj, Grtj, Gsrj, Gssj, Gstj, Gtrj, Gtsj, Gttj;

		for (i=0; i<neln; ++i)
			
			for (j=i; j<neln; ++j)
			{
				Grrj = Grrn[j];
				Grsj = Grsn[j];
				Grtj = Grtn[j];

				Gsrj = Gsrn[j];
				Gssj = Gssn[j];
				Gstj = Gstn[j];

				Gtrj = Gtrn[j];
				Gtsj = Gtsn[j];
				Gttj = Gttn[j];

				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !	
				// Calculate based on Gjk = Ji[l][j]*Glm*Ji[m][k]
				double GXX = Ji[0][0]*Grrj*Ji[0][0] + Ji[0][0]*Grsj*Ji[1][0] + Ji[0][0]*Grtj*Ji[2][0] + Ji[1][0]*Gsrj*Ji[0][0] + Ji[1][0]*Gssj*Ji[1][0] + Ji[1][0]*Gstj*Ji[2][0] + Ji[2][0]*Gtrj*Ji[0][0] + Ji[2][0]*Gtsj*Ji[1][0] + Ji[2][0]*Gttj*Ji[2][0];
				double GXY = Ji[0][0]*Grrj*Ji[0][1] + Ji[0][0]*Grsj*Ji[1][1] + Ji[0][0]*Grtj*Ji[2][1] + Ji[1][0]*Gsrj*Ji[0][1] + Ji[1][0]*Gssj*Ji[1][1] + Ji[1][0]*Gstj*Ji[2][1] + Ji[2][0]*Gtrj*Ji[0][1] + Ji[2][0]*Gtsj*Ji[1][1] + Ji[2][0]*Gttj*Ji[2][1];
				double GXZ = Ji[0][0]*Grrj*Ji[0][2] + Ji[0][0]*Grsj*Ji[1][2] + Ji[0][0]*Grtj*Ji[2][2] + Ji[1][0]*Gsrj*Ji[0][2] + Ji[1][0]*Gssj*Ji[1][2] + Ji[1][0]*Gstj*Ji[2][2] + Ji[2][0]*Gtrj*Ji[0][2] + Ji[2][0]*Gtsj*Ji[1][2] + Ji[2][0]*Gttj*Ji[2][2];
	
				double GYX = Ji[0][1]*Grrj*Ji[0][0] + Ji[0][1]*Grsj*Ji[1][0] + Ji[0][1]*Grtj*Ji[2][0] + Ji[1][1]*Gsrj*Ji[0][0] + Ji[1][1]*Gssj*Ji[1][0] + Ji[1][1]*Gstj*Ji[2][0] + Ji[2][1]*Gtrj*Ji[0][0] + Ji[2][1]*Gtsj*Ji[1][0] + Ji[2][1]*Gttj*Ji[2][0];
				double GYY = Ji[0][1]*Grrj*Ji[0][1] + Ji[0][1]*Grsj*Ji[1][1] + Ji[0][1]*Grtj*Ji[2][1] + Ji[1][1]*Gsrj*Ji[0][1] + Ji[1][1]*Gssj*Ji[1][1] + Ji[1][1]*Gstj*Ji[2][1] + Ji[2][1]*Gtrj*Ji[0][1] + Ji[2][1]*Gtsj*Ji[1][1] + Ji[2][1]*Gttj*Ji[2][1];
				double GYZ = Ji[0][1]*Grrj*Ji[0][2] + Ji[0][1]*Grsj*Ji[1][2] + Ji[0][1]*Grtj*Ji[2][2] + Ji[1][1]*Gsrj*Ji[0][2] + Ji[1][1]*Gssj*Ji[1][2] + Ji[1][1]*Gstj*Ji[2][2] + Ji[2][1]*Gtrj*Ji[0][2] + Ji[2][1]*Gtsj*Ji[1][2] + Ji[2][1]*Gttj*Ji[2][2];

				double GZX = Ji[0][2]*Grrj*Ji[0][0] + Ji[0][2]*Grsj*Ji[1][0] + Ji[0][2]*Grtj*Ji[2][0] + Ji[1][2]*Gsrj*Ji[0][0] + Ji[1][2]*Gssj*Ji[1][0] + Ji[1][2]*Gstj*Ji[2][0] + Ji[2][2]*Gtrj*Ji[0][0] + Ji[2][2]*Gtsj*Ji[1][0] + Ji[2][2]*Gttj*Ji[2][0];
				double GZY = Ji[0][2]*Grrj*Ji[0][1] + Ji[0][2]*Grsj*Ji[1][1] + Ji[0][2]*Grtj*Ji[2][1] + Ji[1][2]*Gsrj*Ji[0][1] + Ji[1][2]*Gssj*Ji[1][1] + Ji[1][2]*Gstj*Ji[2][1] + Ji[2][2]*Gtrj*Ji[0][1] + Ji[2][2]*Gtsj*Ji[1][1] + Ji[2][2]*Gttj*Ji[2][1];
				double GZZ = Ji[0][2]*Grrj*Ji[0][2] + Ji[0][2]*Grsj*Ji[1][2] + Ji[0][2]*Grtj*Ji[2][2] + Ji[1][2]*Gsrj*Ji[0][2] + Ji[1][2]*Gssj*Ji[1][2] + Ji[1][2]*Gstj*Ji[2][2] + Ji[2][2]*Gtrj*Ji[0][2] + Ji[2][2]*Gtsj*Ji[1][2] + Ji[2][2]*Gttj*Ji[2][2];
						
				kab = ((Gx[i]*(s.xx()*Gx[j] + s.xy()*Gy[j] + s.xz()*Gz[j])*detJt + (tau.d[0]*GXX + tau.d[1]*(GXY + GYX) + tau.d[2]*(GXZ + GZX) + tau.d[3]*GYY + tau.d[4]*(GYZ + GZY) + tau.d[5]*GZZ))*detJt*detJt +
					   (Gy[i]*(s.xy()*Gx[j] + s.yy()*Gy[j] + s.yz()*Gz[j])*detJt + (tau.d[1]*GXX + tau.d[3]*(GXY + GYX) + tau.d[4]*(GXZ + GZX) + tau.d[6]*GYY + tau.d[7]*(GYZ + GZY) + tau.d[8]*GZZ))*detJt*detJt + 
					   (Gz[i]*(s.xz()*Gx[j] + s.yz()*Gy[j] + s.zz()*Gz[j])*detJt + (tau.d[2]*GXX + tau.d[4]*(GXY + GYX) + tau.d[5]*(GXZ + GZX) + tau.d[7]*GYY + tau.d[8]*(GYZ + GZY) + tau.d[9]*GZZ))*detJt*detJt);

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix
void FEElasticMultiscaleDomain2O::ElementMaterialStiffness(FESolidElement &el, matrix &ke)
{
	int i, i3, j, j3, n;

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	// Gx = dH/dx
	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];

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

	// calculate element stiffness matrix
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n)*gw[n];
		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		//// LTE 
		double *Grrn = el.Grr(n);
		double *Grsn = el.Grs(n);
		double *Grtn = el.Grt(n);

		double *Gsrn = el.Gsr(n);
		double *Gssn = el.Gss(n);
		double *Gstn = el.Gst(n);

		double *Gtrn = el.Gtr(n);
		double *Gtsn = el.Gts(n);
		double *Gttn = el.Gtt(n);

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());

		// get the 'D' matrix
		tens4ds c; c.zero();
		tens5ds d; d.zero();
		tens6ds e; e.zero();

		FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);
		pmat->Tangent2O(mp, c, d, e);
		
		c.extract(D);

		double Grri, Grsi, Grti, Gsri, Gssi, Gsti, Gtri, Gtsi, Gtti;
		double Grrj, Grsj, Grtj, Gsrj, Gssj, Gstj, Gtrj, Gtsj, Gttj;

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

			Grri = Grrn[i];
			Grsi = Grsn[i];
			Grti = Grtn[i];

			Gsri = Gsrn[i];
			Gssi = Gssn[i];
			Gsti = Gstn[i];

			Gtri = Gtrn[i];
			Gtsi = Gtsn[i];
			Gtti = Gttn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !	
			// Calculate based on Gjk = Ji[l][j]*Glm*Ji[m][k]
			double GXXi = Ji[0][0]*Grri*Ji[0][0] + Ji[0][0]*Grsi*Ji[1][0] + Ji[0][0]*Grti*Ji[2][0] + Ji[1][0]*Gsri*Ji[0][0] + Ji[1][0]*Gssi*Ji[1][0] + Ji[1][0]*Gsti*Ji[2][0] + Ji[2][0]*Gtri*Ji[0][0] + Ji[2][0]*Gtsi*Ji[1][0] + Ji[2][0]*Gtti*Ji[2][0];
			double GXYi = Ji[0][0]*Grri*Ji[0][1] + Ji[0][0]*Grsi*Ji[1][1] + Ji[0][0]*Grti*Ji[2][1] + Ji[1][0]*Gsri*Ji[0][1] + Ji[1][0]*Gssi*Ji[1][1] + Ji[1][0]*Gsti*Ji[2][1] + Ji[2][0]*Gtri*Ji[0][1] + Ji[2][0]*Gtsi*Ji[1][1] + Ji[2][0]*Gtti*Ji[2][1];
			double GXZi = Ji[0][0]*Grri*Ji[0][2] + Ji[0][0]*Grsi*Ji[1][2] + Ji[0][0]*Grti*Ji[2][2] + Ji[1][0]*Gsri*Ji[0][2] + Ji[1][0]*Gssi*Ji[1][2] + Ji[1][0]*Gsti*Ji[2][2] + Ji[2][0]*Gtri*Ji[0][2] + Ji[2][0]*Gtsi*Ji[1][2] + Ji[2][0]*Gtti*Ji[2][2];
		
			double GYXi = Ji[0][1]*Grri*Ji[0][0] + Ji[0][1]*Grsi*Ji[1][0] + Ji[0][1]*Grti*Ji[2][0] + Ji[1][1]*Gsri*Ji[0][0] + Ji[1][1]*Gssi*Ji[1][0] + Ji[1][1]*Gsti*Ji[2][0] + Ji[2][1]*Gtri*Ji[0][0] + Ji[2][1]*Gtsi*Ji[1][0] + Ji[2][1]*Gtti*Ji[2][0];
			double GYYi = Ji[0][1]*Grri*Ji[0][1] + Ji[0][1]*Grsi*Ji[1][1] + Ji[0][1]*Grti*Ji[2][1] + Ji[1][1]*Gsri*Ji[0][1] + Ji[1][1]*Gssi*Ji[1][1] + Ji[1][1]*Gsti*Ji[2][1] + Ji[2][1]*Gtri*Ji[0][1] + Ji[2][1]*Gtsi*Ji[1][1] + Ji[2][1]*Gtti*Ji[2][1];
			double GYZi = Ji[0][1]*Grri*Ji[0][2] + Ji[0][1]*Grsi*Ji[1][2] + Ji[0][1]*Grti*Ji[2][2] + Ji[1][1]*Gsri*Ji[0][2] + Ji[1][1]*Gssi*Ji[1][2] + Ji[1][1]*Gsti*Ji[2][2] + Ji[2][1]*Gtri*Ji[0][2] + Ji[2][1]*Gtsi*Ji[1][2] + Ji[2][1]*Gtti*Ji[2][2];

			double GZXi = Ji[0][2]*Grri*Ji[0][0] + Ji[0][2]*Grsi*Ji[1][0] + Ji[0][2]*Grti*Ji[2][0] + Ji[1][2]*Gsri*Ji[0][0] + Ji[1][2]*Gssi*Ji[1][0] + Ji[1][2]*Gsti*Ji[2][0] + Ji[2][2]*Gtri*Ji[0][0] + Ji[2][2]*Gtsi*Ji[1][0] + Ji[2][2]*Gtti*Ji[2][0];
			double GZYi = Ji[0][2]*Grri*Ji[0][1] + Ji[0][2]*Grsi*Ji[1][1] + Ji[0][2]*Grti*Ji[2][1] + Ji[1][2]*Gsri*Ji[0][1] + Ji[1][2]*Gssi*Ji[1][1] + Ji[1][2]*Gsti*Ji[2][1] + Ji[2][2]*Gtri*Ji[0][1] + Ji[2][2]*Gtsi*Ji[1][1] + Ji[2][2]*Gtti*Ji[2][1];
			double GZZi = Ji[0][2]*Grri*Ji[0][2] + Ji[0][2]*Grsi*Ji[1][2] + Ji[0][2]*Grti*Ji[2][2] + Ji[1][2]*Gsri*Ji[0][2] + Ji[1][2]*Gssi*Ji[1][2] + Ji[1][2]*Gsti*Ji[2][2] + Ji[2][2]*Gtri*Ji[0][2] + Ji[2][2]*Gtsi*Ji[1][2] + Ji[2][2]*Gtti*Ji[2][2];
		
			for (j=i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = Gx[j];
				Gyj = Gy[j];
				Gzj = Gz[j];

				Grrj = Grrn[j];
				Grsj = Grsn[j];
				Grtj = Grtn[j];

				Gsrj = Gsrn[j];
				Gssj = Gssn[j];
				Gstj = Gstn[j];

				Gtrj = Gtrn[j];
				Gtsj = Gtsn[j];
				Gttj = Gttn[j];

				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !	
				// Calculate based on Gjk = Ji[l][j]*Glm*Ji[m][k]
				double GXXj = Ji[0][0]*Grrj*Ji[0][0] + Ji[0][0]*Grsj*Ji[1][0] + Ji[0][0]*Grtj*Ji[2][0] + Ji[1][0]*Gsrj*Ji[0][0] + Ji[1][0]*Gssj*Ji[1][0] + Ji[1][0]*Gstj*Ji[2][0] + Ji[2][0]*Gtrj*Ji[0][0] + Ji[2][0]*Gtsj*Ji[1][0] + Ji[2][0]*Gttj*Ji[2][0];
				double GXYj = Ji[0][0]*Grrj*Ji[0][1] + Ji[0][0]*Grsj*Ji[1][1] + Ji[0][0]*Grtj*Ji[2][1] + Ji[1][0]*Gsrj*Ji[0][1] + Ji[1][0]*Gssj*Ji[1][1] + Ji[1][0]*Gstj*Ji[2][1] + Ji[2][0]*Gtrj*Ji[0][1] + Ji[2][0]*Gtsj*Ji[1][1] + Ji[2][0]*Gttj*Ji[2][1];
				double GXZj = Ji[0][0]*Grrj*Ji[0][2] + Ji[0][0]*Grsj*Ji[1][2] + Ji[0][0]*Grtj*Ji[2][2] + Ji[1][0]*Gsrj*Ji[0][2] + Ji[1][0]*Gssj*Ji[1][2] + Ji[1][0]*Gstj*Ji[2][2] + Ji[2][0]*Gtrj*Ji[0][2] + Ji[2][0]*Gtsj*Ji[1][2] + Ji[2][0]*Gttj*Ji[2][2];
	
				double GYXj = Ji[0][1]*Grrj*Ji[0][0] + Ji[0][1]*Grsj*Ji[1][0] + Ji[0][1]*Grtj*Ji[2][0] + Ji[1][1]*Gsrj*Ji[0][0] + Ji[1][1]*Gssj*Ji[1][0] + Ji[1][1]*Gstj*Ji[2][0] + Ji[2][1]*Gtrj*Ji[0][0] + Ji[2][1]*Gtsj*Ji[1][0] + Ji[2][1]*Gttj*Ji[2][0];
				double GYYj = Ji[0][1]*Grrj*Ji[0][1] + Ji[0][1]*Grsj*Ji[1][1] + Ji[0][1]*Grtj*Ji[2][1] + Ji[1][1]*Gsrj*Ji[0][1] + Ji[1][1]*Gssj*Ji[1][1] + Ji[1][1]*Gstj*Ji[2][1] + Ji[2][1]*Gtrj*Ji[0][1] + Ji[2][1]*Gtsj*Ji[1][1] + Ji[2][1]*Gttj*Ji[2][1];
				double GYZj = Ji[0][1]*Grrj*Ji[0][2] + Ji[0][1]*Grsj*Ji[1][2] + Ji[0][1]*Grtj*Ji[2][2] + Ji[1][1]*Gsrj*Ji[0][2] + Ji[1][1]*Gssj*Ji[1][2] + Ji[1][1]*Gstj*Ji[2][2] + Ji[2][1]*Gtrj*Ji[0][2] + Ji[2][1]*Gtsj*Ji[1][2] + Ji[2][1]*Gttj*Ji[2][2];

				double GZXj = Ji[0][2]*Grrj*Ji[0][0] + Ji[0][2]*Grsj*Ji[1][0] + Ji[0][2]*Grtj*Ji[2][0] + Ji[1][2]*Gsrj*Ji[0][0] + Ji[1][2]*Gssj*Ji[1][0] + Ji[1][2]*Gstj*Ji[2][0] + Ji[2][2]*Gtrj*Ji[0][0] + Ji[2][2]*Gtsj*Ji[1][0] + Ji[2][2]*Gttj*Ji[2][0];
				double GZYj = Ji[0][2]*Grrj*Ji[0][1] + Ji[0][2]*Grsj*Ji[1][1] + Ji[0][2]*Grtj*Ji[2][1] + Ji[1][2]*Gsrj*Ji[0][1] + Ji[1][2]*Gssj*Ji[1][1] + Ji[1][2]*Gstj*Ji[2][1] + Ji[2][2]*Gtrj*Ji[0][1] + Ji[2][2]*Gtsj*Ji[1][1] + Ji[2][2]*Gttj*Ji[2][1];
				double GZZj = Ji[0][2]*Grrj*Ji[0][2] + Ji[0][2]*Grsj*Ji[1][2] + Ji[0][2]*Grtj*Ji[2][2] + Ji[1][2]*Gsrj*Ji[0][2] + Ji[1][2]*Gssj*Ji[1][2] + Ji[1][2]*Gstj*Ji[2][2] + Ji[2][2]*Gtrj*Ji[0][2] + Ji[2][2]*Gtsj*Ji[1][2] + Ji[2][2]*Gttj*Ji[2][2];
				
				// First order component (C*d)
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

				// Second-order components with D*d + D*k
				ke[i3  ][j3  ] += (Gxi*(d.d[0]*GXXj + d.d[1]*(GXYj + GYXj) + d.d[2]*(GXZj + GZXj) + d.d[3]*GYYj + d.d[4]*(GYZj + GZYj) + d.d[5]*GZZj)
					             + Gyi*(d.d[1]*GXXj + d.d[3]*(GXYj + GYXj) + d.d[4]*(GXZj + GZXj) + d.d[10]*GYYj + d.d[11]*(GYZj + GZYj) + d.d[12]*GZZj)
								 + Gzi*(d.d[2]*GXXj + d.d[4]*(GXYj + GYXj) + d.d[5]*(GXZj + GZXj) + d.d[11]*GYYj + d.d[12]*(GYZj + GZYj) + d.d[17]*GZZj))*detJt*detJt;
				
				ke[i3  ][j3+1] += (Gxi*(d.d[1]*GXXj + d.d[3]*(GXYj + GYXj) + d.d[4]*(GXZj + GZXj) + d.d[6]*GYYj + d.d[7]*(GYZj + GZYj) + d.d[8]*GZZj)
					             + Gyi*(d.d[3]*GXXj + d.d[10]*(GXYj + GYXj) + d.d[11]*(GXZj + GZXj) + d.d[13]*GYYj + d.d[14]*(GYZj + GZYj) + d.d[15]*GZZj)
								 + Gzi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[18]*GYYj + d.d[19]*(GYZj + GZYj) + d.d[20]*GZZj))*detJt*detJt;

				ke[i3  ][j3+2] += (Gxi*(d.d[2]*GXXj + d.d[4]*(GXYj + GYXj) + d.d[5]*(GXZj + GZXj) + d.d[7]*GYYj + d.d[8]*(GYZj + GZYj) + d.d[9]*GZZj)
					             + Gyi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[14]*GYYj + d.d[15]*(GYZj + GZYj) + d.d[16]*GZZj)
								 + Gzi*(d.d[5]*GXXj + d.d[12]*(GXYj + GYXj) + d.d[17]*(GXZj + GZXj) + d.d[19]*GYYj + d.d[20]*(GYZj + GZYj) + d.d[21]*GZZj))*detJt*detJt;

				ke[i3+1][j3  ] += (Gxi*(d.d[1]*GXXj + d.d[3]*(GXYj + GYXj) + d.d[4]*(GXZj + GZXj) + d.d[10]*GYYj + d.d[11]*(GYZj + GZYj) + d.d[12]*GZZj)
					             + Gyi*(d.d[3]*GXXj + d.d[10]*(GXYj + GYXj) + d.d[11]*(GXZj + GZXj) + d.d[13]*GYYj + d.d[14]*(GYZj + GZYj) + d.d[19]*GZZj)
								 + Gzi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[14]*GYYj + d.d[15]*(GYZj + GZYj) + d.d[20]*GZZj))*detJt*detJt;
				
				ke[i3+1][j3+1] += (Gxi*(d.d[3]*GXXj + d.d[10]*(GXYj + GYXj) + d.d[11]*(GXZj + GZXj) + d.d[13]*GYYj + d.d[14]*(GYZj + GZYj) + d.d[15]*GZZj)
					             + Gyi*(d.d[10]*GXXj + d.d[13]*(GXYj + GYXj) + d.d[14]*(GXZj + GZXj) + d.d[22]*GYYj + d.d[23]*(GYZj + GZYj) + d.d[24]*GZZj)
								 + Gzi*(d.d[11]*GXXj + d.d[14]*(GXYj + GYXj) + d.d[15]*(GXZj + GZXj) + d.d[23]*GYYj + d.d[24]*(GYZj + GZYj) + d.d[26]*GZZj))*detJt*detJt;

				ke[i3+1][j3+2] += (Gxi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[14]*GYYj + d.d[15]*(GYZj + GZYj) + d.d[16]*GZZj)
					             + Gyi*(d.d[11]*GXXj + d.d[14]*(GXYj + GYXj) + d.d[19]*(GXZj + GZXj) + d.d[23]*GYYj + d.d[24]*(GYZj + GZYj) + d.d[25]*GZZj)
								 + Gzi*(d.d[12]*GXXj + d.d[15]*(GXYj + GYXj) + d.d[20]*(GXZj + GZXj) + d.d[24]*GYYj + d.d[26]*(GYZj + GZYj) + d.d[27]*GZZj))*detJt*detJt;

				ke[i3+2][j3  ] += (Gxi*(d.d[2]*GXXj + d.d[4]*(GXYj + GYXj) + d.d[5]*(GXZj + GZXj) + d.d[11]*GYYj + d.d[12]*(GYZj + GZYj) + d.d[17]*GZZj)
					             + Gyi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[14]*GYYj + d.d[15]*(GYZj + GZYj) + d.d[20]*GZZj)
								 + Gzi*(d.d[5]*GXXj + d.d[12]*(GXYj + GYXj) + d.d[17]*(GXZj + GZXj) + d.d[19]*GYYj + d.d[20]*(GYZj + GZYj) + d.d[21]*GZZj))*detJt*detJt;

				ke[i3+2][j3+1] += (Gxi*(d.d[4]*GXXj + d.d[11]*(GXYj + GYXj) + d.d[12]*(GXZj + GZXj) + d.d[14]*GYYj + d.d[15]*(GYZj + GZYj) + d.d[20]*GZZj)
					             + Gyi*(d.d[11]*GXXj + d.d[14]*(GXYj + GYXj) + d.d[15]*(GXZj + GZXj) + d.d[23]*GYYj + d.d[24]*(GYZj + GZYj) + d.d[26]*GZZj)
								 + Gzi*(d.d[12]*GXXj + d.d[19]*(GXYj + GYXj) + d.d[20]*(GXZj + GZXj) + d.d[24]*GYYj + d.d[26]*(GYZj + GZYj) + d.d[27]*GZZj))*detJt*detJt;

				ke[i3+2][j3+2] += (Gxi*(d.d[5]*GXXj + d.d[12]*(GXYj + GYXj) + d.d[17]*(GXZj + GZXj) + d.d[15]*GYYj + d.d[20]*(GYZj + GZYj) + d.d[21]*GZZj)
					             + Gyi*(d.d[12]*GXXj + d.d[15]*(GXYj + GYXj) + d.d[20]*(GXZj + GZXj) + d.d[24]*GYYj + d.d[26]*(GYZj + GZYj) + d.d[27]*GZZj)
								 + Gzi*(d.d[17]*GXXj + d.d[20]*(GXYj + GYXj) + d.d[21]*(GXZj + GZXj) + d.d[26]*GYYj + d.d[27]*(GYZj + GZYj) + d.d[28]*GZZj))*detJt*detJt;

				
				ke[i3  ][j3  ] += (Gxj*(d.d[0]*GXXi + d.d[1]*(GXYi + GYXi) + d.d[2]*(GXZi + GZXi) + d.d[3]*GYYi + d.d[4]*(GYZi + GZYi) + d.d[5]*GZZi)
					             + Gyj*(d.d[1]*GXXi + d.d[3]*(GXYi + GYXi) + d.d[4]*(GXZi + GZXi) + d.d[10]*GYYi + d.d[11]*(GYZi + GZYi) + d.d[12]*GZZi)
								 + Gzj*(d.d[2]*GXXi + d.d[4]*(GXYi + GYXi) + d.d[5]*(GXZi + GZXi) + d.d[11]*GYYi + d.d[12]*(GYZi + GZYi) + d.d[17]*GZZi))*detJt*detJt;
		
				ke[i3  ][j3+1] += (Gxj*(d.d[1]*GXXi + d.d[3]*(GXYi + GYXi) + d.d[4]*(GXZi + GZXi) + d.d[6]*GYYi + d.d[7]*(GYZi + GZYi) + d.d[8]*GZZi)
					             + Gyj*(d.d[3]*GXXi + d.d[10]*(GXYi + GYXi) + d.d[11]*(GXZi + GZXi) + d.d[13]*GYYi + d.d[14]*(GYZi + GZYi) + d.d[15]*GZZi)
								 + Gzj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[18]*GYYi + d.d[19]*(GYZi + GZYi) + d.d[20]*GZZi))*detJt*detJt;

				ke[i3  ][j3+2] += (Gxj*(d.d[2]*GXXi + d.d[4]*(GXYi + GYXi) + d.d[5]*(GXZi + GZXi) + d.d[7]*GYYi + d.d[8]*(GYZi + GZYi) + d.d[9]*GZZi)
					             + Gyj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[14]*GYYi + d.d[15]*(GYZi + GZYi) + d.d[16]*GZZi)
								 + Gzj*(d.d[5]*GXXi + d.d[12]*(GXYi + GYXi) + d.d[17]*(GXZi + GZXi) + d.d[19]*GYYi + d.d[20]*(GYZi + GZYi) + d.d[21]*GZZi))*detJt*detJt;

				ke[i3+1][j3  ] += (Gxj*(d.d[1]*GXXi + d.d[3]*(GXYi + GYXi) + d.d[4]*(GXZi + GZXi) + d.d[10]*GYYi + d.d[11]*(GYZi + GZYi) + d.d[12]*GZZi)
					             + Gyj*(d.d[3]*GXXi + d.d[10]*(GXYi + GYXi) + d.d[11]*(GXZi + GZXi) + d.d[13]*GYYi + d.d[14]*(GYZi + GZYi) + d.d[19]*GZZi)
								 + Gzj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[14]*GYYi + d.d[15]*(GYZi + GZYi) + d.d[20]*GZZi))*detJt*detJt;
				
				ke[i3+1][j3+1] += (Gxj*(d.d[3]*GXXi + d.d[10]*(GXYi + GYXi) + d.d[11]*(GXZi + GZXi) + d.d[13]*GYYi + d.d[14]*(GYZi + GZYi) + d.d[15]*GZZi)
					             + Gyj*(d.d[10]*GXXi + d.d[13]*(GXYi + GYXi) + d.d[14]*(GXZi + GZXi) + d.d[22]*GYYi + d.d[23]*(GYZi + GZYi) + d.d[24]*GZZi)
								 + Gzj*(d.d[11]*GXXi + d.d[14]*(GXYi + GYXi) + d.d[15]*(GXZi + GZXi) + d.d[23]*GYYi + d.d[24]*(GYZi + GZYi) + d.d[26]*GZZi))*detJt*detJt;

				ke[i3+1][j3+2] += (Gxj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[14]*GYYi + d.d[15]*(GYZi + GZYi) + d.d[16]*GZZi)
					             + Gyj*(d.d[11]*GXXi + d.d[14]*(GXYi + GYXi) + d.d[19]*(GXZi + GZXi) + d.d[23]*GYYi + d.d[24]*(GYZi + GZYi) + d.d[25]*GZZi)
								 + Gzj*(d.d[12]*GXXi + d.d[15]*(GXYi + GYXi) + d.d[20]*(GXZi + GZXi) + d.d[24]*GYYi + d.d[26]*(GYZi + GZYi) + d.d[27]*GZZi))*detJt*detJt;

				ke[i3+2][j3  ] += (Gxj*(d.d[2]*GXXi + d.d[4]*(GXYi + GYXi) + d.d[5]*(GXZi + GZXi) + d.d[11]*GYYi + d.d[12]*(GYZi + GZYi) + d.d[17]*GZZi)
					             + Gyj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[14]*GYYi + d.d[15]*(GYZi + GZYi) + d.d[20]*GZZi)
								 + Gzj*(d.d[5]*GXXi + d.d[12]*(GXYi + GYXi) + d.d[17]*(GXZi + GZXi) + d.d[19]*GYYi + d.d[20]*(GYZi + GZYi) + d.d[21]*GZZi))*detJt*detJt;

				ke[i3+2][j3+1] += (Gxj*(d.d[4]*GXXi + d.d[11]*(GXYi + GYXi) + d.d[12]*(GXZi + GZXi) + d.d[14]*GYYi + d.d[15]*(GYZi + GZYi) + d.d[20]*GZZi)
					             + Gyj*(d.d[11]*GXXi + d.d[14]*(GXYi + GYXi) + d.d[15]*(GXZi + GZXi) + d.d[23]*GYYi + d.d[24]*(GYZi + GZYi) + d.d[26]*GZZi)
								 + Gzj*(d.d[12]*GXXi + d.d[19]*(GXYi + GYXi) + d.d[20]*(GXZi + GZXi) + d.d[24]*GYYi + d.d[26]*(GYZi + GZYi) + d.d[27]*GZZi))*detJt*detJt;

				ke[i3+2][j3+2] += (Gxj*(d.d[5]*GXXi + d.d[12]*(GXYi + GYXi) + d.d[17]*(GXZi + GZXi) + d.d[15]*GYYi + d.d[20]*(GYZi + GZYi) + d.d[21]*GZZi)
					             + Gyj*(d.d[12]*GXXi + d.d[15]*(GXYi + GYXi) + d.d[20]*(GXZi + GZXi) + d.d[24]*GYYi + d.d[26]*(GYZi + GZYi) + d.d[27]*GZZi)
								 + Gzj*(d.d[17]*GXXi + d.d[20]*(GXYi + GYXi) + d.d[21]*(GXZi + GZXi) + d.d[26]*GYYi + d.d[27]*(GYZi + GZYi) + d.d[28]*GZZi))*detJt*detJt;

				// Second-order component E*k
				ke[i3  ][j3  ] += (GXXi*(e.d[0]*GXXj + e.d[1]*(GXYj + GYXj) + e.d[2]*(GXZj + GZXj) + e.d[3]*GYYj + e.d[4]*(GYZj + GZYj) + e.d[5]*GZZj)
					            + (GXYi + GYXi)*(e.d[1]*GXXj + e.d[3]*(GXYj + GYXj) + e.d[4]*(GXZj + GZXj) + e.d[10]*GYYj + e.d[11]*(GYZj + GZYj) + e.d[12]*GZZj)
								+ (GXZi + GZXi)*(e.d[2]*GXXj + e.d[4]*(GXYj + GYXj) + e.d[5]*(GXZj + GZXj) + e.d[17]*GYYj + e.d[18]*(GYZj + GZYj) + e.d[19]*GZZj)
								 + GYYi*(e.d[3]*GXXj + e.d[10]*(GXYj + GYXj) + e.d[11]*(GXZj + GZXj) + e.d[13]*GYYj + e.d[14]*(GYZj + GZYj) + e.d[24]*GZZj)
								+ (GYZi + GZYi)*(e.d[4]*GXXj + e.d[11]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[24]*(GYZj + GZYj) + e.d[29]*GZZj)
								 + GZZi*(e.d[5]*GXXj + e.d[18]*(GXYj + GYXj) + e.d[19]*(GXZj + GZXj) + e.d[24]*GYYj + e.d[22]*(GYZj + GZYj) + e.d[23]*GZZj))*detJt*detJt;

				ke[i3  ][j3+1] += (GXXi*(e.d[1]*GXXj + e.d[3]*(GXYj + GYXj) + e.d[4]*(GXZj + GZXj) + e.d[6]*GYYj + e.d[7]*(GYZj + GZYj) + e.d[8]*GZZj)
					            + (GXYi + GYXi)*(e.d[3]*GXXj + e.d[10]*(GXYj + GYXj) + e.d[11]*(GXZj + GZXj) + e.d[13]*GYYj + e.d[14]*(GYZj + GZYj) + e.d[15]*GZZj)
								+ (GXZi + GZXi)*(e.d[4]*GXXj + e.d[17]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[20]*GYYj + e.d[21]*(GYZj + GZYj) + e.d[22]*GZZj)
								 + GYYi*(e.d[10]*GXXj + e.d[13]*(GXYj + GYXj) + e.d[14]*(GXZj + GZXj) + e.d[25]*GYYj + e.d[26]*(GYZj + GZYj) + e.d[27]*GZZj)
								+ (GYZi + GZYi)*(e.d[11]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[24]*(GXZj + GZXj) + e.d[30]*GYYj + e.d[31]*(GYZj + GZYj) + e.d[32]*GZZj)
								 + GZZi*(e.d[18]*GXXj + e.d[24]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[34]*GYYj + e.d[35]*(GYZj + GZYj) + e.d[36]*GZZj))*detJt*detJt;

				ke[i3  ][j3+2] += (GXXi*(e.d[2]*GXXj + e.d[4]*(GXYj + GYXj) + e.d[5]*(GXZj + GZXj) + e.d[7]*GYYj + e.d[8]*(GYZj + GZYj) + e.d[9]*GZZj)
					            + (GXYi + GYXi)*(e.d[4]*GXXj + e.d[11]*(GXYj + GYXj) + e.d[12]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[15]*(GYZj + GZYj) + e.d[16]*GZZj)
								+ (GXZi + GZXi)*(e.d[5]*GXXj + e.d[18]*(GXYj + GYXj) + e.d[19]*(GXZj + GZXj) + e.d[21]*GYYj + e.d[22]*(GYZj + GZYj) + e.d[23]*GZZj)
								 + GYYi*(e.d[11]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[24]*(GXZj + GZXj) + e.d[26]*GYYj + e.d[27]*(GYZj + GZYj) + e.d[28]*GZZj)
								+ (GYZi + GZYi)*(e.d[18]*GXXj + e.d[24]*(GXYj + GYXj) + e.d[29]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[32]*(GYZj + GZYj) + e.d[33]*GZZj)
								 + GZZi*(e.d[19]*GXXj + e.d[22]*(GXYj + GYXj) + e.d[23]*(GXZj + GZXj) + e.d[35]*GYYj + e.d[36]*(GYZj + GZYj) + e.d[37]*GZZj))*detJt*detJt;

				ke[i3+1][j3  ] += (GXXi*(e.d[1]*GXXj + e.d[3]*(GXYj + GYXj) + e.d[4]*(GXZj + GZXj) + e.d[10]*GYYj + e.d[11]*(GYZj + GZYj) + e.d[18]*GZZj)
					            + (GXYi + GYXi)*(e.d[3]*GXXj + e.d[10]*(GXYj + GYXj) + e.d[11]*(GXZj + GZXj) + e.d[13]*GYYj + e.d[14]*(GYZj + GZYj) + e.d[21]*GZZj)
								+ (GXZi + GZXi)*(e.d[4]*GXXj + e.d[11]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[21]*(GYZj + GZYj) + e.d[22]*GZZj)
								 + GYYi*(e.d[10]*GXXj + e.d[13]*(GXYj + GYXj) + e.d[14]*(GXZj + GZXj) + e.d[25]*GYYj + e.d[26]*(GYZj + GZYj) + e.d[31]*GZZj)
								+ (GYZi + GZYi)*(e.d[17]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[24]*(GXZj + GZXj) + e.d[26]*GYYj + e.d[31]*(GYZj + GZYj) + e.d[35]*GZZj)
								 + GZZi*(e.d[18]*GXXj + e.d[21]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[35]*(GYZj + GZYj) + e.d[36]*GZZj))*detJt*detJt;

				ke[i3+1][j3+1] += (GXXi*(e.d[3]*GXXj + e.d[10]*(GXYj + GYXj) + e.d[11]*(GXZj + GZXj) + e.d[13]*GYYj + e.d[14]*(GYZj + GZYj) + e.d[15]*GZZj)
					            + (GXYi + GYXi)*(e.d[10]*GXXj + e.d[13]*(GXYj + GYXj) + e.d[14]*(GXZj + GZXj) + e.d[25]*GYYj + e.d[26]*(GYZj + GZYj) + e.d[27]*GZZj)
								+ (GXZi + GZXi)*(e.d[11]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[21]*(GXZj + GZXj) + e.d[30]*GYYj + e.d[31]*(GYZj + GZYj) + e.d[32]*GZZj)
								 + GYYi*(e.d[13]*GXXj + e.d[25]*(GXYj + GYXj) + e.d[26]*(GXZj + GZXj) + e.d[38]*GYYj + e.d[39]*(GYZj + GZYj) + e.d[40]*GZZj)
								+ (GYZi + GZYi)*(e.d[14]*GXXj + e.d[26]*(GXYj + GYXj) + e.d[31]*(GXZj + GZXj) + e.d[39]*GYYj + e.d[40]*(GYZj + GZYj) + e.d[42]*GZZj)
								 + GZZi*(e.d[21]*GXXj + e.d[31]*(GXYj + GYXj) + e.d[35]*(GXZj + GZXj) + e.d[40]*GYYj + e.d[42]*(GYZj + GZYj) + e.d[43]*GZZj))*detJt*detJt;

				ke[i3+1][j3+2] += (GXXi*(e.d[4]*GXXj + e.d[11]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[15]*(GYZj + GZYj) + e.d[16]*GZZj)
					            + (GXYi + GYXi)*(e.d[11]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[21]*(GXZj + GZXj) + e.d[26]*GYYj + e.d[27]*(GYZj + GZYj) + e.d[28]*GZZj)
								+ (GXZi + GZXi)*(e.d[18]*GXXj + e.d[21]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[32]*(GYZj + GZYj) + e.d[33]*GZZj)
								 + GYYi*(e.d[14]*GXXj + e.d[26]*(GXYj + GYXj) + e.d[31]*(GXZj + GZXj) + e.d[39]*GYYj + e.d[40]*(GYZj + GZYj) + e.d[41]*GZZj)
								+ (GYZi + GZYi)*(e.d[24]*GXXj + e.d[31]*(GXYj + GYXj) + e.d[35]*(GXZj + GZXj) + e.d[40]*GYYj + e.d[42]*(GYZj + GZYj) + e.d[43]*GZZj)
								 + GZZi*(e.d[22]*GXXj + e.d[35]*(GXYj + GYXj) + e.d[36]*(GXZj + GZXj) + e.d[42]*GYYj + e.d[43]*(GYZj + GZYj) + e.d[44]*GZZj))*detJt*detJt;

				ke[i3+2][j3  ] += (GXXi*(e.d[2]*GXXj + e.d[4]*(GXYj + GYXj) + e.d[5]*(GXZj + GZXj) + e.d[17]*GYYj + e.d[18]*(GYZj + GZYj) + e.d[19]*GZZj)
					            + (GXYi + GYXi)*(e.d[4]*GXXj + e.d[11]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[21]*(GYZj + GZYj) + e.d[22]*GZZj)
								+ (GXZi + GZXi)*(e.d[5]*GXXj + e.d[18]*(GXYj + GYXj) + e.d[19]*(GXZj + GZXj) + e.d[24]*GYYj + e.d[22]*(GYZj + GZYj) + e.d[23]*GZZj)
								 + GYYi*(e.d[17]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[24]*(GXZj + GZXj) + e.d[26]*GYYj + e.d[31]*(GYZj + GZYj) + e.d[35]*GZZj)
								+ (GYZi + GZYi)*(e.d[18]*GXXj + e.d[21]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[35]*(GYZj + GZYj) + e.d[36]*GZZj)
								 + GZZi*(e.d[19]*GXXj + e.d[22]*(GXYj + GYXj) + e.d[23]*(GXZj + GZXj) + e.d[35]*GYYj + e.d[36]*(GYZj + GZYj) + e.d[37]*GZZj))*detJt*detJt;

				ke[i3+2][j3+1] += (GXXi*(e.d[4]*GXXj + e.d[17]*(GXYj + GYXj) + e.d[18]*(GXZj + GZXj) + e.d[14]*GYYj + e.d[21]*(GYZj + GZYj) + e.d[22]*GZZj)
					            + (GXYi + GYXi)*(e.d[11]*GXXj + e.d[14]*(GXYj + GYXj) + e.d[21]*(GXZj + GZXj) + e.d[26]*GYYj + e.d[31]*(GYZj + GZYj) + e.d[35]*GZZj)
								+ (GXZi + GZXi)*(e.d[18]*GXXj + e.d[24]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[35]*(GYZj + GZYj) + e.d[36]*GZZj)
								 + GYYi*(e.d[14]*GXXj + e.d[26]*(GXYj + GYXj) + e.d[31]*(GXZj + GZXj) + e.d[39]*GYYj + e.d[40]*(GYZj + GZYj) + e.d[42]*GZZj)
								+ (GYZi + GZYi)*(e.d[21]*GXXj + e.d[31]*(GXYj + GYXj) + e.d[35]*(GXZj + GZXj) + e.d[40]*GYYj + e.d[42]*(GYZj + GZYj) + e.d[43]*GZZj)
								 + GZZi*(e.d[22]*GXXj + e.d[35]*(GXYj + GYXj) + e.d[36]*(GXZj + GZXj) + e.d[42]*GYYj + e.d[43]*(GYZj + GZYj) + e.d[44]*GZZj))*detJt*detJt;

				ke[i3+2][j3+2] += (GXXi*(e.d[5]*GXXj + e.d[18]*(GXYj + GYXj) + e.d[19]*(GXZj + GZXj) + e.d[21]*GYYj + e.d[22]*(GYZj + GZYj) + e.d[23]*GZZj)
					            + (GXYi + GYXi)*(e.d[18]*GXXj + e.d[21]*(GXYj + GYXj) + e.d[22]*(GXZj + GZXj) + e.d[31]*GYYj + e.d[35]*(GYZj + GZYj) + e.d[33]*GZZj)
								+ (GXZi + GZXi)*(e.d[19]*GXXj + e.d[22]*(GXYj + GYXj) + e.d[23]*(GXZj + GZXj) + e.d[35]*GYYj + e.d[36]*(GYZj + GZYj) + e.d[37]*GZZj)
								 + GYYi*(e.d[24]*GXXj + e.d[31]*(GXYj + GYXj) + e.d[35]*(GXZj + GZXj) + e.d[40]*GYYj + e.d[42]*(GYZj + GZYj) + e.d[43]*GZZj)
								+ (GYZi + GZYi)*(e.d[22]*GXXj + e.d[35]*(GXYj + GYXj) + e.d[36]*(GXZj + GZXj) + e.d[42]*GYYj + e.d[43]*(GYZj + GZYj) + e.d[44]*GZZj)
								 + GZZi*(e.d[23]*GXXj + e.d[36]*(GXYj + GYXj) + e.d[37]*(GXZj + GZXj) + e.d[43]*GYYj + e.d[44]*(GYZj + GZYj) + e.d[45]*GZZj))*detJt*detJt;

			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element el at integration point n.
//! The deformation gradient is returned in F and its determinant is the return
//! value of the function
void FEElasticMultiscaleDomain2O::defhess(FESolidElement &el, tens3drs &G, int n)
{
	int i;

	// number of nodes
	int neln = el.Nodes();

	// shape function derivatives
	double *Grn = el.Gr(n);
	double *Gsn = el.Gs(n);
	double *Gtn = el.Gt(n);

	// nodal points
	vec3d r[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// calculate inverse jacobian
	double Ji[3][3];
	invjac0(el, Ji, n);

	//// LTE 
	// Calculate G (Deformation Hessian)
	double *Grrn = el.Grr(n);
	double *Grsn = el.Grs(n);
	double *Grtn = el.Grt(n);

	double *Gsrn = el.Gsr(n);
	double *Gssn = el.Gss(n);
	double *Gstn = el.Gst(n);

	double *Gtrn = el.Gtr(n);
	double *Gtsn = el.Gts(n);
	double *Gttn = el.Gtt(n);
	
	G.zero();

	for (i=0; i<neln; ++i)
	{
		double Grri = Grrn[i];
		double Grsi = Grsn[i];
		double Grti = Grtn[i];

		double Gsri = Gsrn[i];
		double Gssi = Gssn[i];
		double Gsti = Gstn[i];

		double Gtri = Gtrn[i];
		double Gtsi = Gtsn[i];
		double Gtti = Gttn[i];

		double x = r[i].x;
		double y = r[i].y;
		double z = r[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !	
		// Calculate based on Gjk = Ji[l][j]*Glm*Ji[m][k]
		double GXX = Ji[0][0]*Grri*Ji[0][0] + Ji[0][0]*Grsi*Ji[1][0] + Ji[0][0]*Grti*Ji[2][0] + Ji[1][0]*Gsri*Ji[0][0] + Ji[1][0]*Gssi*Ji[1][0] + Ji[1][0]*Gsti*Ji[2][0] + Ji[2][0]*Gtri*Ji[0][0] + Ji[2][0]*Gtsi*Ji[1][0] + Ji[2][0]*Gtti*Ji[2][0];
		double GXY = Ji[0][0]*Grri*Ji[0][1] + Ji[0][0]*Grsi*Ji[1][1] + Ji[0][0]*Grti*Ji[2][1] + Ji[1][0]*Gsri*Ji[0][1] + Ji[1][0]*Gssi*Ji[1][1] + Ji[1][0]*Gsti*Ji[2][1] + Ji[2][0]*Gtri*Ji[0][1] + Ji[2][0]*Gtsi*Ji[1][1] + Ji[2][0]*Gtti*Ji[2][1];
		double GXZ = Ji[0][0]*Grri*Ji[0][2] + Ji[0][0]*Grsi*Ji[1][2] + Ji[0][0]*Grti*Ji[2][2] + Ji[1][0]*Gsri*Ji[0][2] + Ji[1][0]*Gssi*Ji[1][2] + Ji[1][0]*Gsti*Ji[2][2] + Ji[2][0]*Gtri*Ji[0][2] + Ji[2][0]*Gtsi*Ji[1][2] + Ji[2][0]*Gtti*Ji[2][2];
		
		double GYX = Ji[0][1]*Grri*Ji[0][0] + Ji[0][1]*Grsi*Ji[1][0] + Ji[0][1]*Grti*Ji[2][0] + Ji[1][1]*Gsri*Ji[0][0] + Ji[1][1]*Gssi*Ji[1][0] + Ji[1][1]*Gsti*Ji[2][0] + Ji[2][1]*Gtri*Ji[0][0] + Ji[2][1]*Gtsi*Ji[1][0] + Ji[2][1]*Gtti*Ji[2][0];
		double GYY = Ji[0][1]*Grri*Ji[0][1] + Ji[0][1]*Grsi*Ji[1][1] + Ji[0][1]*Grti*Ji[2][1] + Ji[1][1]*Gsri*Ji[0][1] + Ji[1][1]*Gssi*Ji[1][1] + Ji[1][1]*Gsti*Ji[2][1] + Ji[2][1]*Gtri*Ji[0][1] + Ji[2][1]*Gtsi*Ji[1][1] + Ji[2][1]*Gtti*Ji[2][1];
		double GYZ = Ji[0][1]*Grri*Ji[0][2] + Ji[0][1]*Grsi*Ji[1][2] + Ji[0][1]*Grti*Ji[2][2] + Ji[1][1]*Gsri*Ji[0][2] + Ji[1][1]*Gssi*Ji[1][2] + Ji[1][1]*Gsti*Ji[2][2] + Ji[2][1]*Gtri*Ji[0][2] + Ji[2][1]*Gtsi*Ji[1][2] + Ji[2][1]*Gtti*Ji[2][2];

		double GZX = Ji[0][2]*Grri*Ji[0][0] + Ji[0][2]*Grsi*Ji[1][0] + Ji[0][2]*Grti*Ji[2][0] + Ji[1][2]*Gsri*Ji[0][0] + Ji[1][2]*Gssi*Ji[1][0] + Ji[1][2]*Gsti*Ji[2][0] + Ji[2][2]*Gtri*Ji[0][0] + Ji[2][2]*Gtsi*Ji[1][0] + Ji[2][2]*Gtti*Ji[2][0];
		double GZY = Ji[0][2]*Grri*Ji[0][1] + Ji[0][2]*Grsi*Ji[1][1] + Ji[0][2]*Grti*Ji[2][1] + Ji[1][2]*Gsri*Ji[0][1] + Ji[1][2]*Gssi*Ji[1][1] + Ji[1][2]*Gsti*Ji[2][1] + Ji[2][2]*Gtri*Ji[0][1] + Ji[2][2]*Gtsi*Ji[1][1] + Ji[2][2]*Gtti*Ji[2][1];
		double GZZ = Ji[0][2]*Grri*Ji[0][2] + Ji[0][2]*Grsi*Ji[1][2] + Ji[0][2]*Grti*Ji[2][2] + Ji[1][2]*Gsri*Ji[0][2] + Ji[1][2]*Gssi*Ji[1][2] + Ji[1][2]*Gsti*Ji[2][2] + Ji[2][2]*Gtri*Ji[0][2] + Ji[2][2]*Gtsi*Ji[1][2] + Ji[2][2]*Gtti*Ji[2][2];

		// calculate G
		G.d[0] += GXX*x;   G.d[1] += 0.5*(GXY*x + GYX*x);  G.d[2] += 0.5*(GXZ*x + GZX*x);  G.d[3] += GYY*x;  G.d[4] += 0.5*(GYZ*x + GZY*x);  G.d[5] += GZZ*x; 
		G.d[6] += GXX*y;   G.d[7] += 0.5*(GXY*y + GYX*y);  G.d[8] += 0.5*(GXZ*y + GZX*y);  G.d[9] += GYY*y; G.d[10] += 0.5*(GYZ*y + GZY*y); G.d[11] += GZZ*y;
		G.d[12] += GXX*z; G.d[13] += 0.5*(GXY*z + GYX*z); G.d[14] += 0.5*(GXZ*z + GZX*z); G.d[15] += GYY*z; G.d[16] += 0.5*(GYZ*z + GZY*z); G.d[17] += GZZ*z;
	}
}
