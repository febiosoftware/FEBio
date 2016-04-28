#include "stdafx.h"
#include "FEElasticMultiscaleDomain2O.h"
#include "FEMicroMaterial2O.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"

//-----------------------------------------------------------------------------
// helper function for comparing two facets
bool compare_facets(int* na, int* nb, int nodes)
{
	switch (nodes)
	{
	case 3:
	case 6:
	case 7:
		if ((na[0]!=nb[0])&&(na[0]!=nb[1])&&(na[0]!=nb[2])) return false;
		if ((na[1]!=nb[0])&&(na[1]!=nb[1])&&(na[1]!=nb[2])) return false;
		if ((na[2]!=nb[0])&&(na[2]!=nb[1])&&(na[2]!=nb[2])) return false;
		break;
	default:
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// helper function for mapping surface facet coordinates to element facet coordinates
vec2d map_facet_to_facet(int* na, int* nb, int nodes, double r, double s)
{
	switch (nodes)
	{
	case 3:
	case 6:
	case 7:
		if ((na[0] == nb[0])&&(na[1] == nb[1])&&(na[2] == nb[2])) return vec2d(  r  ,     s);
		if ((na[0] == nb[2])&&(na[1] == nb[0])&&(na[2] == nb[1])) return vec2d(    s, 1-r-s);
		if ((na[0] == nb[1])&&(na[1] == nb[2])&&(na[2] == nb[0])) return vec2d(1-r-s,   r  );
		if ((na[0] == nb[0])&&(na[1] == nb[2])&&(na[2] == nb[1])) return vec2d(    s,   r  );
		if ((na[0] == nb[2])&&(na[1] == nb[1])&&(na[2] == nb[0])) return vec2d(  r  , 1-r-s);
		if ((na[0] == nb[1])&&(na[1] == nb[0])&&(na[2] == nb[2])) return vec2d(1-r-s,     s);
		break;
	}

	// we shouldn't get here
	assert(false);
	return vec2d(0,0);
}

//-----------------------------------------------------------------------------
// helper function that maps natural coordinates from facets to elements
vec3d map_facet_to_solid(FEMesh& mesh, FESurfaceElement& face, FESolidElement& el, double r, double s)
{
	int fn[FEElement::MAX_NODES];
	int nfaces = mesh.Faces(el);
	for (int i=0; i<nfaces; ++i)
	{
		mesh.GetFace(el, i, fn);
		if (compare_facets(&face.m_node[0], fn, face.Nodes()))
		{
			// map the facet coordinates to the element's facet coordinates
			// (faces can be rotated or inverted w.r.t. the element's face)
			vec2d b = map_facet_to_facet(&face.m_node[0], fn, face.Nodes(), r, s);
			double h1 = b.x(), h2 = b.y();

			// convert facet coordinates to volume coordinates
			// TODO: this is hard coded for tets! Need to generalize for hexes
			double g1, g2, g3;
			switch (i)
			{
			case 0: g1 = h1; g2 = h2; g3 = 0.; break;
			case 1: g1 = h1; g2 = 0.; g3 = h2; break;
			case 2: g1 = 0.; g2 = h1; g3 = h2; break;
			case 3: g1 = 1.0-h1-h2; g2 = h1; g3 = h2; break;
			}

			return vec3d(g1, g2, g3);
		}
	}

	// we shouldn't get here
	assert(false);
	return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
FEElasticMultiscaleDomain2O::FEInternalSurface2O::FEInternalSurface2O()
{
}

bool FEElasticMultiscaleDomain2O::FEInternalSurface2O::Initialize(FEElasticMultiscaleDomain2O* dom)
{
	// get the material and make sure it is correct
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(dom->GetMaterial());
	if (pmat == 0) return false;

	// build the inside surface
	FEMesh& mesh = *dom->GetMesh();
	m_ps = mesh.ElementBoundarySurface(false, true);
	if (m_ps == 0) return false;

	// allocate data
	int NF = m_ps->Elements();
	int nnf = 0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = m_ps->Element(i);
		nnf += el.GaussPoints();
	}
	m_data.resize(nnf);

	// allocate material point data
	nnf = 0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& face = m_ps->Element(i);
		int nint = face.GaussPoints();
		for (int n=0; n<nint; ++n, ++nnf)
		{
			m_data[nnf].m_pt[0] = pmat->CreateMaterialPointData();
			m_data[nnf].m_pt[1] = pmat->CreateMaterialPointData();

			m_data[nnf].m_pt[0]->Init(true);
			m_data[nnf].m_pt[1]->Init(true);

			// iso-parametric coordinates in surface element
			double h1 = face.gr(n);
			double h2 = face.gs(n);

			for (int k=0; k<2; ++k)
			{
				// get the adjacent solid element
				FESolidElement& ek = dom->Element(face.m_elem[k]);

				// map the iso-parametric coordinates from the facet to the solid element
				m_data[nnf].ksi[k] = map_facet_to_solid(mesh, face, ek, h1, h2);
			}
		}
	}

	return true;
}

//=============================================================================

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
			defhess(el, j, mmpt2O.m_G);

			mmpt2O.m_F_prev = pt.m_F;	//-> I don't think this does anything since Initialize is only called once
			mmpt2O.m_G_prev = mmpt2O.m_G;
			mmpt2O.m_rve.CopyFrom(rve);
			mmpt2O.m_rve.Init();
		}
	}

	// initialize the internal surface data
	if (m_surf.Initialize(this) == false) return false;

	// initialize surface RVEs
	int nnf = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& face = m_surf.Element(i);
		int nint = face.GaussPoints();
		for (int n=0; n<nint; ++n, ++nnf)
		{
			for (int k=0; k<2; ++k)
			{
				FEMaterialPoint& mp = *m_surf.GetData(i).m_pt[0];
				FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();

				mmpt2O.m_rve.CopyFrom(rve);
				mmpt2O.m_rve.Init();
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
void FEElasticMultiscaleDomain2O::InitElements()
{
	FEElasticSolidDomain::InitElements();

	int NF = m_surf.Elements(), nd = 0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = m_surf.Element(i);
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n, ++nd)
		{
			FEInternalSurface2O::Data& data = m_surf.GetData(nd);
			data.m_pt[0]->Init(false);
			data.m_pt[1]->Init(false);
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::InternalForces(FEGlobalVector& R)
{
	// call base class first
	FEElasticSolidDomain::InternalForces(R);

	// add the discrete-Galerkin contribution
	InternalForcesDG1(R);
	InternalForcesDG2(R);
}

//-----------------------------------------------------------------------------
//! Evaluate contribution of discrete-Galerkin enforcement of stress flux.
void FEElasticMultiscaleDomain2O::InternalForcesDG1(FEGlobalVector& R)
{
	FEMesh& mesh = *GetMesh();
	FESurface& surf = *m_surf.GetSurface();

	mat3d Ji;
	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];

	vector<double> fe;
	vector<int> lm;

	// loop over all internal surfaces elements
	int nd = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		// get the next surface
		FESurfaceElement& face = m_surf.Element(i);
		int nfn = face.Nodes();
		int nint = face.GaussPoints();

		// loop over both sides
		for (int m=0; m<2; ++m)
		{
			// evaluate the spatial gradient of shape functions
			FESolidElement& el = Element(face.m_elem[m]);
			int neln = el.Nodes();

			// sign
			double sgn = (m==0 ? 1.0 : -1.0);

			// allocate force vector
			int ndof = neln*3;
			fe.resize(ndof, 0.0);

			// loop over all the integration points
			double* gw = face.GaussWeights();
			for (int n=0; n<nint; ++n)
			{
				// get facet data
				FEInternalSurface2O::Data& data = m_surf.GetData(nd + n);

				// average stress across interface
				tens3drs& Qavg = data.Qavg;

				// calculate jacobian and surface normal
				vec3d nu(0,0,0);
				double J = surf.jac0(face, n, nu);

				// evaluate element Jacobian and shape function derivatives
				// at this integration point
				vec3d& ksi = data.ksi[m];
				invjac0(el, ksi.x, ksi.y, ksi.z, Ji);
				el.shape_deriv(Gr, Gs, Gt, ksi.x, ksi.y, ksi.z);

				// loop over element nodes
				for (int j=0; j<neln; ++j)
				{
					// calculate global gradient of shape functions
					// note that we need the transposed of Ji, not Ji itself !
					double G[3];
					G[0] = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
					G[1] = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
					G[2] = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];

					// put it all together
					double f[3] = {0}, Nu[3] = {nu.x, nu.y, nu.z};
					for (int k=0; k<3; ++k)
						for (int l=0; l<3; ++l)
						{
							f[0] += Qavg(0,k,l)*G[k]*Nu[l];
							f[1] += Qavg(1,k,l)*G[k]*Nu[l];
							f[2] += Qavg(2,k,l)*G[k]*Nu[l];
						}

					// the negative sign is because we need to subtract the internal forces
					// from the residual
					fe[3*j  ] -= f[0]*gw[n]*J*sgn;
					fe[3*j+1] -= f[1]*gw[n]*J*sgn;
					fe[3*j+2] -= f[2]*gw[n]*J*sgn;
				}
			}

			// unpack the LM values
			UnpackLM(el, lm);

			// assemble 
			R.Assemble(el.m_node, lm, fe);
		}

		// don't forgot to increment data counter
		nd += nint;
	}
}

//-----------------------------------------------------------------------------
//! Evaluate contribution of discrete-Galerkin enforcement of stress flux.
void FEElasticMultiscaleDomain2O::InternalForcesDG2(FEGlobalVector& R)
{
	FEMesh& mesh = *GetMesh();
	FESurface& surf = *m_surf.GetSurface();

	mat3d Ji;
	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];

	vector<double> fe;
	vector<int> lm;

	// loop over all internal surfaces elements
	int nd = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		// get the next surface
		FESurfaceElement& face = m_surf.Element(i);
		int nfn = face.Nodes();
		int nint = face.GaussPoints();

		// loop over both sides
		for (int m=0; m<2; ++m)
		{
			// evaluate the spatial gradient of shape functions
			FESolidElement& el = Element(face.m_elem[m]);
			int neln = el.Nodes();

			// sign
			double sgn = (m==0 ? 1.0 : -1.0);

			// allocate force vector
			int ndof = neln*3;
			fe.resize(ndof, 0.0);

			// loop over all the integration points
			double* gw = face.GaussWeights();
			for (int n=0; n<nint; ++n)
			{
				// get facet data
				FEInternalSurface2O::Data& data = m_surf.GetData(nd + n);

				// average stiffness across interface
				tens6ds& J0avg = data.J0avg;

				// displacement gradient jump
				const mat3d& Du = data.DgradU;

				// calculate jacobian and surface normal
				vec3d nu(0,0,0);
				double J = surf.jac0(face, n, nu);

				// evaluate element Jacobian and shape function derivatives
				// at this integration point
				vec3d& ksi = data.ksi[m];
				invjac0(el, ksi.x, ksi.y, ksi.z, Ji);
				el.shape_deriv(Gr, Gs, Gt, ksi.x, ksi.y, ksi.z);

				// loop over element nodes
				for (int j=0; j<neln; ++j)
				{
					// calculate global gradient of shape functions
					// note that we need the transposed of Ji, not Ji itself !
					double G[3];
					G[0] = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
					G[1] = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
					G[2] = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];

					// put it all together
					double f[3] = {0}, Nu[3] = {nu.x, nu.y, nu.z};
					for (int k=0; k<3; ++k)
						for (int l=0; l<3; ++l)
							for (int p=0; p<3; ++p)
								for (int q=0; q<3; ++q)
									for (int r=0; r<3; ++r)
									{
										f[0] += Du(k,l)*Nu[p]*J0avg(k,l,p,0,q,r)*G[q]*Nu[r];
										f[1] += Du(k,l)*Nu[p]*J0avg(k,l,p,1,q,r)*G[q]*Nu[r];
										f[2] += Du(k,l)*Nu[p]*J0avg(k,l,p,2,q,r)*G[q]*Nu[r];
									}

					// the negative sign is because we need to subtract the internal forces
					// from the residual
					fe[3*j  ] -= f[0]*gw[n]*J*sgn;
					fe[3*j+1] -= f[1]*gw[n]*J*sgn;
					fe[3*j+2] -= f[2]*gw[n]*J*sgn;
				}
			}

			// unpack the LM values
			UnpackLM(el, lm);

			// assemble 
			R.Assemble(el.m_node, lm, fe);
		}

		// don't forgot to increment data counter
		nd += nint;
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements
void FEElasticMultiscaleDomain2O::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	// contribution from stress
	ElementInternalForce_PF(el, fe);

	// contriubtion from higher-order stress
	ElementInternalForce_QG(el, fe);
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::ElementInternalForce_PF(FESolidElement& el, vector<double>& fe)
{
	double Ji[3][3];
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	double*	gw = el.GaussWeights();
	
	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());

		// calculate the jacobian and its derivative
		// (and multiply by integration weight)
		double detJt = invjact(el, Ji, n)*gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.m_s;

		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);
		double* Gt = el.Gt(n);

		// --- first-order contribution ---
		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -=  (Gx*s.xx() + Gy*s.xy() + Gz*s.xz())*detJt;
			fe[3*i+1] -=  (Gx*s.xy() + Gy*s.yy() + Gz*s.yz())*detJt;
			fe[3*i+2] -=  (Gx*s.xz() + Gy*s.yz() + Gz*s.zz())*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::ElementInternalForce_QG(FESolidElement& el, vector<double>& fe)
{
	// get the nodal positions
	vec3d X[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	int neln = el.Nodes();
	for (int i=0; i<neln; ++i) X[i] = mesh.Node(el.m_node[i]).m_r0;

	// repeat for all integration points
	int nint = el.GaussPoints();
	double*	gw = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());

		// get the higher-order stress
		tens3drs& Q = mmpt2O.m_Q;

		// we'll evaluate this term in the material frame
		// so we need the Jacobian with respect to the reference configuration
		double J0 = detJ0(el, n);

		// loop over nodes
		for (int a=0; a<neln; ++a)
		{
			mat3d H;
			shape_gradient2(el, X, n, a, H);

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			for (int j=0; j<3; ++j)
				for (int k=0; k<3; ++k)
				{
					fe[3*a  ] -=  (Q(0,j,k)*H[j][k])*J0*gw[n];
					fe[3*a+1] -=  (Q(1,j,k)*H[j][k])*J0*gw[n];
					fe[3*a+2] -=  (Q(2,j,k)*H[j][k])*J0*gw[n];
				}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticMultiscaleDomain2O::Update()
{
	// call base class first
	// (this will call FEElasticMultiscaleDomain2O::UpdateElementStress)
	FEElasticSolidDomain::Update();

	// update internal surfaces
	UpdateInternalSurfaceStresses();

	// update the kinematic variables
	UpdateKinematics();
}

//-----------------------------------------------------------------------------
// Update some kinematical quantities needed for evaluating the discrete-Galerkin terms
void FEElasticMultiscaleDomain2O::UpdateKinematics()
{
	FEMesh& mesh = *GetMesh();

	// nodal displacements
	vec3d ut[FEElement::MAX_NODES];

	// shape function derivatives
	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];

	// loop over all facets
	int nd = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		// get the next facet
		FESurfaceElement& face = m_surf.Element(i);
		int nfn  = face.Nodes();
		int nint = face.GaussPoints();

		// calculate the displacement gradient jump across this facet
		for (int m=0; m<2; ++m)
		{
			// evaluate the spatial gradient of shape functions
			FESolidElement& el = Element(face.m_elem[m]);
			int neln = el.Nodes();

			// get the nodal displacements
			for (int j=0; j<neln; ++j)
			{
				FENode& node = mesh.Node(el.m_node[j]);
				ut[j] = node.m_rt - node.m_r0;
			}

			// loop over all integration points
			for (int n=0; n<nint; ++n)
			{
				// get the integration point data
				FEInternalSurface2O::Data& data = m_surf.GetData(nd + n);

				// evaluate element Jacobian and shape function derivatives
				// at this integration point
				vec3d& ksi = data.ksi[m];
				mat3d Ji;
				invjac0(el, ksi.x, ksi.y, ksi.z, Ji);
				el.shape_deriv(Gr, Gs, Gt, ksi.x, ksi.y, ksi.z);

				mat3d Gu; Gu.zero();
				for (int j=0; j<neln; ++j)
				{
					// calculate global gradient of shape functions
					// note that we need the transposed of Ji, not Ji itself !
					double Gx = Ji[0][0]*Gr[j]+Ji[1][0]*Gs[j]+Ji[2][0]*Gt[j];
					double Gy = Ji[0][1]*Gr[j]+Ji[1][1]*Gs[j]+Ji[2][1]*Gt[j];
					double Gz = Ji[0][2]*Gr[j]+Ji[1][2]*Gs[j]+Ji[2][2]*Gt[j];

					Gu[0][0] += ut[j].x*Gx; Gu[0][1] += ut[j].x*Gy; Gu[0][2] += ut[j].x*Gz;
					Gu[1][0] += ut[j].y*Gx; Gu[1][1] += ut[j].y*Gy; Gu[1][2] += ut[j].y*Gz;
					Gu[2][0] += ut[j].z*Gx; Gu[2][1] += ut[j].z*Gy; Gu[2][2] += ut[j].z*Gz;
				}

				if (m == 0) data.DgradU = Gu;
				else data.DgradU -= Gu;
			}
		}

		// don't forget to increment data counter
		nd += nint;
	}
}

//-----------------------------------------------------------------------------
// This function evaluates the stresses at either side of the internal surface
// facets.
void FEElasticMultiscaleDomain2O::UpdateInternalSurfaceStresses()
{
	FEMesh& mesh = *GetMesh();
	// calculate the material
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);

	// loop over all the internal surfaces
	int NF = m_surf.Elements(), nd = 0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& face = m_surf.Element(i);
		int nint = face.GaussPoints();
		for (int n=0; n<nint; ++n, ++nd)
		{
			FEInternalSurface2O::Data& data =  m_surf.GetData(nd);
			data.Qavg.zero();
			data.J0avg.zero();

			// get the deformation gradient and determinant
			for (int k=0; k<2; ++k)
			{
				FEMaterialPoint& mp = *data.m_pt[k];
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				FEMicroMaterialPoint2O& pt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();

				vec3d& ksi = data.ksi[k];

				// TODO: Is face.m_elem is a local index?
				FESolidElement& ek = Element(face.m_elem[k]);

				// evaluate deformation gradient, Jacobian and Hessian for this element
				pt.m_J = defgrad(ek, pt.m_F, ksi.x, ksi.y, ksi.z);
				defhess(ek, ksi.x, ksi.y, ksi.z, pt2O.m_G);

				// evaluate stresses at this integration point
				pmat->Stress2O(mp);

				data.Qavg += pt2O.m_Q*0.5;
				data.J0avg += pt2O.m_Ea;		// TODO: I only need to evaluate this on the first iteration
			}
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
		defhess(el, n, mmpt2O.m_G);

		pmat->Stress2O(mp);
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix
void FEElasticMultiscaleDomain2O::ElementStiffness(FEModel& fem, int iel, matrix& ke)
{
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);

	FESolidElement& el = Element(iel);

	// Get the current element's data
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	// Gx = dH/dx
	double Gx[FEElement::MAX_NODES];
	double Gy[FEElement::MAX_NODES];
	double Gz[FEElement::MAX_NODES];

	// H = d2N/dXdX
	mat3d G2[FEElement::MAX_NODES];

	// jacobian
	double Ji[3][3];

	// get the initial nodal coordinates
	vec3d X[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	mesh.GetInitialNodalCoordinates(el, X);
	
	// calculate element stiffness matrix
	const int nint = el.GaussPoints();
	const double *gw = el.GaussWeights();
	for (int ni=0; ni<nint; ++ni)
	{
		// get the material point data
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(ni);
		FEMicroMaterialPoint2O& mmpt2O = *(mp.ExtractData<FEMicroMaterialPoint2O>());

		// get the material tangents
		tens4ds C; C.zero();
		tens5ds L; L.zero();
		tens5ds H; H.zero();	// TODO: Why is this not calculated?
		tens6ds J; J.zero();
		pmat->Tangent2O(mp, C, L, J);

		// calculate jacobian
		double detJ0 = invjac0(el, Ji, ni)*gw[ni];

		// shape function derivatives
		double* Grn = el.Gr(ni);
		double* Gsn = el.Gs(ni);
		double* Gtn = el.Gt(ni);

		for (int i=0; i<neln; ++i)
		{
			double Gr = Grn[i];
			double Gs = Gsn[i];
			double Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;

			// second derivative of shape functions
			shape_gradient2(el, X, ni, i, G2[i]);
		}

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (int a=0; a<neln; ++a)
		{
			double Ga[3] = {Gx[a], Gy[a], Gz[a]};
			mat3d& Ha = G2[a];

			for (int b=0; b<neln; ++b)
			{
				double Gb[3] = {Gx[b], Gy[b], Gz[b]};
				mat3d& Hb = G2[b];

				mat3d Kab; Kab.zero();
				for (int j=0; j<3; ++j)
					for (int l=0; l<3; ++l)
					{
						Kab[0][0] += (Ga[j]*C(0,j,0,l)*Gb[j])*detJ0;
						Kab[0][1] += (Ga[j]*C(0,j,1,l)*Gb[j])*detJ0;
						Kab[0][2] += (Ga[j]*C(0,j,2,l)*Gb[j])*detJ0;

						Kab[1][0] += (Ga[j]*C(1,j,0,l)*Gb[j])*detJ0;
						Kab[1][1] += (Ga[j]*C(1,j,1,l)*Gb[j])*detJ0;
						Kab[1][2] += (Ga[j]*C(1,j,2,l)*Gb[j])*detJ0;

						Kab[2][0] += (Ga[j]*C(2,j,0,l)*Gb[j])*detJ0;
						Kab[2][1] += (Ga[j]*C(2,j,1,l)*Gb[j])*detJ0;
						Kab[2][2] += (Ga[j]*C(2,j,2,l)*Gb[j])*detJ0;

					}

				for (int j=0; j<3; ++j)
					for (int l=0; l<3; ++l)
						for (int m=0; m<3; ++m)
						{
							Kab[0][0] += (Ga[j]*L(0, j, 0, l, m)*Ha(l, m))*detJ0;
							Kab[0][1] += (Ga[j]*L(0, j, 1, l, m)*Ha(l, m))*detJ0;
							Kab[0][2] += (Ga[j]*L(0, j, 2, l, m)*Ha(l, m))*detJ0;

							Kab[1][0] += (Ga[j]*L(1, j, 0, l, m)*Ha(l, m))*detJ0;
							Kab[1][1] += (Ga[j]*L(1, j, 1, l, m)*Ha(l, m))*detJ0;
							Kab[1][2] += (Ga[j]*L(1, j, 2, l, m)*Ha(l, m))*detJ0;

							Kab[2][0] += (Ga[j]*L(2, j, 0, l, m)*Ha(l, m))*detJ0;
							Kab[2][1] += (Ga[j]*L(2, j, 1, l, m)*Ha(l, m))*detJ0;
							Kab[2][2] += (Ga[j]*L(2, j, 2, l, m)*Ha(l, m))*detJ0;
						}

				for (int j=0; j<3; ++j)
					for (int k=0; k<3; ++k)
						for (int m=0; m<3; ++m)
						{
							Kab[0][0] += (Ha(j,k)*H(0, j, k, 0, m)*Ga[m])*detJ0;
							Kab[0][1] += (Ha(j,k)*H(0, j, k, 1, m)*Ga[m])*detJ0;
							Kab[0][2] += (Ha(j,k)*H(0, j, k, 2, m)*Ga[m])*detJ0;

							Kab[1][0] += (Ha(j,k)*H(1, j, k, 0, m)*Ga[m])*detJ0;
							Kab[1][1] += (Ha(j,k)*H(1, j, k, 1, m)*Ga[m])*detJ0;
							Kab[1][2] += (Ha(j,k)*H(1, j, k, 2, m)*Ga[m])*detJ0;

							Kab[2][0] += (Ha(j,k)*H(2, j, k, 0, m)*Ga[m])*detJ0;
							Kab[2][1] += (Ha(j,k)*H(2, j, k, 1, m)*Ga[m])*detJ0;
							Kab[2][2] += (Ha(j,k)*H(2, j, k, 2, m)*Ga[m])*detJ0;
						}

				for (int j=0; j<3; ++j)
					for (int k=0; k<3; ++k)
						for (int m=0; m<3; ++m)
							for (int n=0; n<3; ++n)
							{
								Kab[0][0] += (Ha(j,k)*J(0, j, k, 0, m, n)*Hb(m,n))*detJ0;
								Kab[0][1] += (Ha(j,k)*J(0, j, k, 1, m, n)*Hb(m,n))*detJ0;
								Kab[0][2] += (Ha(j,k)*J(0, j, k, 2, m, n)*Hb(m,n))*detJ0;

								Kab[1][0] += (Ha(j,k)*J(1, j, k, 0, m, n)*Hb(m,n))*detJ0;
								Kab[1][1] += (Ha(j,k)*J(1, j, k, 1, m, n)*Hb(m,n))*detJ0;
								Kab[1][2] += (Ha(j,k)*J(1, j, k, 2, m, n)*Hb(m,n))*detJ0;

								Kab[2][0] += (Ha(j,k)*J(2, j, k, 0, m, n)*Hb(m,n))*detJ0;
								Kab[2][1] += (Ha(j,k)*J(2, j, k, 1, m, n)*Hb(m,n))*detJ0;
								Kab[2][2] += (Ha(j,k)*J(2, j, k, 2, m, n)*Hb(m,n))*detJ0;
							}

				// add it to the element matrix
				ke.add(3*a, 3*b, Kab);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the gradient of deformation gradient of element el at integration point n.
//! The gradient of the deformation gradient is returned in G.
void FEElasticMultiscaleDomain2O::defhess(FESolidElement &el, int n, tens3drs &G)
{
	int neln = el.Nodes();

	// get the nodal positions
	vec3d X[FEElement::MAX_NODES];
	vec3d x[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	for (int i=0; i<neln; ++i)
	{
		X[i] = mesh.Node(el.m_node[i]).m_r0;
		x[i] = mesh.Node(el.m_node[i]).m_rt;
	}

	// loop over nodes
	G.zero();
	for (int a=0; a<neln; ++a)
	{
		mat3d H;
		shape_gradient2(el, X, n, a, H);

		// calculate gradient of deformation gradient
		// Note that k >= j. Since tensdrs has symmetries this
		// prevents overwriting of symmetric components
		for (int j=0; j<3; ++j)
			for (int k=j; k<3; ++k)
			{
				G(0,j,k) += H[j][k]*x[a].x;
				G(1,j,k) += H[j][k]*x[a].y;
				G(2,j,k) += H[j][k]*x[a].z;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the gradient of deformation gradient of element el at integration point n.
//! The gradient of the deformation gradient is returned in G.
void FEElasticMultiscaleDomain2O::defhess(FESolidElement &el, double r, double s, double t, tens3drs &G)
{
	int neln = el.Nodes();

	// get the nodal positions
	vec3d X[FEElement::MAX_NODES];
	vec3d x[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	for (int i=0; i<neln; ++i)
	{
		X[i] = mesh.Node(el.m_node[i]).m_r0;
		x[i] = mesh.Node(el.m_node[i]).m_rt;
	}

	// loop over nodes
	G.zero();
	for (int a=0; a<neln; ++a)
	{
		mat3d H;
		shape_gradient2(el, X, r, s, t, a, H);

		// calculate gradient of deformation gradient
		// Note that k >= j. Since tensdrs has symmetries this
		// prevents overwriting of symmetric components
		for (int j=0; j<3; ++j)
			for (int k=j; k<3; ++k)
			{
				G(0,j,k) += H[j][k]*x[a].x;
				G(1,j,k) += H[j][k]*x[a].y;
				G(2,j,k) += H[j][k]*x[a].z;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the second derivative of shape function N[node] with respect
//! to the material coordinates.
void FEElasticMultiscaleDomain2O::shape_gradient2(const FESolidElement& el, vec3d* X, int n, int node, mat3d& H)
{
	int neln = el.Nodes();
	double*	gw = el.GaussWeights();

	// inverse of Jacobian matrix
	double Ji[3][3];

	// we'll evaluate this term in the material frame
	// so we need the Jacobian with respect to the reference configuration
	invjac0(el, Ji, n);

	// shape function derivatives
	double* Gr = el.Gr(n);
	double* Gs = el.Gs(n);
	double* Gt = el.Gt(n);

	double *Grrn = el.Grr(n); double *Grsn = el.Grs(n); double *Grtn = el.Grt(n);
	double *Gsrn = el.Gsr(n); double *Gssn = el.Gss(n); double *Gstn = el.Gst(n);
	double *Gtrn = el.Gtr(n); double *Gtsn = el.Gts(n); double *Gttn = el.Gtt(n);

	// calculate K = dJ/dr
	double K[3][3][3] = {0};
	for (int a=0; a<neln; ++a)
	{
		// second derivatives of shape functions
		double G2[3][3];
		G2[0][0] = Grrn[a]; G2[0][1] = Grsn[a]; G2[0][2] = Grtn[a];
		G2[1][0] = Gsrn[a]; G2[1][1] = Gssn[a]; G2[1][2] = Gstn[a];
		G2[2][0] = Gtrn[a]; G2[2][1] = Gtsn[a]; G2[2][2] = Gttn[a];

		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
			{
				K[0][j][k] += G2[j][k]*X[a].x;
				K[1][j][k] += G2[j][k]*X[a].y;
				K[2][j][k] += G2[j][k]*X[a].z;
			}
	}

	// calculate A = -J^-1*dJ/drJ^-1
	double A[3][3][3] = {0};
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			for (int p=0; p<3; ++p)
				for (int q=0; q<3; ++q)
				{
					A[i][j][0] -= Ji[j][p]*K[p][q][0]*Ji[q][i];
					A[i][j][1] -= Ji[j][p]*K[p][q][1]*Ji[q][i];
					A[i][j][2] -= Ji[j][p]*K[p][q][2]*Ji[q][i];
				}
		}


	// first derivative of shape functions
	double G1[3];
	G1[0] = Gr[node];
	G1[1] = Gs[node];
	G1[2] = Gt[node];

	// second derivatives of shape functions
	double G2[3][3];
	G2[0][0] = Grrn[node]; G2[0][1] = Grsn[node]; G2[0][2] = Grtn[node];
	G2[1][0] = Gsrn[node]; G2[1][1] = Gssn[node]; G2[1][2] = Gstn[node];
	G2[2][0] = Gtrn[node]; G2[2][1] = Gtsn[node]; G2[2][2] = Gttn[node];

	// calculate dB/dr
	double D[3][3] = {0};
	for (int i=0; i<3; ++i)
		for (int k=0; k<3; ++k)
		{
			for (int j=0; j<3; ++j) D[i][k] += A[i][j][k]*G1[j] + Ji[j][i]*G2[j][k];
		}

	// calculate global gradient of shape functions
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			H[i][j] = D[i][0]*Ji[0][j] + D[i][1]*Ji[1][j] + D[i][2]*Ji[2][j];
		}
}

//-----------------------------------------------------------------------------
//! Calculates the second derivative of shape function N[node] with respect
//! to the material coordinates.
void FEElasticMultiscaleDomain2O::shape_gradient2(const FESolidElement& el, vec3d* X, double r, double s, double t, int node, mat3d& H)
{
	int neln = el.Nodes();

	// we need the Jacobian with respect to the reference configuration
	double Ji[3][3];
	invjac0(el, Ji, r, s, t);

	// shape function derivatives
	const int M = FEElement::MAX_NODES;
	double Gr[M], Gs[M], Gt[M];
	el.shape_deriv(Gr, Gs, Gt, r, s, t);

	double Grr[M], Gss[M], Gtt[M], Grs[M], Gst[M], Grt[M];
	el.shape_deriv2(Grr, Gss, Gtt, Grs, Gst, Grt, r, s, t);

	// calculate K = dJ/dr
	double K[3][3][3] = {0};
	for (int a=0; a<neln; ++a)
	{
		// second derivatives of shape functions
		double G2[3][3];
		G2[0][0] = Grr[a]; G2[0][1] = Grs[a]; G2[0][2] = Grt[a];
		G2[1][0] = Grs[a]; G2[1][1] = Gss[a]; G2[1][2] = Gst[a];
		G2[2][0] = Grt[a]; G2[2][1] = Gst[a]; G2[2][2] = Gtt[a];

		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
			{
				K[0][j][k] += G2[j][k]*X[a].x;
				K[1][j][k] += G2[j][k]*X[a].y;
				K[2][j][k] += G2[j][k]*X[a].z;
			}
	}

	// calculate A = -J^-1*dJ/drJ^-1
	double A[3][3][3] = {0};
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			for (int p=0; p<3; ++p)
				for (int q=0; q<3; ++q)
				{
					A[i][j][0] -= Ji[j][p]*K[p][q][0]*Ji[q][i];
					A[i][j][1] -= Ji[j][p]*K[p][q][1]*Ji[q][i];
					A[i][j][2] -= Ji[j][p]*K[p][q][2]*Ji[q][i];
				}
		}

	// first derivative of shape functions
	double G1[3];
	G1[0] = Gr[node];
	G1[1] = Gs[node];
	G1[2] = Gt[node];

	// second derivatives of shape functions
	double G2[3][3];
	G2[0][0] = Grr[node]; G2[0][1] = Grs[node]; G2[0][2] = Grt[node];
	G2[1][0] = Grs[node]; G2[1][1] = Gss[node]; G2[1][2] = Gst[node];
	G2[2][0] = Grt[node]; G2[2][1] = Gst[node]; G2[2][2] = Gtt[node];

	// calculate dB/dr
	double D[3][3] = {0};
	for (int i=0; i<3; ++i)
		for (int k=0; k<3; ++k)
		{
			for (int j=0; j<3; ++j) D[i][k] += A[i][j][k]*G1[j] + Ji[j][i]*G2[j][k];
		}

	// calculate global gradient of shape functions
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			H[i][j] += D[i][0]*Ji[0][j] + D[i][1]*Ji[1][j] + D[i][2]*Ji[2][j];
		}
}
