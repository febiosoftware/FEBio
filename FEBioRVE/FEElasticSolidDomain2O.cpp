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
#include "FEElasticSolidDomain2O.h"
#include "FEElasticMaterial2O.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// helper function for comparing two facets
bool compare_facets(int* na, int* nb, int nodes)
{
	switch (nodes)
	{
	case 3: // 3-node triangle
	case 6: // 6-node triangle
	case 7: // 7-node triangle
		if ((na[0]!=nb[0])&&(na[0]!=nb[1])&&(na[0]!=nb[2])) return false;
		if ((na[1]!=nb[0])&&(na[1]!=nb[1])&&(na[1]!=nb[2])) return false;
		if ((na[2]!=nb[0])&&(na[2]!=nb[1])&&(na[2]!=nb[2])) return false;
		break;
	case 4: // 4-node quad
	case 8: // 8-node quad
	case 9: // 9-node quad
		if ((na[0]!=nb[0])&&(na[0]!=nb[1])&&(na[0]!=nb[2])&&(na[0]!=nb[3])) return false;
		if ((na[1]!=nb[0])&&(na[1]!=nb[1])&&(na[1]!=nb[2])&&(na[1]!=nb[3])) return false;
		if ((na[2]!=nb[0])&&(na[2]!=nb[1])&&(na[2]!=nb[2])&&(na[2]!=nb[3])) return false;
		if ((na[3]!=nb[0])&&(na[3]!=nb[1])&&(na[3]!=nb[2])&&(na[3]!=nb[3])) return false;
		break;
	default:
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
inline bool compare_tri(int* na, int* nb, int i, int j, int k)
{
	return ((na[0] == nb[i])&&(na[1] == nb[j])&&(na[2] == nb[k]));
}

//-----------------------------------------------------------------------------
inline bool compare_quad(int* na, int* nb, int i, int j, int k, int l)
{
	return ((na[0] == nb[i])&&(na[1] == nb[j])&&(na[2] == nb[k])&&(na[3] == nb[l]));
}

//-----------------------------------------------------------------------------
// helper function for mapping surface facet coordinates to element facet coordinates
vec2d map_facet_to_facet(int* na, int* nb, int nodes, double r, double s)
{
	switch (nodes)
	{
	// triangles
	case 3:
	case 6:
	case 7:
		if (compare_tri(na, nb, 0, 1, 2)) return vec2d(  r  ,     s);
		if (compare_tri(na, nb, 2, 0, 1)) return vec2d(    s, 1-r-s);
		if (compare_tri(na, nb, 1, 2, 0)) return vec2d(1-r-s,   r  );
		if (compare_tri(na, nb, 0, 2, 1)) return vec2d(    s,   r  );
		if (compare_tri(na, nb, 2, 1, 0)) return vec2d(  r  , 1-r-s);
		if (compare_tri(na, nb, 1, 0, 2)) return vec2d(1-r-s,     s);
		break;
	// quads
	case 4:
	case 8:
	case 9:
		if (compare_quad(na, nb, 0, 1, 2, 3)) return vec2d(  r,   s);
		if (compare_quad(na, nb, 3, 0, 1, 2)) return vec2d(  s,  -r);
		if (compare_quad(na, nb, 2, 3, 0, 1)) return vec2d( -r,  -s);
		if (compare_quad(na, nb, 1, 2, 3, 0)) return vec2d( -s,   r);
		if (compare_quad(na, nb, 2, 1, 0, 3)) return vec2d( -s,  -r);
		if (compare_quad(na, nb, 3, 2, 1, 0)) return vec2d(  r,  -s);
		if (compare_quad(na, nb, 0, 3, 2, 1)) return vec2d(  s,   r);
		if (compare_quad(na, nb, 1, 0, 3, 2)) return vec2d( -r,   s);
		break;
	}

	// we shouldn't get here
	assert(false);
	return vec2d(0,0);
}

//-----------------------------------------------------------------------------
// Notice that this depends on how the facet nodes are numbered (see FEElement::GetFace) 
vec3d map_facet_to_volume_coordinates_tet(int nface, const vec2d& q)
{
	double h1 = q.x(), h2 = q.y(), h3 = 1.0 - h1 - h2;
	double g1, g2, g3;
	switch (nface)
	{
	case 0: g1 = h1; g2 = 0.; g3 = h2; break;
	case 1: g1 = h3; g2 = h1; g3 = h2; break;
	case 2: g1 = 0.; g2 = h3; g3 = h2; break;
	case 3: g1 = h1; g2 = h3; g3 = 0.; break;
	default:
		assert(false);
	}
	return vec3d(g1, g2, g3);
}

//-----------------------------------------------------------------------------
// Notice that this depends on how the facet nodes are numbered (see FEElement::GetFace) 
vec3d map_facet_to_volume_coordinates_hex(int nface, const vec2d& q)
{
	double h1 = q.x(), h2 = q.y();
	double g1, g2, g3;
	switch (nface)
	{
	case 0: g1 =  h1; g2 = -1.; g3 =  h2; break;
	case 1: g1 =  1.; g2 =  h1; g3 =  h2; break;
	case 2: g1 = -h1; g2 =  1.; g3 =  h2; break;
	case 3: g1 = -1.; g2 = -h1; g3 =  h2; break;
	case 4: g1 =  h2; g2 =  h1; g3 = -1.; break;
	case 5: g1 =  h1; g2 =  h2; g3 =  1.; break;
	default:
		assert(false);
	}
	return vec3d(g1, g2, g3);
}

//-----------------------------------------------------------------------------
// helper function that maps natural coordinates from facets to elements
vec3d map_facet_to_solid(FEMesh& mesh, FESurfaceElement& face, FESolidElement& el, double r, double s)
{
	int fn[FEElement::MAX_NODES];
	int nfaces = el.Faces();
	for (int i=0; i<nfaces; ++i)
	{
		el.GetFace(i, fn);
		if (compare_facets(&face.m_node[0], fn, face.Nodes()))
		{
			// map the facet coordinates to the element's facet coordinates
			// (faces can be rotated or inverted w.r.t. the element's face)
			vec2d b = map_facet_to_facet(&face.m_node[0], fn, face.Nodes(), r, s);

			// convert facet coordinates to volume coordinates
			switch (el.Nodes())
			{
			// tets
			case 4:
			case 10:
			case 15: return map_facet_to_volume_coordinates_tet(i, b); break;
			// hexes
			case 8:
			case 20:
			case 27: return map_facet_to_volume_coordinates_hex(i, b); break;
			}

			assert(false);
		}
	}

	// we shouldn't get here
	assert(false);
	return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
FEElasticSolidDomain2O::FEInternalSurface2O::FEInternalSurface2O()
{
	m_h = 1.0;
}

bool FEElasticSolidDomain2O::FEInternalSurface2O::Initialize(FEElasticSolidDomain2O* dom)
{
	// get the material
	// We need this for allocating material point data on the internal facets
	FEMaterial* pmat = dom->GetMaterial();
	assert(pmat);
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
			m_data[nnf].m_pt[0] = new FEMaterialPoint(pmat->CreateMaterialPointData());
			m_data[nnf].m_pt[1] = new FEMaterialPoint(pmat->CreateMaterialPointData());

			m_data[nnf].m_pt[0]->Init();
			m_data[nnf].m_pt[1]->Init();

			// iso-parametric coordinates in surface element
			double h1 = face.gr(n);
			double h2 = face.gs(n);

			for (int k=0; k<2; ++k)
			{
				// get the adjacent solid element
				FESolidElement& ek = static_cast<FESolidElement&>(*face.m_elem[k]);

				// map the iso-parametric coordinates from the facet to the solid element
				m_data[nnf].ksi[k] = map_facet_to_solid(mesh, face, ek, h1, h2);
			}

#ifdef _DEBUG
			// This checks that the two integration points coincide physically at the same point
			FESolidElement& ea = static_cast<FESolidElement&>(*face.m_elem[0]);
			FESolidElement& eb = static_cast<FESolidElement&>(*face.m_elem[1]);
			vec3d xa[FEElement::MAX_NODES];
			vec3d xb[FEElement::MAX_NODES];
			for (int a=0; a<ea.Nodes(); a++) xa[a] = mesh.Node(ea.m_node[a]).m_r0;
			for (int b=0; b<eb.Nodes(); b++) xb[b] = mesh.Node(eb.m_node[b]).m_r0;
			vec3d ksia = m_data[nnf].ksi[0];
			vec3d ksib = m_data[nnf].ksi[1];

			vec3d ra = ea.evaluate(xa, ksia.x, ksia.y, ksia.z);
			vec3d rb = eb.evaluate(xb, ksib.x, ksib.y, ksib.z);
			vec3d dr = ra - rb;
			double Dr = dr.norm();
			assert(Dr <1e-12);
#endif
		}
	}

	// calculate the element size
	m_h = 0.0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = m_ps->Element(i);
		int neln = el.Nodes();
		FEBoundingBox box(mesh.Node(el.m_node[0]).m_r0);
		for (int j=1; j<neln; ++j) box.add(mesh.Node(el.m_node[j]).m_r0);

		double R = 0.5*box.radius();
		if (R > m_h) m_h = R;
	}

	return true;
}

//-----------------------------------------------------------------------------
FEElasticSolidDomain2O::FEElasticSolidDomain2O(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
	m_binitJ0 = false;
}

//-----------------------------------------------------------------------------
//! Initialize element data
bool FEElasticSolidDomain2O::Init()
{
	// initialize base class first
	if (FEElasticSolidDomain::Init() == false) return false;

	// initialize the internal surface data
	if (m_surf.Initialize(this) == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// call base class first
	FEElasticSolidDomain::BuildMatrixProfile(M);

	// do the interface elements
	FEMesh& mesh = *GetMesh();
	vector<int> lm;

	// loop over all internal facets
	int NF = m_surf.Elements();
	for (int n=0; n<NF; ++n)
	{
		// get the next surface element
		FESurfaceElement& el = m_surf.Element(n);

		// get the solid elements from this interface element
		FESolidElement& ela = static_cast<FESolidElement&>(*el.m_elem[1]);
		FESolidElement& elb = static_cast<FESolidElement&>(*el.m_elem[0]);
		int nelna = ela.Nodes();
		int nelnb = elb.Nodes();

		// setup the LM vector
		// TODO: this assumes that the X,Y,Z degrees of freedom are 0,1,2 respectively
		int ndof = 3*(nelna + nelnb);
		lm.resize(ndof);
		for (int i=0; i<nelna; ++i)
		{
			vector<int>& id = mesh.Node(ela.m_node[i]).m_ID;
			lm[3*i  ] = id[0];
			lm[3*i+1] = id[1];
			lm[3*i+2] = id[2];
		}
		for (int i=0; i<nelnb; ++i)
		{
			vector<int>& id = mesh.Node(elb.m_node[i]).m_ID;
			lm[3*(nelna+i)  ] = id[0];
			lm[3*(nelna+i)+1] = id[1];
			lm[3*(nelna+i)+2] = id[2];
		}

		// add it to the profile
		M.build_add(lm);
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	FEElasticSolidDomain::PreSolveUpdate(timeInfo);

	int NF = m_surf.Elements(), nd = 0;
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = m_surf.Element(i);
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n, ++nd)
		{
			FEInternalSurface2O::Data& data = m_surf.GetData(nd);
			data.m_pt[0]->Update(timeInfo);
			data.m_pt[1]->Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::Update(const FETimeInfo& tp)
{
	// call base class first
	// (this will call FEElasticMultiscaleDomain2O::UpdateElementStress)
	FEElasticSolidDomain::Update(tp);

	// update internal surfaces
	UpdateInternalSurfaceStresses();

	// update the kinematic variables
	UpdateKinematics();
}

//-----------------------------------------------------------------------------
// Update some kinematical quantities needed for evaluating the discontinuous-Galerkin terms
void FEElasticSolidDomain2O::UpdateKinematics()
{
	FEMesh& mesh = *GetMesh();

	// nodal displacements
	vec3d ut[FEElement::MAX_NODES];

	// shape function derivatives
	vec3d G[FEElement::MAX_NODES];

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
			FESolidElement& el = static_cast<FESolidElement&>(*face.m_elem[m]);
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
				ShapeGradient0(el, ksi.x, ksi.y, ksi.z, G);

				// calculate displacement gradient
				mat3d Gu; Gu.zero();
				for (int j=0; j<neln; ++j) Gu += ut[j] & G[j];

				if (m == 0) data.DgradU = Gu;
				else data.DgradU = Gu - data.DgradU;
			}
		}

		// don't forget to increment data counter
		nd += nint;
	}
}

//-----------------------------------------------------------------------------
// This function evaluates the stresses at either side of the internal surface
// facets.
void FEElasticSolidDomain2O::UpdateInternalSurfaceStresses()
{
	FEMesh& mesh = *GetMesh();
	// calculate the material
	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(m_pMat);

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
			
			if (m_binitJ0 == false) data.J0avg.zero();

			// get the deformation gradient and determinant
			for (int k=0; k<2; ++k)
			{
				FEMaterialPoint& mp = *data.m_pt[k];
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				FEElasticMaterialPoint2O& pt2O = *mp.ExtractData<FEElasticMaterialPoint2O>();

				vec3d& ksi = data.ksi[k];

				// TODO: Is face.m_elem is a local index?
				FESolidElement& ek = static_cast<FESolidElement&>(*face.m_elem[k]);

				// evaluate deformation gradient, Jacobian and Hessian for this element
				pt.m_J = defgrad(ek, pt.m_F, ksi.x, ksi.y, ksi.z);
				defhess(ek, ksi.x, ksi.y, ksi.z, pt2O.m_G);

				// evaluate stresses at this integration point
				pmat->Stress(mp, pt2O.m_PK1, pt2O.m_Q);

				data.Qavg += pt2O.m_Q*0.5;

				// TODO: Can I evaluate this elsewhere?
				tens4d C;
				tens5d L, H;
				tens6d J;
				pmat->Tangent(mp, C, L, H, J);
				data.H[k] = H;
				data.J[k] = J;

				if (m_binitJ0 == false)
				{
					// we need to evaluate the initial stiffnesses as well
					data.J0[k] = J;
					data.H0[k] = H;
					data.J0avg += J*0.5;
				}
			}
		}
	}

	// set flag indicating J0 has been initialized
	if (pmat->m_buseJ0) m_binitJ0 = true;
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FEElasticSolidDomain2O::UpdateElementStress(int iel, const FETimeInfo& tp)
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
	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(m_pMat);

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEElasticMaterialPoint2O& pt2O = *(mp.ExtractData<FEElasticMaterialPoint2O>());
		
		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		mp.m_r0 = el.Evaluate(r0, n);
		mp.m_rt = el.Evaluate(rt, n);

		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
		defhess(el, n, pt2O.m_G);

		// evaluate stresses
		pmat->Stress(mp, pt2O.m_PK1, pt2O.m_Q);

		// store the Cauchy stress as well
		double J = pt.m_J;
		const mat3d Ft = pt.m_F.transpose();
		pt.m_s = (pt2O.m_PK1*Ft).sym()/J;	
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::InternalForces(FEGlobalVector& R)
{
	// call base class first
	FEElasticSolidDomain::InternalForces(R);

	// add the discontinuous-Galerkin contribution
	InternalForcesDG1(R);
	InternalForcesDG2(R);
	InternalForcesDG3(R);
}

//-----------------------------------------------------------------------------
//! Evaluate contribution of discontinuous-Galerkin enforcement of stress flux.
void FEElasticSolidDomain2O::InternalForcesDG1(FEGlobalVector& R)
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
			FESolidElement& el = static_cast<FESolidElement&>(*face.m_elem[m]);
			int neln = el.Nodes();

			// sign
			double sgn = (m==0 ? -1.0 : 1.0);

			// allocate force vector
			int ndof = neln*3;
			fe.assign(ndof, 0.0);

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
void FEElasticSolidDomain2O::InternalForcesDG2(FEGlobalVector& R)
{
	FEMesh& mesh = *GetMesh();
	FESurface& surf = *m_surf.GetSurface();

	vec3d X[FEElement::MAX_NODES];
	vec3d G[FEElement::MAX_NODES];
	mat3d H[FEElement::MAX_NODES];

	vector<double> fe;
	vector<int> lm;

	// loop over all internal surface elements
	int nd = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		// get the next surface element
		FESurfaceElement& face = m_surf.Element(i);
		int nfn = face.Nodes();
		int nint = face.GaussPoints();

		// loop over both sides
		for (int m=0; m<2; ++m)
		{
			FESolidElement& el = static_cast<FESolidElement&>(*face.m_elem[m]);
			int neln = el.Nodes();

			// get the initial nodal positions
			for (int j=0; j<neln; ++j) X[j] = mesh.Node(el.m_node[j]).m_r0;

			// allocate force vector
			int ndof = neln*3;
			fe.assign(ndof, 0.0);

			// loop over integration points
			double* gw = face.GaussWeights();
			for  (int n=0; n<nint; ++n)
			{
				// get facet data
				FEInternalSurface2O::Data& data = m_surf.GetData(nd + n);
				const tens6d& J0 = data.J0[m];
				const tens5d& H0 = data.H0[m];
				const mat3d&  Du = data.DgradU;

				// calculate jacobian and surface normal
				vec3d nu(0,0,0);
				double J = surf.jac0(face, n, nu);
				double Nu[3] = {nu.x, nu.y, nu.z};

				// evaluate element shape function derivatives
				// at this integration point
				vec3d& ksi = data.ksi[m];
				ShapeGradient0(el, ksi.x, ksi.y, ksi.z, G);
				shape_gradient2(el, X, ksi.x, ksi.y, ksi.z, H);

				for (int a=0; a<neln; ++a)
				{
					double Ax[3][3][3] = {0};
					double Ay[3][3][3] = {0};
					double Az[3][3][3] = {0};
					for (int j=0; j<3; ++j)
						for (int k=0; k<3; ++k)
						{
							for (int q=0; q<3; ++q)
							{
								for (int r=0; r<3; ++r)
								{
									Ax[0][j][k] += J0(0, j, k, 0, q, r)*H[a](q ,r);
									Ax[1][j][k] += J0(1, j, k, 0, q, r)*H[a](q ,r);
									Ax[2][j][k] += J0(2, j, k, 0, q, r)*H[a](q ,r);

									Ay[0][j][k] += J0(0, j, k, 1, q, r)*H[a](q ,r);
									Ay[1][j][k] += J0(1, j, k, 1, q, r)*H[a](q ,r);
									Ay[2][j][k] += J0(2, j, k, 1, q, r)*H[a](q ,r);

									Az[0][j][k] += J0(0, j, k, 2, q, r)*H[a](q ,r);
									Az[1][j][k] += J0(1, j, k, 2, q, r)*H[a](q ,r);
									Az[2][j][k] += J0(2, j, k, 2, q, r)*H[a](q ,r);
								}
							}

							Ax[0][j][k] += H0(0, j, k, 0, 0)*G[a].x + H0(0, j, k, 0, 1)*G[a].y + H0(0, j, k, 0, 2)*G[a].z;
							Ax[1][j][k] += H0(1, j, k, 0, 0)*G[a].x + H0(1, j, k, 0, 1)*G[a].y + H0(1, j, k, 0, 2)*G[a].z;
							Ax[2][j][k] += H0(2, j, k, 0, 0)*G[a].x + H0(2, j, k, 0, 1)*G[a].y + H0(2, j, k, 0, 2)*G[a].z;

							Ay[0][j][k] += H0(0, j, k, 1, 0)*G[a].x + H0(0, j, k, 1, 1)*G[a].y + H0(0, j, k, 1, 2)*G[a].z;
							Ay[1][j][k] += H0(1, j, k, 1, 0)*G[a].x + H0(1, j, k, 1, 1)*G[a].y + H0(1, j, k, 1, 2)*G[a].z;
							Ay[2][j][k] += H0(2, j, k, 1, 0)*G[a].x + H0(2, j, k, 1, 1)*G[a].y + H0(2, j, k, 1, 2)*G[a].z;

							Az[0][j][k] += H0(0, j, k, 2, 0)*G[a].x + H0(0, j, k, 2, 1)*G[a].y + H0(0, j, k, 2, 2)*G[a].z;
							Az[1][j][k] += H0(1, j, k, 2, 0)*G[a].x + H0(1, j, k, 2, 1)*G[a].y + H0(1, j, k, 2, 2)*G[a].z;
							Az[2][j][k] += H0(2, j, k, 2, 0)*G[a].x + H0(2, j, k, 2, 1)*G[a].y + H0(2, j, k, 2, 2)*G[a].z;
						}

						double fa[3] = {0};
						for (int j=0; j<3; ++j)
							for (int k=0; k<3; ++k)
							{
								fa[0] += Du(0,j)*Ax[0][j][k]*Nu[k] + Du(1,j)*Ax[1][j][k]*Nu[k] + Du(2,j)*Ax[2][j][k]*Nu[k];
								fa[1] += Du(0,j)*Ay[0][j][k]*Nu[k] + Du(1,j)*Ay[1][j][k]*Nu[k] + Du(2,j)*Ay[2][j][k]*Nu[k];
								fa[2] += Du(0,j)*Az[0][j][k]*Nu[k] + Du(1,j)*Az[1][j][k]*Nu[k] + Du(2,j)*Az[2][j][k]*Nu[k];
							}

					// the negative sign is because we need to subtract the internal forces
					// from the residual
					fe[3*a  ] -= fa[0]*J*gw[n]*0.5;
					fe[3*a+1] -= fa[1]*J*gw[n]*0.5;
					fe[3*a+2] -= fa[2]*J*gw[n]*0.5;
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
//! Evaluate contribution of discontinuous-Galerkin enforcement of stress flux.
void FEElasticSolidDomain2O::InternalForcesDG3(FEGlobalVector& R)
{
	FEMesh& mesh = *GetMesh();
	FESurface& surf = *m_surf.GetSurface();

	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(GetMaterial());
	assert(pmat);
	double beta = pmat->m_beta;
	double h = m_surf.GetElementSize();

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
			FESolidElement& el = static_cast<FESolidElement&>(*face.m_elem[m]);
			int neln = el.Nodes();

			// sign
			double sgn = (m==0 ? -1.0 : 1.0);

			// allocate force vector
			int ndof = neln*3;
			fe.assign(ndof, 0.0);

			// loop over all the integration points
			double* gw = face.GaussWeights();
			for (int n=0; n<nint; ++n)
			{
				// get facet data
				FEInternalSurface2O::Data& data = m_surf.GetData(nd + n);

				// average stiffness across interface
				tens6d& J0avg = data.J0avg;

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
					fe[3*j  ] -= f[0]*gw[n]*J*sgn*beta/h;
					fe[3*j+1] -= f[1]*gw[n]*J*sgn*beta/h;
					fe[3*j+2] -= f[2]*gw[n]*J*sgn*beta/h;
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
void FEElasticSolidDomain2O::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	// contribution from stress
	ElementInternalForce_PF(el, fe);

	// contriubtion from higher-order stress
	ElementInternalForce_QG(el, fe);
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::ElementInternalForce_PF(FESolidElement& el, vector<double>& fe)
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
void FEElasticSolidDomain2O::ElementInternalForce_QG(FESolidElement& el, vector<double>& fe)
{
	// get the nodal positions
	vec3d X[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	int neln = el.Nodes();
	for (int i=0; i<neln; ++i) X[i] = mesh.Node(el.m_node[i]).m_r0;

	mat3d H[FEElement::MAX_NODES];

	// repeat for all integration points
	int nint = el.GaussPoints();
	double*	gw = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint2O& pt2O = *(mp.ExtractData<FEElasticMaterialPoint2O>());

		// get the higher-order stress
		const tens3drs& Q = pt2O.m_Q;

		// we'll evaluate this term in the material frame
		// so we need the Jacobian with respect to the reference configuration
		double J0 = detJ0(el, n);

		// shape function derivatives
		shape_gradient2(el, X, n, H);

		// loop over nodes
		for (int a=0; a<neln; ++a)
		{
			mat3d& Ha = H[a];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			for (int j=0; j<3; ++j)
				for (int k=0; k<3; ++k)
				{
					fe[3*a  ] -=  (Q(0,j,k)*Ha[j][k])*J0*gw[n];
					fe[3*a+1] -=  (Q(1,j,k)*Ha[j][k])*J0*gw[n];
					fe[3*a+2] -=  (Q(2,j,k)*Ha[j][k])*J0*gw[n];
				}
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::StiffnessMatrix(FELinearSystem& LS)
{
	// repeat over all solid elements
	int NE = (int)m_Elem.size();
	FETimeInfo tp = GetFEModel()->GetTime();
	
	#pragma omp parallel for shared (NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);
		
		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate element stiffness
		ElementStiffness(tp, iel, ke);

		// get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}

	// stiffness matrix from discontinuous Galerkin
	StiffnessMatrixDG(LS);
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::StiffnessMatrixDG(FELinearSystem& LS)
{
	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(GetMaterial());
	assert(pmat);

	// get stiffness flags
	bool bKDG1 = pmat->m_bKDG1;
	bool bKDG2 = pmat->m_bKDG2;
	bool bKDG3 = pmat->m_bKDG3;
	if ((bKDG1==false)&&(bKDG2==false)&&(bKDG3==false)) return;

	FEElementMatrix ke;
	vector<int> lm;
	vector<int> en;

	FEMesh& mesh = *GetMesh();

	// loop over all internal facets
	int nd = 0;
	int NF = m_surf.Elements();
	for (int n=0; n<NF; ++n)
	{
		// get the next surface element
		FESurfaceElement& el = m_surf.Element(n);

		// get the solid elements from this interface element
		FESolidElement& ela = static_cast<FESolidElement&>(*el.m_elem[1]);
		FESolidElement& elb = static_cast<FESolidElement&>(*el.m_elem[0]);
		int nelna = ela.Nodes();
		int nelnb = elb.Nodes();

		// init element stiffness
		int ndof = 3*(nelna + nelnb);
		ke.resize(ndof, ndof);
		ke.zero();

		// get the element stiffness matrix
		if (bKDG1) ElementStiffnessMatrixDG1(el, &m_surf.GetData(nd), ke);
		if (bKDG2) ElementStiffnessMatrixDG2(el, &m_surf.GetData(nd), ke);
		if (bKDG3) ElementStiffnessMatrixDG3(el, &m_surf.GetData(nd), ke);

		// setup the LM vector
		// TODO: this assumes that the X,Y,Z degrees of freedom are 0,1,2 respectively
		lm.resize(ndof);
		for (int i=0; i<nelna; ++i)
		{
			vector<int>& id = mesh.Node(ela.m_node[i]).m_ID;
			lm[3*i  ] = id[0];
			lm[3*i+1] = id[1];
			lm[3*i+2] = id[2];
		}
		for (int i=0; i<nelnb; ++i)
		{
			vector<int>& id = mesh.Node(elb.m_node[i]).m_ID;
			lm[3*(nelna+i)  ] = id[0];
			lm[3*(nelna+i)+1] = id[1];
			lm[3*(nelna+i)+2] = id[2];
		}
		ke.SetIndices(lm);

		// setup the en vector
		en.resize(nelna + nelnb);
		for (int i=0; i<nelna; i++) en[i        ] = ela.m_node[i];
		for (int i=0; i<nelnb; i++) en[i + nelna] = elb.m_node[i];
		ke.SetNodes(en);

		// assemble into global matrix
		LS.Assemble(ke);

		// don't forget to increment data counter
		nd += el.GaussPoints();
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::ElementStiffnessMatrixDG1(FESurfaceElement& face, FEInternalSurface2O::Data* pdata, matrix& ke)
{
	FESurface& surf = *m_surf.GetSurface();
	FEMesh& mesh = *GetMesh();

	// get the solid elements of this interface element
	FESolidElement& ela = static_cast<FESolidElement&>(*face.m_elem[1]);
	FESolidElement& elb = static_cast<FESolidElement&>(*face.m_elem[0]);

	// get the number nodes of each element
	int nelna = ela.Nodes();
	int nelnb = elb.Nodes();

	// initialize the element stiffness matrix
	int ndof = 3*(nelna + nelnb);

	// nodal coordinates
	vec3d Xa[FEElement::MAX_NODES];
	vec3d Xb[FEElement::MAX_NODES];
	for (int i=0; i<nelna; ++i) Xa[i] = mesh.Node(ela.m_node[i]).m_r0;
	for (int i=0; i<nelnb; ++i) Xb[i] = mesh.Node(elb.m_node[i]).m_r0;

	// shape function derivatives
	vec3d Ga[FEElement::MAX_NODES];
	vec3d Gb[FEElement::MAX_NODES];
	mat3d Ha[FEElement::MAX_NODES];
	mat3d Hb[FEElement::MAX_NODES];

	// loop over all integration points
	int nint = face.GaussPoints();
	double* gw = face.GaussWeights();
	for (int ng=0; ng<nint; ++ng)
	{
		// get the integration point data
		FEInternalSurface2O::Data& data = pdata[ng];

		// Jacobian and normal at this integration point
		vec3d nu;
		double J0 = surf.jac0(face, ng, nu);
		double Nu[3] = {nu.x, nu.y, nu.z};

		// shape function gradients at this integration point
		ShapeGradient0(ela, data.ksi[1].x, data.ksi[1].y, data.ksi[1].z, Ga);
		ShapeGradient0(elb, data.ksi[0].x, data.ksi[0].y, data.ksi[0].z, Gb);

		shape_gradient2(ela, Xa, data.ksi[1].x, data.ksi[1].y, data.ksi[1].z, Ha);
		shape_gradient2(elb, Xb, data.ksi[0].x, data.ksi[0].y, data.ksi[0].z, Hb);

		// The a-a term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelna; ++b)
			{
				vec3d& GA = Ga[a];
				vec3d& GB = Ga[b];
				mat3d& G2B = Ha[b];

				const tens5d& H = data.H[1];
				const tens6d& J = data.J[1];

				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
						{
							kik += GA.x*H(i,0,l,k,0)*GB.x + GA.x*H(i,0,l,k,1)*GB.y + GA.x*H(i,0,l,k,2)*GB.z;
							kik += GA.y*H(i,1,l,k,0)*GB.x + GA.y*H(i,1,l,k,1)*GB.y + GA.y*H(i,1,l,k,2)*GB.z;
							kik += GA.z*H(i,2,l,k,0)*GB.x + GA.z*H(i,2,l,k,1)*GB.y + GA.z*H(i,2,l,k,2)*GB.z;
						}
						kab[i][k] = kik*J0*gw[ng]*0.5;
					}

				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
							for (int m=0; m<3; ++m)
								for (int n=0; n<3; ++n)
								{
									kik += GA.x*J(i,0,k,l,m,n)*G2B[l][m] + GA.y*J(i,1,k,l,m,n)*G2B[l][m] + GA.z*J(i,2,k,l,m,n)*G2B[l][m];
								}
						kab[i][k] += kik*J0*gw[ng]*0.5;
					}

				// add it to the element matrix
				ke.add(3*a, 3*b, kab);
			}

		// The a-b term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelnb; ++b)
			{
				vec3d& GA = Ga[a];
				vec3d& GB = Gb[b];
				mat3d& G2B = Hb[b];

				const tens5d& H = data.H[0];
				const tens6d& J = data.J[0];

				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
						{
							kik += GA.x*H(i,0,l,k,0)*GB.x + GA.x*H(i,0,l,k,1)*GB.y + GA.x*H(i,0,l,k,2)*GB.z;
							kik += GA.y*H(i,1,l,k,0)*GB.x + GA.y*H(i,1,l,k,1)*GB.y + GA.y*H(i,1,l,k,2)*GB.z;
							kik += GA.z*H(i,2,l,k,0)*GB.x + GA.z*H(i,2,l,k,1)*GB.y + GA.z*H(i,2,l,k,2)*GB.z;
						}
						kab[i][k] = kik*J0*gw[ng]*0.5;
					}

				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
							for (int m=0; m<3; ++m)
								for (int n=0; n<3; ++n)
								{
									kik += GA.x*J(i,0,k,l,m,n)*G2B[l][m] + GA.y*J(i,1,k,l,m,n)*G2B[l][m] + GA.z*J(i,2,k,l,m,n)*G2B[l][m];
								}
						kab[i][k] += kik*J0*gw[ng]*0.5;
					}

				// add it to the element matrix
				ke.add(3*a, 3*(nelna + b), kab);
			}

		// The b-a term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelna; ++b)
			{
				vec3d& GA = Gb[a];
				vec3d& GB = Ga[b];
				mat3d& G2B = Ha[b];

				const tens5d& H = data.H[1];
				const tens6d& J = data.J[1];

				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
						{
							kik += GA.x*H(i,0,l,k,0)*GB.x + GA.x*H(i,0,l,k,1)*GB.y + GA.x*H(i,0,l,k,2)*GB.z;
							kik += GA.y*H(i,1,l,k,0)*GB.x + GA.y*H(i,1,l,k,1)*GB.y + GA.y*H(i,1,l,k,2)*GB.z;
							kik += GA.z*H(i,2,l,k,0)*GB.x + GA.z*H(i,2,l,k,1)*GB.y + GA.z*H(i,2,l,k,2)*GB.z;
						}
						kab[i][k] = kik*J0*gw[ng]*0.5;
					}

				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
							for (int m=0; m<3; ++m)
								for (int n=0; n<3; ++n)
								{
									kik += GA.x*J(i,0,k,l,m,n)*G2B[l][m] + GA.y*J(i,1,k,l,m,n)*G2B[l][m] + GA.z*J(i,2,k,l,m,n)*G2B[l][m];
								}
						kab[i][k] += kik*J0*gw[ng]*0.5;
					}

				// add it to the element matrix
				ke.sub(3*(nelna + a), 3*b, kab);
			}

		// The b-b term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelnb; ++b)
			{
				vec3d& GA = Gb[a];
				vec3d& GB = Gb[b];
				mat3d& G2B = Hb[b];

				const tens5d& H = data.H[0];
				const tens6d& J = data.J[0];

				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
						{
							kik += GA.x*H(i,0,l,k,0)*GB.x + GA.x*H(i,0,l,k,1)*GB.y + GA.x*H(i,0,l,k,2)*GB.z;
							kik += GA.y*H(i,1,l,k,0)*GB.x + GA.y*H(i,1,l,k,1)*GB.y + GA.y*H(i,1,l,k,2)*GB.z;
							kik += GA.z*H(i,2,l,k,0)*GB.x + GA.z*H(i,2,l,k,1)*GB.y + GA.z*H(i,2,l,k,2)*GB.z;
						}
						kab[i][k] = kik*J0*gw[ng]*0.5;
					}

				for (int i=0; i<3; ++i)
					for (int k=0; k<3; ++k)
					{
						double kik = 0.0;
						for (int l=0; l<3; ++l)
							for (int m=0; m<3; ++m)
								for (int n=0; n<3; ++n)
								{
									kik += GA.x*J(i,0,k,l,m,n)*G2B[l][m] + GA.y*J(i,1,k,l,m,n)*G2B[l][m] + GA.z*J(i,2,k,l,m,n)*G2B[l][m];
								}
						kab[i][k] += kik*J0*gw[ng]*0.5;
					}

				// add it to the element matrix
				ke.sub(3*(nelna + a), 3*(nelna + b), kab);
			}
	}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::ElementStiffnessMatrixDG2(FESurfaceElement& face, FEInternalSurface2O::Data* pdata, matrix& ke)
{
	// get the mesh and surface
	FEMesh& mesh = *GetMesh();
	FESurface& surf = *m_surf.GetSurface();

	// get the solid elements of this interface element
	FESolidElement& ela = static_cast<FESolidElement&>(*face.m_elem[1]);
	FESolidElement& elb = static_cast<FESolidElement&>(*face.m_elem[0]);

	// get the number nodes of each element
	int nelna = ela.Nodes();
	int nelnb = elb.Nodes();

	// initialize the element stiffness matrix
	int ndof = 3*(nelna + nelnb);

	// shape function derivatives
	vec3d G1a[FEElement::MAX_NODES];
	vec3d G1b[FEElement::MAX_NODES];
	mat3d G2a[FEElement::MAX_NODES];
	mat3d G2b[FEElement::MAX_NODES];

	// nodal coordinates
	vec3d Xa[FEElement::MAX_NODES];
	vec3d Xb[FEElement::MAX_NODES];
	for (int i=0; i<nelna; ++i) Xa[i] = mesh.Node(ela.m_node[i]).m_r0;
	for (int i=0; i<nelnb; ++i) Xb[i] = mesh.Node(elb.m_node[i]).m_r0;

	// loop over all integration points
	int nint = face.GaussPoints();
	double* gw = face.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		// get the integration point data
		FEInternalSurface2O::Data& data = pdata[n];

		// Jacobian and normal at this integration point
		vec3d nu;
		double J0 = surf.jac0(face, n, nu);
		double Nu[3] = {nu.x, nu.y, nu.z};

		// shape function gradients at this integration point
		ShapeGradient0(ela, data.ksi[1].x, data.ksi[1].y, data.ksi[1].z, G1a);
		ShapeGradient0(elb, data.ksi[0].x, data.ksi[0].y, data.ksi[0].z, G1b);
		shape_gradient2(ela, Xa, data.ksi[1].x, data.ksi[1].y, data.ksi[1].z, G2a);
		shape_gradient2(elb, Xb, data.ksi[0].x, data.ksi[0].y, data.ksi[0].z, G2b);

		// stiffnesses
		const tens5d& Ha = data.H0[1];
		const tens5d& Hb = data.H0[0];

		const tens6d& Ja = data.J0[1];
		const tens6d& Jb = data.J0[0];

		// The a-a term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelna; ++b)
			{
				double Ax[3][3][3] = {0};
				double Ay[3][3][3] = {0};
				double Az[3][3][3] = {0};
				for (int i=0; i<3; ++i)
					for (int j=0; j<3; ++j)
						for (int k=0; k<3; ++k)
						{
							for (int q=0; q<3; ++q)
								for (int r=0; r<3; ++r)
								{
									Ax[i][j][k] += Ja(i, j, k, 0, q, r)*G2a[a](q ,r);
									Ay[i][j][k] += Ja(i, j, k, 1, q, r)*G2a[a](q ,r);
									Az[i][j][k] += Ja(i, j, k, 2, q, r)*G2a[a](q ,r);
								}

							Ax[i][j][k] += Ha(i, j, k, 0, 0)*G1a[a].x + Ha(i, j, k, 0, 1)*G1a[a].y + Ha(i, j, k, 0, 2)*G1a[a].z;
							Ay[i][j][k] += Ha(i, j, k, 1, 0)*G1a[a].x + Ha(i, j, k, 1, 1)*G1a[a].y + Ha(i, j, k, 1, 2)*G1a[a].z;
							Az[i][j][k] += Ha(i, j, k, 2, 0)*G1a[a].x + Ha(i, j, k, 2, 1)*G1a[a].y + Ha(i, j, k, 2, 2)*G1a[a].z;
						}

				mat3d kab; kab.zero();
				for (int k=0; k<3; ++k)
				{
					kab(0,0) += G1a[b].x*Ax[0][0][k]*Nu[k] + G1a[b].y*Ax[0][1][k]*Nu[k] + G1a[b].z*Ax[0][2][k]*Nu[k];
					kab(0,1) += G1a[b].x*Ax[1][0][k]*Nu[k] + G1a[b].y*Ax[1][1][k]*Nu[k] + G1a[b].z*Ax[1][2][k]*Nu[k];
					kab(0,2) += G1a[b].x*Ax[2][0][k]*Nu[k] + G1a[b].y*Ax[2][1][k]*Nu[k] + G1a[b].z*Ax[2][2][k]*Nu[k];

					kab(1,0) += G1a[b].x*Ay[0][0][k]*Nu[k] + G1a[b].y*Ay[0][1][k]*Nu[k] + G1a[b].z*Ay[0][2][k]*Nu[k];
					kab(1,1) += G1a[b].x*Ay[1][0][k]*Nu[k] + G1a[b].y*Ay[1][1][k]*Nu[k] + G1a[b].z*Ay[1][2][k]*Nu[k];
					kab(1,2) += G1a[b].x*Ay[2][0][k]*Nu[k] + G1a[b].y*Ay[2][1][k]*Nu[k] + G1a[b].z*Ay[2][2][k]*Nu[k];

					kab(2,0) += G1a[b].x*Az[0][0][k]*Nu[k] + G1a[b].y*Az[0][1][k]*Nu[k] + G1a[b].z*Az[0][2][k]*Nu[k];
					kab(2,1) += G1a[b].x*Az[1][0][k]*Nu[k] + G1a[b].y*Az[1][1][k]*Nu[k] + G1a[b].z*Az[1][2][k]*Nu[k];
					kab(2,2) += G1a[b].x*Az[2][0][k]*Nu[k] + G1a[b].y*Az[2][1][k]*Nu[k] + G1a[b].z*Az[2][2][k]*Nu[k];
				}
				// multiply by weights
				kab *= J0*gw[n]*0.5;

				ke.add(3*a, 3*b, kab);
			}

		// The a-b term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelnb; ++b)
			{
				double Ax[3][3][3] = {0};
				double Ay[3][3][3] = {0};
				double Az[3][3][3] = {0};
				for (int i=0; i<3; ++i)
					for (int j=0; j<3; ++j)
						for (int k=0; k<3; ++k)
						{
							for (int q=0; q<3; ++q)
								for (int r=0; r<3; ++r)
								{
									Ax[i][j][k] += Ja(i, j, k, 0, q, r)*G2a[a](q ,r);
									Ay[i][j][k] += Ja(i, j, k, 1, q, r)*G2a[a](q ,r);
									Az[i][j][k] += Ja(i, j, k, 2, q, r)*G2a[a](q ,r);
								}

							Ax[i][j][k] += Ha(i, j, k, 0, 0)*G1a[a].x + Ha(i, j, k, 0, 1)*G1a[a].y + Ha(i, j, k, 0, 2)*G1a[a].z;
							Ay[i][j][k] += Ha(i, j, k, 1, 0)*G1a[a].x + Ha(i, j, k, 1, 1)*G1a[a].y + Ha(i, j, k, 1, 2)*G1a[a].z;
							Az[i][j][k] += Ha(i, j, k, 2, 0)*G1a[a].x + Ha(i, j, k, 2, 1)*G1a[a].y + Ha(i, j, k, 2, 2)*G1a[a].z;
						}

				mat3d kab; kab.zero();
				for (int k=0; k<3; ++k)
				{
					kab(0,0) += G1b[b].x*Ax[0][0][k]*Nu[k] + G1b[b].y*Ax[0][1][k]*Nu[k] + G1b[b].z*Ax[0][2][k]*Nu[k];
					kab(0,1) += G1b[b].x*Ax[1][0][k]*Nu[k] + G1b[b].y*Ax[1][1][k]*Nu[k] + G1b[b].z*Ax[1][2][k]*Nu[k];
					kab(0,2) += G1b[b].x*Ax[2][0][k]*Nu[k] + G1b[b].y*Ax[2][1][k]*Nu[k] + G1b[b].z*Ax[2][2][k]*Nu[k];

					kab(1,0) += G1b[b].x*Ay[0][0][k]*Nu[k] + G1b[b].y*Ay[0][1][k]*Nu[k] + G1b[b].z*Ay[0][2][k]*Nu[k];
					kab(1,1) += G1b[b].x*Ay[1][0][k]*Nu[k] + G1b[b].y*Ay[1][1][k]*Nu[k] + G1b[b].z*Ay[1][2][k]*Nu[k];
					kab(1,2) += G1b[b].x*Ay[2][0][k]*Nu[k] + G1b[b].y*Ay[2][1][k]*Nu[k] + G1b[b].z*Ay[2][2][k]*Nu[k];

					kab(2,0) += G1b[b].x*Az[0][0][k]*Nu[k] + G1b[b].y*Az[0][1][k]*Nu[k] + G1b[b].z*Az[0][2][k]*Nu[k];
					kab(2,1) += G1b[b].x*Az[1][0][k]*Nu[k] + G1b[b].y*Az[1][1][k]*Nu[k] + G1b[b].z*Az[1][2][k]*Nu[k];
					kab(2,2) += G1b[b].x*Az[2][0][k]*Nu[k] + G1b[b].y*Az[2][1][k]*Nu[k] + G1b[b].z*Az[2][2][k]*Nu[k];
				}
				// multiply by weights
				kab *= J0*gw[n]*0.5;

				ke.sub(3*a, 3*(nelna + b), kab);
			}

		// The b-a term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelna; ++b)
			{
				double Ax[3][3][3] = {0};
				double Ay[3][3][3] = {0};
				double Az[3][3][3] = {0};
				for (int i=0; i<3; ++i)
					for (int j=0; j<3; ++j)
						for (int k=0; k<3; ++k)
						{
							for (int q=0; q<3; ++q)
								for (int r=0; r<3; ++r)
								{
									Ax[i][j][k] += Jb(i, j, k, 0, q, r)*G2b[a](q ,r);
									Ay[i][j][k] += Jb(i, j, k, 1, q, r)*G2b[a](q ,r);
									Az[i][j][k] += Jb(i, j, k, 2, q, r)*G2b[a](q ,r);
								}

							Ax[i][j][k] += Hb(i, j, k, 0, 0)*G1b[a].x + Hb(i, j, k, 0, 1)*G1b[a].y + Hb(i, j, k, 0, 2)*G1b[a].z;
							Ay[i][j][k] += Hb(i, j, k, 1, 0)*G1b[a].x + Hb(i, j, k, 1, 1)*G1b[a].y + Hb(i, j, k, 1, 2)*G1b[a].z;
							Az[i][j][k] += Hb(i, j, k, 2, 0)*G1b[a].x + Hb(i, j, k, 2, 1)*G1b[a].y + Hb(i, j, k, 2, 2)*G1b[a].z;
						}

				mat3d kab; kab.zero();
				for (int k=0; k<3; ++k)
				{
					kab(0,0) += G1a[b].x*Ax[0][0][k]*Nu[k] + G1a[b].y*Ax[0][1][k]*Nu[k] + G1a[b].z*Ax[0][2][k]*Nu[k];
					kab(0,1) += G1a[b].x*Ax[1][0][k]*Nu[k] + G1a[b].y*Ax[1][1][k]*Nu[k] + G1a[b].z*Ax[1][2][k]*Nu[k];
					kab(0,2) += G1a[b].x*Ax[2][0][k]*Nu[k] + G1a[b].y*Ax[2][1][k]*Nu[k] + G1a[b].z*Ax[2][2][k]*Nu[k];

					kab(1,0) += G1a[b].x*Ay[0][0][k]*Nu[k] + G1a[b].y*Ay[0][1][k]*Nu[k] + G1a[b].z*Ay[0][2][k]*Nu[k];
					kab(1,1) += G1a[b].x*Ay[1][0][k]*Nu[k] + G1a[b].y*Ay[1][1][k]*Nu[k] + G1a[b].z*Ay[1][2][k]*Nu[k];
					kab(1,2) += G1a[b].x*Ay[2][0][k]*Nu[k] + G1a[b].y*Ay[2][1][k]*Nu[k] + G1a[b].z*Ay[2][2][k]*Nu[k];

					kab(2,0) += G1a[b].x*Az[0][0][k]*Nu[k] + G1a[b].y*Az[0][1][k]*Nu[k] + G1a[b].z*Az[0][2][k]*Nu[k];
					kab(2,1) += G1a[b].x*Az[1][0][k]*Nu[k] + G1a[b].y*Az[1][1][k]*Nu[k] + G1a[b].z*Az[1][2][k]*Nu[k];
					kab(2,2) += G1a[b].x*Az[2][0][k]*Nu[k] + G1a[b].y*Az[2][1][k]*Nu[k] + G1a[b].z*Az[2][2][k]*Nu[k];
				}
				// multiply by weights
				kab *= J0*gw[n]*0.5;

				ke.add(3*(nelna + a), 3*b, kab);
			}

		// The b-b term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelnb; ++b)
			{
				double Ax[3][3][3] = {0};
				double Ay[3][3][3] = {0};
				double Az[3][3][3] = {0};
				for (int i=0; i<3; ++i)
					for (int j=0; j<3; ++j)
						for (int k=0; k<3; ++k)
						{
							for (int q=0; q<3; ++q)
								for (int r=0; r<3; ++r)
								{
									Ax[i][j][k] += Jb(i, j, k, 0, q, r)*G2b[a](q ,r);
									Ay[i][j][k] += Jb(i, j, k, 1, q, r)*G2b[a](q ,r);
									Az[i][j][k] += Jb(i, j, k, 2, q, r)*G2b[a](q ,r);
								}

							Ax[i][j][k] += Hb(i, j, k, 0, 0)*G1b[a].x + Hb(i, j, k, 0, 1)*G1b[a].y + Hb(i, j, k, 0, 2)*G1b[a].z;
							Ay[i][j][k] += Hb(i, j, k, 1, 0)*G1b[a].x + Hb(i, j, k, 1, 1)*G1b[a].y + Hb(i, j, k, 1, 2)*G1b[a].z;
							Az[i][j][k] += Hb(i, j, k, 2, 0)*G1b[a].x + Hb(i, j, k, 2, 1)*G1b[a].y + Hb(i, j, k, 2, 2)*G1b[a].z;
						}

				mat3d kab; kab.zero();
				for (int k=0; k<3; ++k)
				{
					kab(0,0) += G1b[b].x*Ax[0][0][k]*Nu[k] + G1b[b].y*Ax[0][1][k]*Nu[k] + G1b[b].z*Ax[0][2][k]*Nu[k];
					kab(0,1) += G1b[b].x*Ax[1][0][k]*Nu[k] + G1b[b].y*Ax[1][1][k]*Nu[k] + G1b[b].z*Ax[1][2][k]*Nu[k];
					kab(0,2) += G1b[b].x*Ax[2][0][k]*Nu[k] + G1b[b].y*Ax[2][1][k]*Nu[k] + G1b[b].z*Ax[2][2][k]*Nu[k];

					kab(1,0) += G1b[b].x*Ay[0][0][k]*Nu[k] + G1b[b].y*Ay[0][1][k]*Nu[k] + G1b[b].z*Ay[0][2][k]*Nu[k];
					kab(1,1) += G1b[b].x*Ay[1][0][k]*Nu[k] + G1b[b].y*Ay[1][1][k]*Nu[k] + G1b[b].z*Ay[1][2][k]*Nu[k];
					kab(1,2) += G1b[b].x*Ay[2][0][k]*Nu[k] + G1b[b].y*Ay[2][1][k]*Nu[k] + G1b[b].z*Ay[2][2][k]*Nu[k];

					kab(2,0) += G1b[b].x*Az[0][0][k]*Nu[k] + G1b[b].y*Az[0][1][k]*Nu[k] + G1b[b].z*Az[0][2][k]*Nu[k];
					kab(2,1) += G1b[b].x*Az[1][0][k]*Nu[k] + G1b[b].y*Az[1][1][k]*Nu[k] + G1b[b].z*Az[1][2][k]*Nu[k];
					kab(2,2) += G1b[b].x*Az[2][0][k]*Nu[k] + G1b[b].y*Az[2][1][k]*Nu[k] + G1b[b].z*Az[2][2][k]*Nu[k];
				}
				// multiply by weights
				kab *= J0*gw[n]*0.5;

				ke.sub(3*(nelna + a), 3*(nelna + b), kab);
			}
		}
}

//-----------------------------------------------------------------------------
void FEElasticSolidDomain2O::ElementStiffnessMatrixDG3(FESurfaceElement& face, FEInternalSurface2O::Data* pdata, matrix& ke)
{
	// get the beta parameter
	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(GetMaterial());
	double beta = pmat->m_beta;

	double h = m_surf.GetElementSize();

	// get the surface
	FESurface& surf = *m_surf.GetSurface();

	// get the solid elements of this interface element
	FESolidElement& ela = static_cast<FESolidElement&>(*face.m_elem[1]);
	FESolidElement& elb = static_cast<FESolidElement&>(*face.m_elem[0]);

	// get the number nodes of each element
	int nelna = ela.Nodes();
	int nelnb = elb.Nodes();

	// initialize the element stiffness matrix
	int ndof = 3*(nelna + nelnb);

	// shape function derivatives
	vec3d Ga[FEElement::MAX_NODES];
	vec3d Gb[FEElement::MAX_NODES];

	// loop over all integration points
	int nint = face.GaussPoints();
	double* gw = face.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		// get the integration point data
		FEInternalSurface2O::Data& data = pdata[n];

		// Jacobian and normal at this integration point
		vec3d nu;
		double J0 = surf.jac0(face, n, nu);
		double Nu[3] = {nu.x, nu.y, nu.z};

		// shape function gradients at this integration point
		ShapeGradient0(ela, data.ksi[1].x, data.ksi[1].y, data.ksi[1].z, Ga);
		ShapeGradient0(elb, data.ksi[0].x, data.ksi[0].y, data.ksi[0].z, Gb);

		// average stiffness 
		tens6d& Javg = data.J0avg;

		// The ++ term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelna; b++)
			{
				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int p=0; p<3; ++p)
					{
						double kip = 0.0;
						for (int k=0; k<3; ++k)
							for (int r=0; r<3; ++r)
							{
								double NN = Nu[k]*Nu[r];
								kip += Ga[a].x*Ga[b].x*Javg(i,0,k,p,0,r)*NN + Ga[a].x*Ga[b].y*Javg(i,0,k,p,1,r)*NN + Ga[a].x*Ga[b].z*Javg(i,0,k,p,2,r)*NN;
								kip += Ga[a].y*Ga[b].x*Javg(i,1,k,p,0,r)*NN + Ga[a].y*Ga[b].y*Javg(i,1,k,p,1,r)*NN + Ga[a].y*Ga[b].z*Javg(i,1,k,p,2,r)*NN;
								kip += Ga[a].z*Ga[b].x*Javg(i,2,k,p,0,r)*NN + Ga[a].z*Ga[b].y*Javg(i,2,k,p,1,r)*NN + Ga[a].z*Ga[b].z*Javg(i,2,k,p,2,r)*NN;
							}

						kab[i][p] = kip*(beta/h)*J0*gw[n];
					}

				ke.add(3*a, 3*b, kab);
			}

		// The +- term
		for (int a=0; a<nelna; ++a)
			for (int b=0; b<nelnb; b++)
			{
				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int p=0; p<3; ++p)
					{
						double kip = 0.0;
						for (int k=0; k<3; ++k)
							for (int r=0; r<3; ++r)
							{
								double NN = Nu[k]*Nu[r];
								kip += Ga[a].x*Gb[b].x*Javg(i,0,k,p,0,r)*NN + Ga[a].x*Gb[b].y*Javg(i,0,k,p,1,r)*NN + Ga[a].x*Gb[b].z*Javg(i,0,k,p,2,r)*NN;
								kip += Ga[a].y*Gb[b].x*Javg(i,1,k,p,0,r)*NN + Ga[a].y*Gb[b].y*Javg(i,1,k,p,1,r)*NN + Ga[a].y*Gb[b].z*Javg(i,1,k,p,2,r)*NN;
								kip += Ga[a].z*Gb[b].x*Javg(i,2,k,p,0,r)*NN + Ga[a].z*Gb[b].y*Javg(i,2,k,p,1,r)*NN + Ga[a].z*Gb[b].z*Javg(i,2,k,p,2,r)*NN;
							}

						kab[i][p] = kip*(beta/h)*J0*gw[n];
					}

				ke.sub(3*a, 3*(nelna + b), kab);
			}

		// The -+ term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelna; b++)
			{
				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int p=0; p<3; ++p)
					{
						double kip = 0.0;
						for (int k=0; k<3; ++k)
							for (int r=0; r<3; ++r)
							{
								double NN = Nu[k]*Nu[r];
								kip += Gb[a].x*Ga[b].x*Javg(i,0,k,p,0,r)*NN + Gb[a].x*Ga[b].y*Javg(i,0,k,p,1,r)*NN + Gb[a].x*Ga[b].z*Javg(i,0,k,p,2,r)*NN;
								kip += Gb[a].y*Ga[b].x*Javg(i,1,k,p,0,r)*NN + Gb[a].y*Ga[b].y*Javg(i,1,k,p,1,r)*NN + Gb[a].y*Ga[b].z*Javg(i,1,k,p,2,r)*NN;
								kip += Gb[a].z*Ga[b].x*Javg(i,2,k,p,0,r)*NN + Gb[a].z*Ga[b].y*Javg(i,2,k,p,1,r)*NN + Gb[a].z*Ga[b].z*Javg(i,2,k,p,2,r)*NN;
							}

						kab[i][p] = kip*(beta/h)*J0*gw[n];
					}

				ke.sub(3*(a+nelna), 3*b, kab);
			}

		// The -- term
		for (int a=0; a<nelnb; ++a)
			for (int b=0; b<nelnb; b++)
			{
				mat3d kab;
				for (int i=0; i<3; ++i)
					for (int p=0; p<3; ++p)
					{
						double kip = 0.0;
						for (int k=0; k<3; ++k)
							for (int r=0; r<3; ++r)
							{
								double NN = Nu[k]*Nu[r];
								kip += Gb[a].x*Gb[b].x*Javg(i,0,k,p,0,r)*NN + Gb[a].x*Gb[b].y*Javg(i,0,k,p,1,r)*NN + Gb[a].x*Gb[b].z*Javg(i,0,k,p,2,r)*NN;
								kip += Gb[a].y*Gb[b].x*Javg(i,1,k,p,0,r)*NN + Gb[a].y*Gb[b].y*Javg(i,1,k,p,1,r)*NN + Gb[a].y*Gb[b].z*Javg(i,1,k,p,2,r)*NN;
								kip += Gb[a].z*Gb[b].x*Javg(i,2,k,p,0,r)*NN + Gb[a].z*Gb[b].y*Javg(i,2,k,p,1,r)*NN + Gb[a].z*Gb[b].z*Javg(i,2,k,p,2,r)*NN;
							}

						kab[i][p] = kip*(beta/h)*J0*gw[n];
					}

				ke.add(3*(nelna + a), 3*(nelna + b), kab);
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix
void FEElasticSolidDomain2O::ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke)
{
	FEElasticMaterial2O* pmat = dynamic_cast<FEElasticMaterial2O*>(m_pMat);

	FESolidElement& el = Element(iel);

	// Get the current element's data
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	// Gx = dN/dx
	vec3d G[FEElement::MAX_NODES];

	// H = d2N/dXdX
	mat3d G2[FEElement::MAX_NODES];

	// get the initial nodal coordinates
	vec3d X[FEElement::MAX_NODES];
	FEMesh& mesh = *GetMesh();
	mesh.GetInitialNodalCoordinates(el, X);

	ke.zero();
	
	// calculate element stiffness matrix
	const int nint = el.GaussPoints();
	const double *gw = el.GaussWeights();
	for (int ni=0; ni<nint; ++ni)
	{
		// get the material point data
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(ni);
		FEElasticMaterialPoint2O& pt2O = *(mp.ExtractData<FEElasticMaterialPoint2O>());

		// get the material tangents
		tens4d C;
		tens5d L;
		tens5d H;
		tens6d J;
		pmat->Tangent(mp, C, L, H, J);

		// get the first derivative of shape functions
		double J0 = ShapeGradient0(el, ni, G);

		// second derivative of shape functions
		shape_gradient2(el, X, ni, G2);

		for (int a=0; a<neln; ++a)
		{
			double Ga[3] = {G[a].x, G[a].y, G[a].z};
			mat3d& Ha = G2[a];

			for (int b=0; b<neln; ++b)
			{
				double Gb[3] = {G[b].x, G[b].y, G[b].z};
				mat3d& Hb = G2[b];

				// gradN*C*gradN
				mat3d Kab; Kab.zero();
				for (int j=0; j<3; ++j)
					for (int l=0; l<3; ++l)
					{
						Kab[0][0] += (Ga[j]*C(0,j,0,l)*Gb[l]);
						Kab[0][1] += (Ga[j]*C(0,j,1,l)*Gb[l]);
						Kab[0][2] += (Ga[j]*C(0,j,2,l)*Gb[l]);

						Kab[1][0] += (Ga[j]*C(1,j,0,l)*Gb[l]);
						Kab[1][1] += (Ga[j]*C(1,j,1,l)*Gb[l]);
						Kab[1][2] += (Ga[j]*C(1,j,2,l)*Gb[l]);

						Kab[2][0] += (Ga[j]*C(2,j,0,l)*Gb[l]);
						Kab[2][1] += (Ga[j]*C(2,j,1,l)*Gb[l]);
						Kab[2][2] += (Ga[j]*C(2,j,2,l)*Gb[l]);

					}

				// gradN*L*grad2N
				for (int j=0; j<3; ++j)
					for (int l=0; l<3; ++l)
						for (int m=0; m<3; ++m)
						{
							Kab[0][0] += (Ga[j]*L(0, j, 0, l, m)*Ha(l, m));
							Kab[0][1] += (Ga[j]*L(0, j, 1, l, m)*Ha(l, m));
							Kab[0][2] += (Ga[j]*L(0, j, 2, l, m)*Ha(l, m));

							Kab[1][0] += (Ga[j]*L(1, j, 0, l, m)*Ha(l, m));
							Kab[1][1] += (Ga[j]*L(1, j, 1, l, m)*Ha(l, m));
							Kab[1][2] += (Ga[j]*L(1, j, 2, l, m)*Ha(l, m));

							Kab[2][0] += (Ga[j]*L(2, j, 0, l, m)*Ha(l, m));
							Kab[2][1] += (Ga[j]*L(2, j, 1, l, m)*Ha(l, m));
							Kab[2][2] += (Ga[j]*L(2, j, 2, l, m)*Ha(l, m));
						}

				// grad2N*H*gradN
				for (int j=0; j<3; ++j)
					for (int k=0; k<3; ++k)
						for (int m=0; m<3; ++m)
						{
							Kab[0][0] += (Ha(j,k)*H(0, j, k, 0, m)*Ga[m]);
							Kab[0][1] += (Ha(j,k)*H(0, j, k, 1, m)*Ga[m]);
							Kab[0][2] += (Ha(j,k)*H(0, j, k, 2, m)*Ga[m]);

							Kab[1][0] += (Ha(j,k)*H(1, j, k, 0, m)*Ga[m]);
							Kab[1][1] += (Ha(j,k)*H(1, j, k, 1, m)*Ga[m]);
							Kab[1][2] += (Ha(j,k)*H(1, j, k, 2, m)*Ga[m]);

							Kab[2][0] += (Ha(j,k)*H(2, j, k, 0, m)*Ga[m]);
							Kab[2][1] += (Ha(j,k)*H(2, j, k, 1, m)*Ga[m]);
							Kab[2][2] += (Ha(j,k)*H(2, j, k, 2, m)*Ga[m]);
						}

				// grad2N*J*grad2N
				for (int j=0; j<3; ++j)
					for (int k=0; k<3; ++k)
						for (int m=0; m<3; ++m)
							for (int n=0; n<3; ++n)
							{
								Kab[0][0] += (Ha(j,k)*J(0, j, k, 0, m, n)*Hb(m,n));
								Kab[0][1] += (Ha(j,k)*J(0, j, k, 1, m, n)*Hb(m,n));
								Kab[0][2] += (Ha(j,k)*J(0, j, k, 2, m, n)*Hb(m,n));

								Kab[1][0] += (Ha(j,k)*J(1, j, k, 0, m, n)*Hb(m,n));
								Kab[1][1] += (Ha(j,k)*J(1, j, k, 1, m, n)*Hb(m,n));
								Kab[1][2] += (Ha(j,k)*J(1, j, k, 2, m, n)*Hb(m,n));

								Kab[2][0] += (Ha(j,k)*J(2, j, k, 0, m, n)*Hb(m,n));
								Kab[2][1] += (Ha(j,k)*J(2, j, k, 1, m, n)*Hb(m,n));
								Kab[2][2] += (Ha(j,k)*J(2, j, k, 2, m, n)*Hb(m,n));
							}

				// multiply by jacobian and integration weight
				Kab *= J0*gw[ni];

				// add it to the element matrix
				ke.add(3*a, 3*b, Kab);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the gradient of deformation gradient of element el at integration point n.
//! The gradient of the deformation gradient is returned in G.
void FEElasticSolidDomain2O::defhess(FESolidElement &el, int n, tens3drs &G)
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

	mat3d H[FEElement::MAX_NODES];
	shape_gradient2(el, X, n, H);

	// loop over nodes
	G.zero();
	for (int a=0; a<neln; ++a)
	{
		const mat3d& Ha = H[a];

		// calculate gradient of deformation gradient
		// Note that k >= j. Since tensdrs has symmetries this
		// prevents overwriting of symmetric components
		for (int j=0; j<3; ++j)
			for (int k=j; k<3; ++k)
			{
				G(0,j,k) += Ha(j,k)*x[a].x;
				G(1,j,k) += Ha(j,k)*x[a].y;
				G(2,j,k) += Ha(j,k)*x[a].z;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the gradient of deformation gradient of element el at integration point n.
//! The gradient of the deformation gradient is returned in G.
void FEElasticSolidDomain2O::defhess(FESolidElement &el, double r, double s, double t, tens3drs &G)
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

	mat3d H[FEElement::MAX_NODES];
	shape_gradient2(el, X, r, s, t, H);

	// loop over nodes
	G.zero();
	for (int a=0; a<neln; ++a)
	{
		mat3d& Ha = H[a];

		// calculate gradient of deformation gradient
		// Note that k >= j. Since tensdrs has symmetries this
		// prevents overwriting of symmetric components
		for (int j=0; j<3; ++j)
			for (int k=j; k<3; ++k)
			{
				G(0,j,k) += Ha[j][k]*x[a].x;
				G(1,j,k) += Ha[j][k]*x[a].y;
				G(2,j,k) += Ha[j][k]*x[a].z;
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the second derivative of shape function N[node] with respect
//! to the material coordinates.
void FEElasticSolidDomain2O::shape_gradient2(const FESolidElement& el, vec3d* X, int n, mat3d* H)
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
	for (int a=0; a<neln; ++a)
	{
		double G1[3];
		G1[0] = Gr[a];
		G1[1] = Gs[a];
		G1[2] = Gt[a];

		// second derivatives of shape functions
		double G2[3][3];
		G2[0][0] = Grrn[a]; G2[0][1] = Grsn[a]; G2[0][2] = Grtn[a];
		G2[1][0] = Gsrn[a]; G2[1][1] = Gssn[a]; G2[1][2] = Gstn[a];
		G2[2][0] = Gtrn[a]; G2[2][1] = Gtsn[a]; G2[2][2] = Gttn[a];

		// calculate dB/dr
		double D[3][3] = {0};
		for (int i=0; i<3; ++i)
			for (int k=0; k<3; ++k)
			{
				for (int j=0; j<3; ++j) D[i][k] += A[i][j][k]*G1[j] + Ji[j][i]*G2[j][k];
			}

		// calculate global gradient of shape functions
		mat3d& Ha = H[a];
		for (int i=0; i<3; ++i)
			for (int j=0; j<3; ++j)
			{
				Ha[i][j] = D[i][0]*Ji[0][j] + D[i][1]*Ji[1][j] + D[i][2]*Ji[2][j];
			}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the second derivative of shape function N[node] with respect
//! to the material coordinates.
void FEElasticSolidDomain2O::shape_gradient2(const FESolidElement& el, vec3d* X, double r, double s, double t, mat3d* H)
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
	double G2[3][3];
	for (int a=0; a<neln; ++a)
	{
		G1[0] = Gr[a];
		G1[1] = Gs[a];
		G1[2] = Gt[a];

		// second derivatives of shape functions
		G2[0][0] = Grr[a]; G2[0][1] = Grs[a]; G2[0][2] = Grt[a];
		G2[1][0] = Grs[a]; G2[1][1] = Gss[a]; G2[1][2] = Gst[a];
		G2[2][0] = Grt[a]; G2[2][1] = Gst[a]; G2[2][2] = Gtt[a];

		// calculate dB/dr
		double D[3][3] = {0};
		for (int i=0; i<3; ++i)
			for (int k=0; k<3; ++k)
			{
				for (int j=0; j<3; ++j) D[i][k] += A[i][j][k]*G1[j] + Ji[j][i]*G2[j][k];
			}

		// calculate global gradient of shape functions
		mat3d& Ha = H[a];
		for (int i=0; i<3; ++i)
			for (int j=0; j<3; ++j)
			{
				Ha[i][j] = D[i][0]*Ji[0][j] + D[i][1]*Ji[1][j] + D[i][2]*Ji[2][j];
			}
	}
}
