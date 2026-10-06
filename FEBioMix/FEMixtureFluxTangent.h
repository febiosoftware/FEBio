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

#pragma once
#include <FECore/mat3d.h>
#include <FECore/tens4d.h>
#include <vector>

//-----------------------------------------------------------------------------
// Linearization of the fluid flux w and of the solute molar fluxes j^a of a
// mixture (biphasic-solute, triphasic, multiphasic) with respect to the solid
// displacement, for the displacement interpolated by a shape function whose
// spatial gradient is gradN.
//
//   w   = -ke . ( grad p + R T sum_a (kappa_a/d0_a) d_a . grad c_a )
//   ke  = [ k^-1 + (R T/phiw) sum_a (kappa_a c_a/d0_a) (I - d_a/d0_a) ]^-1
//   j_a = kappa_a d_a . ( -phiw grad c_a + (c_a/d0_a) w )
//
// where c_a are the effective concentrations, kappa_a(J) the partition
// coefficients, and phiw(J) = 1 - phisr/J the porosity.
//
// On output, column m of wu (and of ju[a]) is J^-1 F . D(J F^-1 w)[e_m N],
// i.e., the push-forward of the directional derivative of the referential flux
// along a displacement increment e_m N, so that the contribution of the flux
// to the stiffness matrix is (wu^T . grad N_a) for the test function N_a.
//
// The derivatives are evaluated one displacement component at a time, with the
// spatial velocity gradient L = e_m (x) gradN. The spatial rates of change of k
// and d_a are recovered from the material tangents dKdE and dDdE, which are the
// push-forwards of 2 dK_ref/dC (with K_ref = J F^-1 k F^-T), as returned by the
// FEHydraulicPermeability and FESoluteDiffusivity classes:
//   Dk = dKdE : sym(L) - k tr(L) + L k + k L^T
//
// This function was verified against a finite-difference approximation of the
// referential fluxes.
inline void MixtureFluxTangent(
	const vec3d& gradN,							// spatial gradient of shape function
	double J, double phiw, double RT,
	const mat3ds& K, const tens4dmm& dKdE,		// hydraulic permeability and its strain tangent
	const mat3ds& Ke,							// effective permeability ke
	const vec3d& gradp,							// gradient of effective fluid pressure
	const vec3d& w,								// fluid flux
	const std::vector<mat3ds>& D,				// solute diffusivities
	const std::vector<tens4dmm>& dDdE,			// diffusivity strain tangents
	const std::vector<double>& D0,				// free diffusivities
	const std::vector<double>& kappa,			// partition coefficients
	const std::vector<double>& dkdJ,			// derivatives of kappa w.r.t. J
	const std::vector<double>& c,				// effective concentrations
	const std::vector<vec3d>& gradc,			// gradients of effective concentrations
	const std::vector<vec3d>& j,				// solute fluxes
	mat3d& wu,									// output: fluid flux tangent
	std::vector<mat3d>& ju)						// output: solute flux tangents
{
	const int nsol = (int)D.size();
	const double phis = 1.0 - phiw;
	mat3dd I(1);

	mat3d Km(K), Kim(K.inverse()), Kem(Ke);

	// g = grad p + R T sum_a (kappa_a/d0_a) d_a . grad c_a
	vec3d g = gradp;
	for (int a = 0; a < nsol; ++a) g += (D[a] * gradc[a])*(RT*kappa[a] / D0[a]);

	wu.zero();
	ju.resize(nsol);
	for (int a = 0; a < nsol; ++a) ju[a].zero();

	std::vector<mat3d> Dd(nsol);
	for (int m = 0; m < 3; ++m)
	{
		vec3d em((m == 0 ? 1 : 0), (m == 1 ? 1 : 0), (m == 2 ? 1 : 0));
		mat3d L = em & gradN;
		mat3d LT = L.transpose();
		mat3ds Ls = L.sym();
		double trL = L.trace();

		// tens4dmm::dot expects the off-diagonal components in engineering (Voigt) form
		mat3ds Le(Ls.xx(), Ls.yy(), Ls.zz(), 2 * Ls.xy(), 2 * Ls.yz(), 2 * Ls.xz());

		// derivative of k^-1
		mat3d Dk = mat3d(dKdE.dot(Le)) - Km*trL + L*Km + Km*LT;
		mat3d DKinv = -(Kim*Dk*Kim);

		// derivative of g (at fixed nodal values of p and c)
		vec3d Dg = -(LT*gradp);
		for (int a = 0; a < nsol; ++a)
		{
			mat3d Da(D[a]);
			Dd[a] = mat3d(dDdE[a].dot(Le)) - Da*trL + L*Da + Da*LT;
			mat3d ImD(mat3ds(I) - D[a] / D0[a]);
			DKinv += ImD*(RT*c[a] / D0[a] / phiw*(J*dkdJ[a] - kappa[a] * phis / phiw)*trL)
				- Dd[a] * (RT*kappa[a] * c[a] / (D0[a] * D0[a] * phiw));
			Dg += ((D[a] * gradc[a])*(J*dkdJ[a] * trL) + (Dd[a] * gradc[a])*kappa[a]
				- (Da*(LT*gradc[a]))*kappa[a])*(RT / D0[a]);
		}

		// derivative of w
		mat3d DKe = -(Kem*DKinv*Kem);
		vec3d Dw = -(DKe*g) - Kem*Dg;
		vec3d Qw = Dw + w*trL - L*w;
		wu[0][m] = Qw.x; wu[1][m] = Qw.y; wu[2][m] = Qw.z;

		// derivatives of j_a
		for (int a = 0; a < nsol; ++a)
		{
			mat3d Da(D[a]);
			vec3d ga = -gradc[a] * phiw + w*(c[a] / D0[a]);
			vec3d Dga = -gradc[a] * (phis*trL) + (LT*gradc[a])*phiw + Dw*(c[a] / D0[a]);
			vec3d Dj = (Da*ga)*(J*dkdJ[a] * trL) + (Dd[a] * ga)*kappa[a] + (Da*Dga)*kappa[a];
			vec3d Qj = Dj + j[a] * trL - L*j[a];
			ju[a][0][m] = Qj.x; ju[a][1][m] = Qj.y; ju[a][2][m] = Qj.z;
		}
	}
}
