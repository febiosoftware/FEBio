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
#include "FEMindlinElastic2O.h"

BEGIN_FECORE_CLASS(FEMindlinElastic2O, FEElasticMaterial2O)
	ADD_PARAMETER(m_lam, "lam");
	ADD_PARAMETER(m_mu , "mu" );
	ADD_PARAMETER(m_a1 , "a1" );
	ADD_PARAMETER(m_a2 , "a2" );
	ADD_PARAMETER(m_a3 , "a3" );
	ADD_PARAMETER(m_a4 , "a4" );
	ADD_PARAMETER(m_a5 , "a5" );
END_FECORE_CLASS();

FEMindlinElastic2O::FEMindlinElastic2O(FEModel* pfem) : FEElasticMaterial2O(pfem)
{
	m_mu = 0.0;
	m_lam = 0.0;
	m_a1 = 0.0;
	m_a2 = 0.0;
	m_a3 = 0.0;
	m_a4 = 0.0;
	m_a5 = 0.0;
}

//! Calculate PK1 stress and higher order stress Q
void FEMindlinElastic2O::Stress(FEMaterialPoint& mp, mat3d& P, tens3drs& Q)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEElasticMaterialPoint2O& pt2O = *mp.ExtractData<FEElasticMaterialPoint2O>();

	// deformation gradient and transpose
	const mat3d& F = pt.m_F;
	mat3d Ft = F.transpose();

	// gradient of deformation gradient
	const tens3drs& G = pt2O.m_G;

	// identity tensor
	mat3dd I(1.0);

	// Lagrange strain
	mat3d E = (Ft*F - I)*0.5;
	double trE = E.trace();

	// PK1 stress
	P = F*(m_lam*trE) + (F*E)*(2.0*m_mu);

	// higher order stress
	for (int p=0; p<3; ++p)
		for (int q=0; q<3; ++q)
			for (int r=0; r<3; ++r)
			{
				double a = 0.0;
				for (int i=0; i<3; ++i)
				{
					a += m_a1*(I(p,r)*G(i,q,i) + I(p,q)*G(i,i,r));
					a += 0.5*m_a2*(2.0*I(q,r)*G(i,p,i) + G(q,i,i)*I(p,r) + G(r,i,i)*I(p,q));
					a += 2.0*m_a3*G(p,i,i)*I(q,r);
				}
				a += 2.0*m_a4*G(p,q,r);
				a += m_a5*(G(q,p,r) + G(r,p,q));

				Q(p,q,r) = a;
			}
}

//! Calculate material tangents
//! C = dP/dF
//! L = dP/dG
//! H = dQ/dF
//! J = dQ/dG
void FEMindlinElastic2O::Tangent(FEMaterialPoint& mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEElasticMaterialPoint2O& pt2O = *mp.ExtractData<FEElasticMaterialPoint2O>();

	// deformation gradient and transpose
	const mat3d& F = pt.m_F;
	mat3d Ft = F.transpose();

	// identity tensor
	mat3dd I(1.0);

	// Lagrange strain
	mat3d E = (Ft*F - I)*0.5;
	double trE = E.trace();

	mat3d B = F*Ft;

	// C-tensor
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
				{
					double c = 0.0;
					c += I(i,k)*(m_lam*trE*I(j,l) + 2.0*m_mu*E(l,j));
					c += m_lam*F(i,j)*F(k,l);
					c += m_mu*(F(i,l)*F(k,j) + B(k,i)*I(j,l));

					C(i,j,k,l) = c;
				}

	// L and H are zero
	L.zero();
	H.zero();

	// J-tensor
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
					for (int m=0; m<3; ++m)
						for (int n=0; n<3; ++n)
						{
							double a = 0.0;
							a += m_a1*(I(i,k)*I(l,n)*I(j,m) + I(i,j)*I(l,m)*I(k,n));
							a += 0.5*m_a2*(2.0*I(j,k)*I(l,n)*I(i,m) + I(j,l)*I(m,n)*I(i,k) + I(i,j)*I(k,l)*I(m,n));
							a += 2.0*m_a3*I(i,l)*I(m,n)*I(j,k);
							a += 2.0*m_a4*I(i,l)*I(j,m)*I(k,n);
							a += m_a5*(I(j,l)*I(i,m)*I(k,n) + I(k,l)*I(i,m)*I(j,n));

							J(i,j,k,l,m,n) = a;
						}

}
