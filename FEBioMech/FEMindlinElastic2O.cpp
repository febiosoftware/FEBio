#include "stdafx.h"
#include "FEMindlinElastic2O.h"

BEGIN_PARAMETER_LIST(FEMindlinElastic2O, FEElasticMaterial2O)
	ADD_PARAMETER(m_lam, FE_PARAM_DOUBLE, "lam");
	ADD_PARAMETER(m_mu , FE_PARAM_DOUBLE, "mu" );
	ADD_PARAMETER(m_a1 , FE_PARAM_DOUBLE, "a1" );
	ADD_PARAMETER(m_a2 , FE_PARAM_DOUBLE, "a2" );
	ADD_PARAMETER(m_a3 , FE_PARAM_DOUBLE, "a3" );
	ADD_PARAMETER(m_a4 , FE_PARAM_DOUBLE, "a4" );
	ADD_PARAMETER(m_a5 , FE_PARAM_DOUBLE, "a5" );
END_PARAMETER_LIST();

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
				a += 0.5*m_a5*(3.0*G(q,p,r) + G(r,p,q));

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
							a += 0.5*m_a5*(3.0*I(j,l)*I(i,m)*I(k,n) + I(k,l)*I(i,m)*I(j,n));

							J(i,j,k,l,m,n) = a;
						}

}
