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
#include "FEElasticMaterialPoint.h"

//-----------------------------------------------------------------------------
FEElasticMaterialPoint::FEElasticMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
	m_F.unit();
	m_J = 1;
	m_J_star = 1;
	m_s.zero();
    m_v = m_a = m_gradJ = vec3d(0, 0, 0);
    m_buncoupled = false;
    m_Wt = m_Wp = 0;
    m_p = 0;
    m_Lk = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEElasticMaterialPoint::Copy()
{
	FEElasticMaterialPoint* pt = new FEElasticMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEElasticMaterialPoint::Init()
{
	m_F.unit();

	m_J = 1;
	m_J_star = 1;

	m_s.zero();

    m_v = m_a = m_gradJ = vec3d(0, 0, 0);
    m_L.zero();
    
    m_Wt = m_Wp = 0;

    m_p = 0;
    m_Lk = 0;
    
	// don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEElasticMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_F & m_J & m_J_star & m_s & m_v & m_a & m_gradJ & m_L & m_Wt & m_Wp & m_p & m_Lk;
	ar & m_buncoupled;
}

//-----------------------------------------------------------------------------
//! Calculates the right Cauchy-Green tensor at the current material point

mat3ds FEElasticMaterialPoint::RightCauchyGreen() const
{
	// get the right Cauchy-Green tensor
	// C = Ft*F
	const mat3d& F = m_F;
	mat3ds C;
	C.xx() = F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]; // = C[0][0]
	C.yy() = F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]; // = C[1][1]
	C.zz() = F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]; // = C[2][2]
	C.xy() = F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]; // = C[0][1]
	C.yz() = F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]; // = C[1][2]
	C.xz() = F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]; // = C[0][2]

	return C;
}

//-----------------------------------------------------------------------------
//! Calculates the left Cauchy-Green tensor at the current material point

mat3ds FEElasticMaterialPoint::LeftCauchyGreen() const
{
	// get the left Cauchy-Green tensor
	// b = F*Ft
	const mat3d& F = m_F;
	mat3ds b;
	b.xx() = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]; // = b[0][0]
	b.yy() = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]; // = b[1][1]
	b.zz() = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]; // = b[2][2]
	b.xy() = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]; // = b[0][1]
	b.yz() = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]; // = b[1][2]
	b.xz() = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]; // = b[0][2]

	return b;
}

//-----------------------------------------------------------------------------
//! Calculates the right Cauchy-Green tensor at the current material point

mat3ds FEElasticMaterialPoint::DevRightCauchyGreen() const
{
	double Jm23 = pow(m_J, -2.0/3.0);

	// get the deviatoric right Cauchy-Green tensor
	// C = Ft*F
	const mat3d& F = m_F;
	mat3ds C;
	C.xx() = Jm23*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]); // = C[0][0]
	C.yy() = Jm23*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]); // = C[1][1]
	C.zz() = Jm23*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]); // = C[2][2]
	C.xy() = Jm23*(F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]); // = C[0][1]
	C.yz() = Jm23*(F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]); // = C[1][2]
	C.xz() = Jm23*(F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]); // = C[0][2]

	return C;
}

//-----------------------------------------------------------------------------
//! Calculates the left Cauchy-Green tensor at the current material point

mat3ds FEElasticMaterialPoint::DevLeftCauchyGreen() const
{
	double Jm23 = pow(m_J, -2.0/3.0);

	// get the left Cauchy-Green tensor
	// b = F*Ft
	const mat3d& F = m_F;
	mat3ds b;
	b.xx() = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]); // = b[0][0]
	b.yy() = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]); // = b[1][1]
	b.zz() = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]); // = b[2][2]
	b.xy() = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]); // = b[0][1]
	b.yz() = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]); // = b[1][2]
	b.xz() = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]); // = b[0][2]

	return b;
}

//-----------------------------------------------------------------------------
//! Calculates the right stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::RightStretch() const
{
    // get the right stretch tensor
    mat3ds C = RightCauchyGreen();
    double l2[3];
    vec3d v[3];
    C.eigen2(l2, v);
    mat3ds U = dyad(v[0])*sqrt(l2[0]) + dyad(v[1])*sqrt(l2[1]) + dyad(v[2])*sqrt(l2[2]);
    
    return U;
}

//-----------------------------------------------------------------------------
//! Calculates the left stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::LeftStretch() const
{
    // get the left stretch tensor
    mat3ds B = LeftCauchyGreen();
    double l2[3];
    vec3d v[3];
    B.eigen2(l2, v);
    mat3ds V = dyad(v[0])*sqrt(l2[0]) + dyad(v[1])*sqrt(l2[1]) + dyad(v[2])*sqrt(l2[2]);
    
    return V;
}

//-----------------------------------------------------------------------------
//! Calculates the right stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::RightStretchInverse() const
{
    // get the right stretch tensor
    mat3ds C = RightCauchyGreen();
    double l2[3];
    vec3d v[3];
    C.eigen2(l2, v);
    mat3ds U = dyad(v[0])/sqrt(l2[0]) + dyad(v[1])/sqrt(l2[1]) + dyad(v[2])/sqrt(l2[2]);
    
    return U;
}

//-----------------------------------------------------------------------------
//! Calculates the left stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::LeftStretchInverse() const
{
    // get the left stretch tensor
    mat3ds B = LeftCauchyGreen();
    double l2[3];
    vec3d v[3];
    B.eigen2(l2, v);
    mat3ds V = dyad(v[0])/sqrt(l2[0]) + dyad(v[1])/sqrt(l2[1]) + dyad(v[2])/sqrt(l2[2]);
    
    return V;
}

//-----------------------------------------------------------------------------
//! Calculates the right stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::RightHencky() const
{
    // get the right stretch tensor
    mat3ds C = RightCauchyGreen();
    double l2[3];
    vec3d v[3];
    C.eigen2(l2, v);
    mat3ds H = dyad(v[0])*log(l2[0])/2 + dyad(v[1])*log(l2[1])/2 + dyad(v[2])*log(l2[2])/2;
    
    return H;
}

//-----------------------------------------------------------------------------
//! Calculates the left stretch tensor at the current material point

mat3ds FEElasticMaterialPoint::LeftHencky() const
{
    // get the left stretch tensor
    mat3ds B = LeftCauchyGreen();
    double l2[3];
    vec3d v[3];
    B.eigen2(l2, v);
    mat3ds h = dyad(v[0])*log(l2[0])/2 + dyad(v[1])*log(l2[1])/2 + dyad(v[2])*log(l2[2])/2;
    
    return h;
}

//-----------------------------------------------------------------------------
//! Calculates the rotation tensor at the current material point

mat3d FEElasticMaterialPoint::Rotation() const
{
    return m_F*RightStretchInverse();
}

//-----------------------------------------------------------------------------
//! Calculates the Lagrangian strain at the current material point

mat3ds FEElasticMaterialPoint::Strain() const
{
	// get the right Cauchy-Green tensor
	// C = Ft*F
	const mat3d& F = m_F;
	mat3ds C;
	C.xx() = F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]; // = C[0][0]
	C.yy() = F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]; // = C[1][1]
	C.zz() = F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]; // = C[2][2]
	C.xy() = F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]; // = C[0][1]
	C.yz() = F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]; // = C[1][2]
	C.xz() = F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]; // = C[0][2]

	// calculate the Lagrangian strain
	// E = 1/2*(C - 1)
	mat3ds E = (C - mat3dd(1))*0.5;

	return E;
}

//-----------------------------------------------------------------------------
//! Calculates the small-strain tensor from the deformation gradient
mat3ds FEElasticMaterialPoint::SmallStrain() const
{
	// caculate small strain tensor
	const mat3d& F = m_F;
	return mat3ds(F[0][0] - 1.0, F[1][1] - 1.0, F[2][2] - 1.0, 0.5*(F[0][1] + F[1][0]), 0.5*(F[0][2] + F[2][0]), 0.5*(F[1][2] + F[2][1]));
}

//-----------------------------------------------------------------------------
//! Calculates the 2nd PK stress from the Cauchy stress stored in the point

mat3ds FEElasticMaterialPoint::pull_back(const mat3ds& A) const
{
	const mat3d& F = m_F;
	double J = m_J;
	mat3d Fi = F.inverse();
	mat3d P = Fi*A;

	return mat3ds(J*(P[0][0]*Fi[0][0]+P[0][1]*Fi[0][1]+P[0][2]*Fi[0][2]),
				  J*(P[1][0]*Fi[1][0]+P[1][1]*Fi[1][1]+P[1][2]*Fi[1][2]),
				  J*(P[2][0]*Fi[2][0]+P[2][1]*Fi[2][1]+P[2][2]*Fi[2][2]),
				  J*(P[0][0]*Fi[1][0]+P[0][1]*Fi[1][1]+P[0][2]*Fi[1][2]),
				  J*(P[1][0]*Fi[2][0]+P[1][1]*Fi[2][1]+P[1][2]*Fi[2][2]),
				  J*(P[0][0]*Fi[2][0]+P[0][1]*Fi[2][1]+P[0][2]*Fi[2][2]));
}

//-----------------------------------------------------------------------------
mat3ds FEElasticMaterialPoint::push_forward(const mat3ds& A) const
{
	const mat3d& F = m_F;
	mat3d P = F*A;
	double Ji = 1 / m_J;

	return mat3ds(Ji*(P[0][0]*F[0][0]+P[0][1]*F[0][1]+P[0][2]*F[0][2]),
				  Ji*(P[1][0]*F[1][0]+P[1][1]*F[1][1]+P[1][2]*F[1][2]),
				  Ji*(P[2][0]*F[2][0]+P[2][1]*F[2][1]+P[2][2]*F[2][2]),
				  Ji*(P[0][0]*F[1][0]+P[0][1]*F[1][1]+P[0][2]*F[1][2]),
				  Ji*(P[1][0]*F[2][0]+P[1][1]*F[2][1]+P[1][2]*F[2][2]),
				  Ji*(P[0][0]*F[2][0]+P[0][1]*F[2][1]+P[0][2]*F[2][2]));
}

//-----------------------------------------------------------------------------
// This function converts the spatial tangent to the material tangent
tens4ds FEElasticMaterialPoint::pull_back(const tens4ds& c) const
{
	const mat3d& F = m_F;
	mat3d Fi = F.inverse();
	double J = F.det();

	double C[6][6] = {0};

	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
				{
					C[0][0] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[0][l]*c(i,j,k,l);
					C[0][1] += J*Fi[0][i]*Fi[0][j]*Fi[1][k]*Fi[1][l]*c(i,j,k,l);
					C[0][2] += J*Fi[0][i]*Fi[0][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[0][3] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[0][4] += J*Fi[0][i]*Fi[0][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[0][5] += J*Fi[0][i]*Fi[0][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[1][1] += J*Fi[1][i]*Fi[1][j]*Fi[1][k]*Fi[1][l]*c(i,j,k,l);
					C[1][2] += J*Fi[1][i]*Fi[1][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[1][3] += J*Fi[1][i]*Fi[1][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[1][4] += J*Fi[1][i]*Fi[1][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[1][5] += J*Fi[1][i]*Fi[1][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[2][2] += J*Fi[2][i]*Fi[2][j]*Fi[2][k]*Fi[2][l]*c(i,j,k,l);
					C[2][3] += J*Fi[2][i]*Fi[2][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[2][4] += J*Fi[2][i]*Fi[2][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[2][5] += J*Fi[2][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[3][3] += J*Fi[0][i]*Fi[1][j]*Fi[0][k]*Fi[1][l]*c(i,j,k,l);
					C[3][4] += J*Fi[0][i]*Fi[1][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[3][5] += J*Fi[0][i]*Fi[1][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[4][4] += J*Fi[1][i]*Fi[2][j]*Fi[1][k]*Fi[2][l]*c(i,j,k,l);
					C[4][5] += J*Fi[1][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);

					C[5][5] += J*Fi[0][i]*Fi[2][j]*Fi[0][k]*Fi[2][l]*c(i,j,k,l);
				}

	return tens4ds(C);
}

//-----------------------------------------------------------------------------
// This function converts the material tangent to the spatial tangent
tens4ds FEElasticMaterialPoint::push_forward(const tens4ds& C) const
{
	const mat3d& F = m_F;
	double Ji = 1/F.det();
	double c[6][6] = {0};

	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
				{
					c[0][0] += Ji*F[0][i]*F[0][j]*F[0][k]*F[0][l]*C(i,j,k,l);
					c[0][1] += Ji*F[0][i]*F[0][j]*F[1][k]*F[1][l]*C(i,j,k,l);
					c[0][2] += Ji*F[0][i]*F[0][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[0][3] += Ji*F[0][i]*F[0][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[0][4] += Ji*F[0][i]*F[0][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[0][5] += Ji*F[0][i]*F[0][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[1][1] += Ji*F[1][i]*F[1][j]*F[1][k]*F[1][l]*C(i,j,k,l);
					c[1][2] += Ji*F[1][i]*F[1][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[1][3] += Ji*F[1][i]*F[1][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[1][4] += Ji*F[1][i]*F[1][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[1][5] += Ji*F[1][i]*F[1][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[2][2] += Ji*F[2][i]*F[2][j]*F[2][k]*F[2][l]*C(i,j,k,l);
					c[2][3] += Ji*F[2][i]*F[2][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[2][4] += Ji*F[2][i]*F[2][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[2][5] += Ji*F[2][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[3][3] += Ji*F[0][i]*F[1][j]*F[0][k]*F[1][l]*C(i,j,k,l);
					c[3][4] += Ji*F[0][i]*F[1][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[3][5] += Ji*F[0][i]*F[1][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[4][4] += Ji*F[1][i]*F[2][j]*F[1][k]*F[2][l]*C(i,j,k,l);
					c[4][5] += Ji*F[1][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);

					c[5][5] += Ji*F[0][i]*F[2][j]*F[0][k]*F[2][l]*C(i,j,k,l);
				}

	return tens4ds(c);
}
