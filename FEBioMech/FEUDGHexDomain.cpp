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
#include "FEUDGHexDomain.h"
#include "FEElasticMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEUDGHexDomain, FEElasticSolidDomain)
	ADD_PARAMETER(m_hg, "hg");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEUDGHexDomain::FEUDGHexDomain(FEModel* pfem) : FEElasticSolidDomain(pfem)
{ 
	m_hg = 1.0; 
}

//-----------------------------------------------------------------------------
//! Set the hourglass parameter
void FEUDGHexDomain::SetHourGlassParameter(double hg)
{
	m_hg = hg;
}

//-----------------------------------------------------------------------------
bool FEUDGHexDomain::Create(int nelems, FE_Element_Spec spec)
{
	// make sure these solid, hex8 elements are requested
	if (spec.eclass == FE_Element_Class::FE_ELEM_INVALID_CLASS) spec.eclass = FE_Element_Class::FE_ELEM_SOLID;
	else if (spec.eclass != FE_Element_Class::FE_ELEM_SOLID) return false;

	if (spec.eshape == FE_Element_Shape::FE_ELEM_INVALID_SHAPE) spec.eshape = FE_Element_Shape::ET_HEX8;
	else if (spec.eshape != FE_Element_Shape::ET_HEX8) return false;

	// we need to enforce HEX8G1 integration rule
	spec.etype = FE_HEX8G1;

	// now allocate the domain
	return FESolidDomain::Create(nelems, spec);
}

//-----------------------------------------------------------------------------
void FEUDGHexDomain::InternalForces(FEGlobalVector& R)
{
	int NE = (int)m_Elem.size();
#pragma omp parallel for
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

		// get the element force vector and initialize it to zero
		int ndof = 3*el.Nodes();

		// element force vector
		vector<double> fe;
		fe.assign(ndof, 0);

		// calculate internal force vector
		UDGInternalForces(el, fe);

		// get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		R.Assemble(el.m_node, lm, fe);
	}
}


//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for enhanced strain
//! solid elements.

void FEUDGHexDomain::UDGInternalForces(FESolidElement& el, vector<double>& fe)
{
	// get the stress data
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	mat3ds& s = pt.m_s;

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx, Gy, Gz;

	// calculate the internal force
	for (int i=0; i<8; ++i)
	{
		Gx = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];

		fe[3*i  ] -= ve*(Gx*s.xx() + Gy*s.xy() + Gz*s.xz());
		fe[3*i+1] -= ve*(Gx*s.xy() + Gy*s.yy() + Gz*s.yz());
		fe[3*i+2] -= ve*(Gx*s.xz() + Gy*s.yz() + Gz*s.zz());
	}

	// add hourglass forces
	UDGHourglassForces(el, fe);
}


//-----------------------------------------------------------------------------
//! calculates the hourglass forces

void FEUDGHexDomain::UDGHourglassForces(FESolidElement &el, vector<double> &fe)
{
	int i;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	int neln = el.Nodes();

	vec3d r0[8], rt[8];
	for (i=0; i<neln; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	mat3d F;
	AvgDefGrad(el, F, GX, GY, GZ);

	double u4 = 0, u5 = 0, u6 = 0, u7 = 0;
	double v4 = 0, v5 = 0, v6 = 0, v7 = 0;
	double w4 = 0, w5 = 0, w6 = 0, w7 = 0;

	u4 = x4 - (F[0][0]*X4 + F[0][1]*Y4 + F[0][2]*Z4);
	v4 = y4 - (F[1][0]*X4 + F[1][1]*Y4 + F[1][2]*Z4);
	w4 = z4 - (F[2][0]*X4 + F[2][1]*Y4 + F[2][2]*Z4);

	u5 = x5 - (F[0][0]*X5 + F[0][1]*Y5 + F[0][2]*Z5);
	v5 = y5 - (F[1][0]*X5 + F[1][1]*Y5 + F[1][2]*Z5);
	w5 = z5 - (F[2][0]*X5 + F[2][1]*Y5 + F[2][2]*Z5);

	u6 = x6 - (F[0][0]*X6 + F[0][1]*Y6 + F[0][2]*Z6);
	v6 = y6 - (F[1][0]*X6 + F[1][1]*Y6 + F[1][2]*Z6);
	w6 = z6 - (F[2][0]*X6 + F[2][1]*Y6 + F[2][2]*Z6);

	u7 = x7 - (F[0][0]*X7 + F[0][1]*Y7 + F[0][2]*Z7);
	v7 = y7 - (F[1][0]*X7 + F[1][1]*Y7 + F[1][2]*Z7);
	w7 = z7 - (F[2][0]*X7 + F[2][1]*Y7 + F[2][2]*Z7);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass forces
	for (i=0; i<8; ++i)
	{
		fe[3*i  ] -= m_hg*(g4[i]*u4 + g5[i]*u5 + g6[i]*u6 + g7[i]*u7);
		fe[3*i+1] -= m_hg*(g4[i]*v4 + g5[i]*v5 + g6[i]*v6 + g7[i]*v7);
		fe[3*i+2] -= m_hg*(g4[i]*w4 + g5[i]*w5 + g6[i]*w6 + g7[i]*w7);
	}
}

void FEUDGHexDomain::StiffnessMatrix(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();

	vector<int> lm;

	// repeat over all solid elements
	int NE = (int)m_Elem.size();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate material stiffness
		UDGMaterialStiffness(el, ke);

		// calculate geometrical stiffness
		UDGGeometricalStiffness(el, ke);

		// add hourglass stiffness
		UDGHourglassStiffness(fem, el, ke);

		// assign symmetic parts
		// TODO: Can this be omitted by changing the Assemble routine so that it only
		// grabs elements from the upper diagonal matrix?
		for (int i=0; i<ndof; ++i)
			for (int j=i+1; j<ndof; ++j)
				ke[j][i] = ke[i][j];

		// get the element's LM vector
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates the hourglass stiffness for UDG hex elements

void FEUDGHexDomain::UDGHourglassStiffness(FEModel& fem, FESolidElement& el, matrix& ke)
{
	int i, j;

	const double h4[8] = { 1,-1, 1,-1, 1,-1, 1,-1 };
	const double h5[8] = { 1,-1,-1, 1,-1, 1, 1,-1 };
	const double h6[8] = { 1, 1,-1,-1,-1,-1, 1, 1 };
	const double h7[8] = {-1, 1,-1, 1, 1,-1, 1,-1 };

	int neln = el.Nodes();

	vec3d r0[8], rt[8];
	for (i=0; i<neln; ++i)
	{
		r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
		rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	double x4 = 0, x5 = 0, x6 = 0, x7 = 0;
	double y4 = 0, y5 = 0, y6 = 0, y7 = 0;
	double z4 = 0, z5 = 0, z6 = 0, z7 = 0;

	double X4 = 0, X5 = 0, X6 = 0, X7 = 0;
	double Y4 = 0, Y5 = 0, Y6 = 0, Y7 = 0;
	double Z4 = 0, Z5 = 0, Z6 = 0, Z7 = 0;

	for (i=0; i<8; ++i)
	{
		X4 += h4[i]*r0[i].x; Y4 += h4[i]*r0[i].y; Z4 += h4[i]*r0[i].z;
		X5 += h5[i]*r0[i].x; Y5 += h5[i]*r0[i].y; Z5 += h5[i]*r0[i].z;
		X6 += h6[i]*r0[i].x; Y6 += h6[i]*r0[i].y; Z6 += h6[i]*r0[i].z;
		X7 += h7[i]*r0[i].x; Y7 += h7[i]*r0[i].y; Z7 += h7[i]*r0[i].z;

		x4 += h4[i]*rt[i].x; y4 += h4[i]*rt[i].y; z4 += h4[i]*rt[i].z;
		x5 += h5[i]*rt[i].x; y5 += h5[i]*rt[i].y; z5 += h5[i]*rt[i].z;
		x6 += h6[i]*rt[i].x; y6 += h6[i]*rt[i].y; z6 += h6[i]*rt[i].z;
		x7 += h7[i]*rt[i].x; y7 += h7[i]*rt[i].y; z7 += h7[i]*rt[i].z;
	}

	FEMesh& mesh = *GetMesh();

	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	double g4[8] = {0}, g5[8] = {0}, g6[8] = {0}, g7[8] = {0};

	for (i=0; i<8; ++i)
	{
		g4[i] = h4[i] - (GX[i]*X4 + GY[i]*Y4 + GZ[i]*Z4);
		g5[i] = h5[i] - (GX[i]*X5 + GY[i]*Y5 + GZ[i]*Z5);
		g6[i] = h6[i] - (GX[i]*X6 + GY[i]*Y6 + GZ[i]*Z6);
		g7[i] = h7[i] - (GX[i]*X7 + GY[i]*Y7 + GZ[i]*Z7);
	}

	// calculate hourglass stiffness
	for (i=0; i<8; ++i)
	{
		for (j=i; j<8; ++j)
		{
			double kab = m_hg*(g4[i]*g4[j] + g5[i]*g5[j] + g6[i]*g6[j] + g7[i]*g7[j]);

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
	}
}


//-----------------------------------------------------------------------------

void FEUDGHexDomain::UDGGeometricalStiffness(FESolidElement& el, matrix& ke)
{
	int i, j;

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// get the material point data
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// element's Cauchy-stress tensor at gauss point n
	// s is the voight vector
	mat3ds& s = pt.m_s;

	FEMesh& mesh = *GetMesh();

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	for (i=0; i<8; ++i)
		for (j=i; j<8; ++j)
		{
			kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
				   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
				   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*ve;

			ke[3*i  ][3*j  ] += kab;
			ke[3*i+1][3*j+1] += kab;
			ke[3*i+2][3*j+2] += kab;
		}
}


//-----------------------------------------------------------------------------

void FEUDGHexDomain::UDGMaterialStiffness(FESolidElement &el, matrix &ke)
{
	// make sure we have the right element type
	assert(el.Type() == FE_HEX8G1);

	int i, i3, j, j3;

	// Get the current element's data
	const int neln = el.Nodes();

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	FEMesh& mesh = *GetMesh();

	// calculate the average cartesian derivatives
	double GX[8], GY[8], GZ[8];
	AvgCartDerivs(el, GX, GY, GZ);

	// calculate average deformation gradient Fbar
	mat3d Fb;
	AvgDefGrad(el, Fb, GX, GY, GZ);

	// calculate the transposed inverse of Fbar
	mat3d Fti = Fb.transinv();

	// calculate current element volume
	double ve = HexVolume(el, 1);

	// current averaged shape derivatives
	double Gx[8], Gy[8], Gz[8];
	for (i=0; i<8; ++i)
	{
		Gx[i] = Fti(0,0)*GX[i]+Fti(0,1)*GY[i]+Fti(0,2)*GZ[i];
		Gy[i] = Fti(1,0)*GX[i]+Fti(1,1)*GY[i]+Fti(1,2)*GZ[i];
		Gz[i] = Fti(2,0)*GX[i]+Fti(2,1)*GY[i]+Fti(2,2)*GZ[i];
	}

	// setup the material point
	// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
	FEMaterialPoint& mp = *el.GetMaterialPoint(0);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

	// get the 'D' matrix
	tens4ds C = m_pMat->Tangent(mp);
	C.extract(D);

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

			ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*ve;
			ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*ve;
			ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*ve;

			ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*ve;
			ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*ve;
			ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*ve;

			ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*ve;
			ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*ve;
			ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*ve;
		}
	}
}

//-----------------------------------------------------------------------------
void FEUDGHexDomain::Update(const FETimeInfo& tp)
{
	int nint, neln;
	double* gw;
	vec3d r0[8], rt[8];

	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the number of integration points
		nint = el.GaussPoints();

		// number of nodes
		neln = el.Nodes();

		// nodal coordinates
		for (int j=0; j<neln; ++j)
		{
			r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
			rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
		}

		// get the integration weights
		gw = el.GaussWeights();

		// for the enhanced strain hex we need a slightly different procedure
		// for calculating the element's stress. For this element, the stress
		// is evaluated using an average deformation gradient.

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(0);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		mp.m_r0 = el.Evaluate(r0, 0);
		mp.m_rt = el.Evaluate(rt, 0);

		// get the average cartesian derivatives
		double GX[8], GY[8], GZ[8];
		AvgCartDerivs(el, GX, GY, GZ);

		// get the average deformation gradient and determinant
		AvgDefGrad(el, pt.m_F, GX, GY, GZ);
		pt.m_J = pt.m_F.det();

		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the average Cartesian derivatives
//! Note that we assume that the GX, GY and GX contain the averaged 
//! Cartesian derivatives

void FEUDGHexDomain::AvgDefGrad(FESolidElement& el, mat3d& F, double GX[8], double GY[8], double GZ[8])
{
	vec3d rt[8];
	for (int j=0; j<8; ++j) rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;

	F.zero();
	for (int i=0; i<8; ++i)
	{
		F[0][0] += rt[i].x*GX[i];
		F[0][1] += rt[i].x*GY[i];
		F[0][2] += rt[i].x*GZ[i];

		F[1][0] += rt[i].y*GX[i];
		F[1][1] += rt[i].y*GY[i];
		F[1][2] += rt[i].y*GZ[i];

		F[2][0] += rt[i].z*GX[i];
		F[2][1] += rt[i].z*GY[i];
		F[2][2] += rt[i].z*GZ[i];
	}
}

//-----------------------------------------------------------------------------
//! Calculates the average Cartesian derivatives

void FEUDGHexDomain::AvgCartDerivs(FESolidElement& el, double GX[8], double GY[8], double GZ[8], int nstate)
{
	// get the nodal coordinates
	int neln = el.Nodes();
	vec3d r[8];
	if (nstate == 0)
	{
		for (int i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	}
	else
	{
		for (int i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	double x1 = r[0].x, y1 = r[0].y, z1 = r[0].z;
	double x2 = r[1].x, y2 = r[1].y, z2 = r[1].z;
	double x3 = r[2].x, y3 = r[2].y, z3 = r[2].z;
	double x4 = r[3].x, y4 = r[3].y, z4 = r[3].z;
	double x5 = r[4].x, y5 = r[4].y, z5 = r[4].z;
	double x6 = r[5].x, y6 = r[5].y, z6 = r[5].z;
	double x7 = r[6].x, y7 = r[6].y, z7 = r[6].z;
	double x8 = r[7].x, y8 = r[7].y, z8 = r[7].z;

	const double f12 = 1.0/12.0;

	// set up the B-matrix
	// we use the G arrays to store the B-matrix
	GX[0] = f12*(y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-(z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+y8*(z4-z5));
	GY[0] = f12*(z2*((x6-x3)-(x4-x5))+z3*(x2-x4)+z4*((x3-x8)-(x5-x2))+z5*((x8-x6)-(x2-x4))+z6*(x5-x2)+z8*(x4-x5));
	GZ[0] = f12*(x2*((y6-y3)-(y4-y5))+x3*(y2-y4)+x4*((y3-y8)-(y5-y2))+x5*((y8-y6)-(y2-y4))+x6*(y5-y2)+x8*(y4-y5));

	GX[1] = f12*(y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-(z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+y5*(z1-z6));
	GY[1] = f12*(z3*((x7-x4)-(x1-x6))+z4*(x3-x1)+z1*((x4-x5)-(x6-x3))+z6*((x5-x7)-(x3-x1))+z7*(x6-x3)+z5*(x1-x6));
	GZ[1] = f12*(x3*((y7-y4)-(y1-y6))+x4*(y3-y1)+x1*((y4-y5)-(y6-y3))+x6*((y5-y7)-(y3-y1))+x7*(y6-y3)+x5*(y1-y6));

	GX[2] = f12*(y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-(z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+y6*(z2-z7));
	GY[2] = f12*(z4*((x8-x1)-(x2-x7))+z1*(x4-x2)+z2*((x1-x6)-(x7-x4))+z7*((x6-x8)-(x4-x2))+z8*(x7-x4)+z6*(x2-x7));
	GZ[2] = f12*(x4*((y8-y1)-(y2-y7))+x1*(y4-y2)+x2*((y1-y6)-(y7-y4))+x7*((y6-y8)-(y4-y2))+x8*(y7-y4)+x6*(y2-y7));

	GX[3] = f12*(y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-(z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+y7*(z3-z8));
	GY[3] = f12*(z1*((x5-x2)-(x3-x8))+z2*(x1-x3)+z3*((x2-x7)-(x8-x1))+z8*((x7-x5)-(x1-x3))+z5*(x8-x1)+z7*(x3-x8));
	GZ[3] = f12*(x1*((y5-y2)-(y3-y8))+x2*(y1-y3)+x3*((y2-y7)-(y8-y1))+x8*((y7-y5)-(y1-y3))+x5*(y8-y1)+x7*(y3-y8));

	GX[4] = f12*(y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-(z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+y2*(z6-z1));
	GY[4] = f12*(z8*((x4-x7)-(x6-x1))+z7*(x8-x6)+z6*((x7-x2)-(x1-x8))+z1*((x2-x4)-(x8-x6))+z4*(x1-x8)+z2*(x6-x1));
	GZ[4] = f12*(x8*((y4-y7)-(y6-y1))+x7*(y8-y6)+x6*((y7-y2)-(y1-y8))+x1*((y2-y4)-(y8-y6))+x4*(y1-y8)+x2*(y6-y1));

	GX[5] = f12*(y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-(z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+y3*(z7-z2));
	GY[5] = f12*(z5*((x1-x8)-(x7-x2))+z8*(x5-x7)+z7*((x8-x3)-(x2-x5))+z2*((x3-x1)-(x5-x7))+z1*(x2-x5)+z3*(x7-x2));
	GZ[5] = f12*(x5*((y1-y8)-(y7-y2))+x8*(y5-y7)+x7*((y8-y3)-(y2-y5))+x2*((y3-y1)-(y5-y7))+x1*(y2-y5)+x3*(y7-y2));

	GX[6] = f12*(y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-(z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+y4*(z8-z3));
	GY[6] = f12*(z6*((x2-x5)-(x8-x3))+z5*(x6-x8)+z8*((x5-x4)-(x3-x6))+z3*((x4-x2)-(x6-x8))+z2*(x3-x6)+z4*(x8-x3));
	GZ[6] = f12*(x6*((y2-y5)-(y8-y3))+x5*(y6-y8)+x8*((y5-y4)-(y3-y6))+x3*((y4-y2)-(y6-y8))+x2*(y3-y6)+x4*(y8-y3));

	GX[7] = f12*(y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-(z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+y1*(z5-z4));
	GY[7] = f12*(z7*((x3-x6)-(x5-x4))+z6*(x7-x5)+z5*((x6-x1)-(x4-x7))+z4*((x1-x3)-(x7-x5))+z3*(x4-x7)+z1*(x5-x4));
	GZ[7] = f12*(x7*((y3-y6)-(y5-y4))+x6*(y7-y5)+x5*((y6-y1)-(y4-y7))+x4*((y1-y3)-(y7-y5))+x3*(y4-y7)+x1*(y5-y4));

	// calculate the volume
	double Vi = 1./(x1*GX[0]+x2*GX[1]+x3*GX[2]+x4*GX[3]+x5*GX[4]+x6*GX[5]+x7*GX[6]+x8*GX[7]);

	// divide the B-matrix by the volume
	GX[0] *= Vi; GY[0] *= Vi; GZ[0] *= Vi;
	GX[1] *= Vi; GY[1] *= Vi; GZ[1] *= Vi;
	GX[2] *= Vi; GY[2] *= Vi; GZ[2] *= Vi;
	GX[3] *= Vi; GY[3] *= Vi; GZ[3] *= Vi;
	GX[4] *= Vi; GY[4] *= Vi; GZ[4] *= Vi;
	GX[5] *= Vi; GY[5] *= Vi; GZ[5] *= Vi;
	GX[6] *= Vi; GY[6] *= Vi; GZ[6] *= Vi;
	GX[7] *= Vi; GY[7] *= Vi; GZ[7] *= Vi;
}
//-----------------------------------------------------------------------------
//! Calculates the exact volume of a hexahedral element

double FEUDGHexDomain::HexVolume(FESolidElement& el, int state)
{
	// let's make sure this is indeed a hex element
//	assert(el.Type() == FE_HEX8G8);

	int neln = el.Nodes();
	vec3d r[8];
	if (state == 0)
	{
		for (int i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_r0;
	}
	else
	{
		for (int i=0; i<neln; ++i) r[i] = m_pMesh->Node(el.m_node[i]).m_rt;
	}

	// get the nodal coordinates
	double x1 = r[0].x, y1 = r[0].y, z1 = r[0].z;
	double x2 = r[1].x, y2 = r[1].y, z2 = r[1].z;
	double x3 = r[2].x, y3 = r[2].y, z3 = r[2].z;
	double x4 = r[3].x, y4 = r[3].y, z4 = r[3].z;
	double x5 = r[4].x, y5 = r[4].y, z5 = r[4].z;
	double x6 = r[5].x, y6 = r[5].y, z6 = r[5].z;
	double x7 = r[6].x, y7 = r[6].y, z7 = r[6].z;
	double x8 = r[7].x, y8 = r[7].y, z8 = r[7].z;

	// set up the B-matrix
	double B1[8];
//	double B2[8];
//	double B3[8];

	const double f12 = 1.0/12.0;

	B1[0] = f12*(y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-(z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+y8*(z4-z5));
//	B2[0] = f12*(z2*((x6-x3)-(x4-x5))+z3*(x2-x4)+z4*((x3-x8)-(x5-x2))+z5*((x8-x6)-(x2-x4))+z6*(x5-x2)+z8*(x4-x5));
//	B3[0] = f12*(x2*((y6-y3)-(y4-y5))+x3*(y2-y4)+x4*((y3-y8)-(y5-y2))+x5*((y8-y6)-(y2-y4))+x6*(y5-y2)+x8*(y4-y5));

	B1[1] = f12*(y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-(z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+y5*(z1-z6));
//	B2[1] = f12*(z3*((x7-x4)-(x1-x6))+z4*(x3-x1)+z1*((x4-x5)-(x6-x3))+z6*((x5-x7)-(x3-x1))+z7*(x6-x3)+z5*(x1-x6));
//	B3[1] = f12*(x3*((y7-y4)-(y1-y6))+x4*(y3-y1)+x1*((y4-y5)-(y6-y3))+x6*((y5-y7)-(y3-y1))+x7*(y6-y3)+x5*(y1-y6));

	B1[2] = f12*(y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-(z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+y6*(z2-z7));
//	B2[2] = f12*(z4*((x8-x1)-(x2-x7))+z1*(x4-x2)+z2*((x1-x6)-(x7-x4))+z7*((x6-x8)-(x4-x2))+z8*(x7-x4)+z6*(x2-x7));
//	B3[2] = f12*(x4*((y8-y1)-(y2-y7))+x1*(y4-y2)+x2*((y1-y6)-(y7-y4))+x7*((y6-y8)-(y4-y2))+x8*(y7-y4)+x6*(y2-y7));

	B1[3] = f12*(y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-(z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+y7*(z3-z8));
//	B2[3] = f12*(z1*((x5-x2)-(x3-x8))+z2*(x1-x3)+z3*((x2-x7)-(x8-x1))+z8*((x7-x5)-(x1-x3))+z5*(x8-x1)+z7*(x3-x8));
//	B3[3] = f12*(x1*((y5-y2)-(y3-y8))+x2*(y1-y3)+x3*((y2-y7)-(y8-y1))+x8*((y7-y5)-(y1-y3))+x5*(y8-y1)+x7*(y3-y8));

	B1[4] = f12*(y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-(z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+y2*(z6-z1));
//	B2[4] = f12*(z8*((x4-x7)-(x6-x1))+z7*(x8-x6)+z6*((x7-x2)-(x1-x8))+z1*((x2-x4)-(x8-x6))+z4*(x1-x8)+z2*(x6-x1));
//	B3[4] = f12*(x8*((y4-y7)-(y6-y1))+x7*(y8-y6)+x6*((y7-y2)-(y1-y8))+x1*((y2-y4)-(y8-y6))+x4*(y1-y8)+x2*(y6-y1));

	B1[5] = f12*(y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-(z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+y3*(z7-z2));
//	B2[5] = f12*(z5*((x1-x8)-(x7-x2))+z8*(x5-x7)+z7*((x8-x3)-(x2-x5))+z2*((x3-x1)-(x5-x7))+z1*(x2-x5)+z3*(x7-x2));
//	B3[5] = f12*(x5*((y1-y8)-(y7-y2))+x8*(y5-y7)+x7*((y8-y3)-(y2-y5))+x2*((y3-y1)-(y5-y7))+x1*(y2-y5)+x3*(y7-y2));

	B1[6] = f12*(y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-(z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+y4*(z8-z3));
//	B2[6] = f12*(z6*((x2-x5)-(x8-x3))+z5*(x6-x8)+z8*((x5-x4)-(x3-x6))+z3*((x4-x2)-(x6-x8))+z2*(x3-x6)+z4*(x8-x3));
//	B3[6] = f12*(x6*((y2-y5)-(y8-y3))+x5*(y6-y8)+x8*((y5-y4)-(y3-y6))+x3*((y4-y2)-(y6-y8))+x2*(y3-y6)+x4*(y8-y3));

	B1[7] = f12*(y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-(z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+y1*(z5-z4));
//	B2[7] = f12*(z7*((x3-x6)-(x5-x4))+z6*(x7-x5)+z5*((x6-x1)-(x4-x7))+z4*((x1-x3)-(x7-x5))+z3*(x4-x7)+z1*(x5-x4));
//	B3[7] = f12*(x7*((y3-y6)-(y5-y4))+x6*(y7-y5)+x5*((y6-y1)-(y4-y7))+x4*((y1-y3)-(y7-y5))+x3*(y4-y7)+x1*(y5-y4));

	// calculate the volume V= xi*B1[i] = yi*B2[i] = zi*B3[i] (sum over i)
	return (x1*B1[0]+x2*B1[1]+x3*B1[2]+x4*B1[3]+x5*B1[4]+x6*B1[5]+x7*B1[6]+x8*B1[7]);
}
