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
#include "FENeoHookeanTransIso.h"


// define the material parameters
BEGIN_FECORE_CLASS(FENeoHookeanTransIso, FEElasticMaterial)
	ADD_PARAMETER(m_Ep, "Ep");
	ADD_PARAMETER(m_Ez, "Ez");
	ADD_PARAMETER(m_vz, "vz");
	ADD_PARAMETER(m_vp, "vp");
	ADD_PARAMETER(m_gz, "gz");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// CompNeoHookean_Transiso
//////////////////////////////////////////////////////////////////////

mat3ds FENeoHookeanTransIso::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	//define jacobian
	double J = detF;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();

	// calculate left Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double B[3][3];
	double B2[3][3];

	B[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	B[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	B[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

 	B[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	B[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	B[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

  	B[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
	B[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	B[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// calculate square of B
	// (we commented out the matrix components we do not need)
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// calculate Ba
	vec3d Ba;
	Ba.x = B[0][0]*a.x + B[0][1]*a.y + B[0][2]*a.z;
	Ba.y = B[1][0]*a.x + B[1][1]*a.y + B[1][2]*a.z;
	Ba.z = B[2][0]*a.x + B[2][1]*a.y + B[2][2]*a.z;

	//compute invariants
	double I1,I3,I4,I5;

	I1 = B[0][0]+B[1][1]+B[2][2];
	I3=J*J;
	I4 = lam*lam;;
	I5 = I4*(a*Ba);

	// compute lame parameters
	double m_gp=m_Ep/(2*(1+m_vp));
	double m_lamp=-((m_Ep*(m_Ez*m_vp + m_Ep*m_vz*m_vz))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz)));
	double m_lamz=(m_Ez*m_Ez*(-1 + m_vp*m_vp) - m_Ep*(m_Ep + 4*(m_gp + m_gz)*(1 + m_vp))*m_vz*m_vz + m_Ez*(2*m_gz - m_Ep*m_vp - 2*m_gz*m_vp*m_vp - 2*m_gp*(-1 + m_vp*m_vp) + 2*m_Ep*m_vz + 2*m_Ep*m_vp*m_vz))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz));
	double m_alpha=(m_Ep*(m_Ep*m_vz*m_vz - m_Ez*(m_vp*(-1 + m_vz) + m_vz)))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz));


	// calculate stress
	mat3ds s;

	//trans iso
	s.xx()=(2*(B[0][0]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + a.x*(a.x*B[0][0] + a.y*B[0][1] + a.z*B[0][2])*(-m_gp + m_gz)*I4 + I3*(-m_gp/(2.*I3) + ((-1 + J)*m_lamp)/(2.*J)) + a.x*a.x*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;
	s.xy()=(2*(B[0][1]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + ((a.y*(a.x*B[0][0] + a.y*B[0][1] + a.z*B[0][2]) + a.x*(a.x*B[0][1] + a.y*B[1][1] + a.z*B[1][2]))*(-m_gp + m_gz)*I4)/2. + a.x*a.y*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;
	s.xz()=(2*(B[0][2]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + ((a.z*(a.x*B[0][0] + a.y*B[0][1] + a.z*B[0][2]) + a.x*(a.x*B[0][2] + a.y*B[1][2] + a.z*B[2][2]))*(-m_gp + m_gz)*I4)/2. + a.x*a.z*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;
	s.yy()=(2*(B[1][1]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + a.y*(a.x*B[0][1] + a.y*B[1][1] + a.z*B[1][2])*(-m_gp + m_gz)*I4 + I3*(-m_gp/(2.*I3) + ((-1 + J)*m_lamp)/(2.*J)) + a.y*a.y*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;
	s.yz()=(2*(B[1][2]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + ((a.z*(a.x*B[0][1] + a.y*B[1][1] + a.z*B[1][2]) + a.y*(a.x*B[0][2] + a.y*B[1][2] + a.z*B[2][2]))*(-m_gp + m_gz)*I4)/2. + a.y*a.z*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;
	s.zz()=(2*(B[2][2]*(m_gp/2. + (m_alpha*(-1 + I4))/4.) + a.z*(a.x*B[0][2] + a.y*B[1][2] + a.z*B[2][2])*(-m_gp + m_gz)*I4 + I3*(-m_gp/(2.*I3) + ((-1 + J)*m_lamp)/(2.*J)) + a.z*a.z*I4*(m_gz/2. + (m_alpha*(-3 + I1))/4. - m_gz/(2.*I4) - (-m_gp + m_gz)*I4 + ((-1 + lam)*m_lamz)/(2.*lam))))/J;


	//compressible neohookean (for testing--comment out)
	//s.x=(2*((B[0][0]*m_gp)/2. + I3*(-m_gp/(2.*I3) + (m_lamp*log(J))/(2.*I3))))/J;
	//s.xy=(B[0][1]*m_gp)/J;
	//s.xz=(B[0][2]*m_gp)/J;
	//s.y=(2*((B[1][1]*m_gp)/2. + I3*(-m_gp/(2.*I3) + (m_lamp*log(J))/(2.*I3))))/J;
	//s.yz=(B[1][2]*m_gp)/J;
	//s.z=(2*((B[2][2]*m_gp)/2. + I3*(-m_gp/(2.*I3) + (m_lamp*log(J))/(2.*I3))))/J;


	return s;

}

tens4ds FENeoHookeanTransIso::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double detF = pt.m_J;

	//define jacobian
	double J = detF;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();

	// calculate left Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double B[3][3];
	double B2[3][3];

	B[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	B[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	B[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

 	B[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	B[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	B[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

  	B[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
	B[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	B[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// calculate square of B
	// (we commented out the matrix components we do not need)
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// calculate Ba
	vec3d Ba;
	Ba.x = B[0][0]*a.x + B[0][1]*a.y + B[0][2]*a.z;
	Ba.y = B[1][0]*a.x + B[1][1]*a.y + B[1][2]*a.z;
	Ba.z = B[2][0]*a.x + B[2][1]*a.y + B[2][2]*a.z;

	//compute invariants
	double I1,I3,I4,I5;

	I1 = B[0][0]+B[1][1]+B[2][2];
	I3=J*J;
	I4 = lam*lam;;
	I5 = I4*(a*Ba);

	// lame parameters
	double m_gp=m_Ep/(2*(1+m_vp));
	double m_lamp=-((m_Ep*(m_Ez*m_vp + m_Ep*m_vz*m_vz))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz)));
	double m_lamz=(m_Ez*m_Ez*(-1 + m_vp*m_vp) - m_Ep*(m_Ep + 4*(m_gp + m_gz)*(1 + m_vp))*m_vz*m_vz + m_Ez*(2*m_gz - m_Ep*m_vp - 2*m_gz*m_vp*m_vp - 2*m_gp*(-1 + m_vp*m_vp) + 2*m_Ep*m_vz + 2*m_Ep*m_vp*m_vz))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz));
	double m_alpha=(m_Ep*(m_Ep*m_vz*m_vz - m_Ez*(m_vp*(-1 + m_vz) + m_vz)))/((1 + m_vp)*(m_Ez*(-1 + m_vp) + 2*m_Ep*m_vz*m_vz));


	//trans iso
	double D[6][6] = {0};
	D[0][0]=(2*a.x*a.x*B[0][0]*(m_alpha + 2*m_gz)*I4 + m_gp*(2 - 4*a.x*a.x*B[0][0]*I4 + 4*a.x*a.x*a.x*a.x*I4*I4) + J*m_lamp + a.x*a.x*a.x*a.x*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz))/J;
	D[0][1]=(a.y*a.y*m_alpha*B[0][0]*I4 + 4*a.x*a.y*B[0][1]*(-m_gp + m_gz)*I4 - J*m_lamp + 2*I3*m_lamp + a.x*a.x*(m_alpha*B[1][1]*I4 + a.y*a.y*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[0][2]=(a.z*a.z*m_alpha*B[0][0]*I4 + 4*a.x*a.z*B[0][2]*(-m_gp + m_gz)*I4 - J*m_lamp + 2*I3*m_lamp + a.x*a.x*(m_alpha*B[2][2]*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[0][3]=(a.x*(a.y*B[0][0]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*B[0][1]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*a.x*a.y*(4*m_gp*I4*I4 + m_gz*(2 - 4*I4*I4) + lam*m_lamz)))/J;
	D[0][4]=(a.x*(a.z*B[0][0]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*B[0][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*a.x*a.z*(4*m_gp*I4*I4 + m_gz*(2 - 4*I4*I4) + lam*m_lamz)))/J;
	D[0][5]=(a.y*a.z*m_alpha*B[0][0]*I4 - 2*a.x*(a.z*B[0][1] + a.y*B[0][2])*(m_gp - m_gz)*I4 + a.x*a.x*(m_alpha*B[1][2]*I4 + a.y*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[1][1]=(2*a.y*a.y*B[1][1]*(m_alpha + 2*m_gz)*I4 + m_gp*(2 - 4*a.y*a.y*B[1][1]*I4 + 4*a.y*a.y*a.y*a.y*I4*I4) + J*m_lamp + a.y*a.y*a.y*a.y*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz))/J;
	D[1][2]=(a.z*a.z*m_alpha*B[1][1]*I4 + 4*a.y*a.z*B[1][2]*(-m_gp + m_gz)*I4 - J*m_lamp + 2*I3*m_lamp + a.y*a.y*(m_alpha*B[2][2]*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[1][3]=(a.y*(a.y*B[0][1]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*(B[1][1]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.y*a.y*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz))))/J;
	D[1][4]=(a.y*(a.y*m_alpha*B[0][2] + 2*a.z*B[0][1]*(-m_gp + m_gz))*I4 + a.x*(a.z*m_alpha*B[1][1]*I4 + 2*a.y*B[1][2]*(-m_gp + m_gz)*I4 + a.y*a.y*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[1][5]=(a.y*(a.z*B[1][1]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.y*B[1][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.y*a.y*a.z*(4*m_gp*I4*I4 + m_gz*(2 - 4*I4*I4) + lam*m_lamz)))/J;
	D[2][2]=(2*a.z*a.z*B[2][2]*(m_alpha + 2*m_gz)*I4 + m_gp*(2 - 4*a.z*a.z*B[2][2]*I4 + 4*a.z*a.z*a.z*a.z*I4*I4) + J*m_lamp + a.z*a.z*a.z*a.z*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz))/J;
	D[2][3]=(a.z*(a.z*m_alpha*B[0][1] + 2*a.y*B[0][2]*(-m_gp + m_gz))*I4 + a.x*(2*a.z*B[1][2]*(-m_gp + m_gz)*I4 + a.y*(m_alpha*B[2][2]*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz))))/J;
	D[2][4]=(a.z*(a.z*B[0][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.x*(B[2][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz))))/J;
	D[2][5]=(a.z*(a.z*B[1][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.y*(B[2][2]*(m_alpha - 2*m_gp + 2*m_gz)*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz))))/J;
	D[3][3]=(a.y*a.y*B[0][0]*m_gz*I4 + 2*a.x*a.y*B[0][1]*(m_alpha + m_gz)*I4 - m_gp*(-1 + 2*a.x*a.y*B[0][1]*I4 + a.x*a.x*B[1][1]*I4 + a.y*a.y*I4*(B[0][0] - 4*a.x*a.x*I4)) + J*m_lamp - I3*m_lamp + a.x*a.x*(B[1][1]*m_gz*I4 + a.y*a.y*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[3][4]=(a.y*a.z*B[0][0]*(-m_gp + m_gz)*I4 + a.x*(a.z*B[0][1] + a.y*B[0][2])*(m_alpha - m_gp + m_gz)*I4 + a.x*a.x*(B[1][2]*(-m_gp + m_gz)*I4 + a.y*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[3][5]=(a.y*(a.y*B[0][2]*(-m_gp + m_gz) + a.z*B[0][1]*(m_alpha - m_gp + m_gz))*I4 + a.x*(a.z*B[1][1]*(-m_gp + m_gz)*I4 + a.y*B[1][2]*(m_alpha - m_gp + m_gz)*I4 + a.y*a.y*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[4][4]=(a.z*a.z*B[0][0]*m_gz*I4 + 2*a.x*a.z*B[0][2]*(m_alpha + m_gz)*I4 - m_gp*(-1 + 2*a.x*a.z*B[0][2]*I4 + a.x*a.x*B[2][2]*I4 + a.z*a.z*I4*(B[0][0] - 4*a.x*a.x*I4)) + J*m_lamp - I3*m_lamp + a.x*a.x*(B[2][2]*m_gz*I4 + a.z*a.z*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz)))/J;
	D[4][5]=(a.z*(a.z*B[0][1]*(-m_gp + m_gz) + a.y*B[0][2]*(m_alpha - m_gp + m_gz))*I4 + a.x*(a.z*B[1][2]*(m_alpha - m_gp + m_gz)*I4 + a.y*(B[2][2]*(-m_gp + m_gz)*I4 + a.z*a.z*(2*m_gz + 4*m_gp*I4*I4 - 4*m_gz*I4*I4 + lam*m_lamz))))/J;
	D[5][5]=(a.z*a.z*B[1][1]*m_gz*I4 + 2*a.y*a.z*B[1][2]*(m_alpha + m_gz)*I4 - m_gp*(-1 + 2*a.y*a.z*B[1][2]*I4 + a.y*a.y*B[2][2]*I4 + a.z*a.z*I4*(B[1][1] - 4*a.y*a.y*I4)) + J*m_lamp - I3*m_lamp + a.y*a.y*(B[2][2]*m_gz*I4 + a.z*a.z*(2*m_gz - 4*m_gz*I4*I4 + lam*m_lamz)))/J;


	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];

/*
	FILE* fp = fopen("C.txt", "wt");
	int i, j;
	for (i=0; i<6; ++i)
	{
		for (j=0; j<6; ++j) fprintf(fp, "%15.7lg ", D[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
*/

	//the following code is used for testing purposes and replaces the above stiffness when testing
	//comment out for normal use

	//constant trans iso
	//D[0][0]=2*m_gp + m_lamp;
	//D[0][1]=m_lamp;
	//D[0][2]=m_alpha + m_lamp;
	//D[0][3]=0;
	//D[0][4]=0;
	//D[0][5]=0;
	//D[1][1]=2*m_gp + m_lamp;
	//D[1][2]=m_alpha + m_lamp;
	//D[1][3]=0;
	//D[1][4]=0;
	//D[1][5]=0;
	//D[2][2]=2*m_alpha + 2*m_gp + 2*m_gz + m_lamp + m_lamz;
	//D[2][3]=0;
	//D[2][4]=0;
	//D[2][5]=0;
	//D[3][3]=m_gp;
	//D[3][4]=0;
	//D[3][5]=0;
	//D[4][4]=m_gz;
	//D[4][5]=0;
	//D[5][5]=m_gz;


	//compressible neohookean
	//D[0][0]=(2*m_gp + m_lamp - m_lamp*log(I3))/J;
	//D[0][1]=m_lamp/J;
	//D[0][2]=m_lamp/J;
	//D[0][3]=0;
	//D[0][4]=0;
	//D[0][5]=0;
	//D[1][1]=(2*m_gp + m_lamp - m_lamp*log(I3))/J;
	//D[1][2]=m_lamp/J;
	//D[1][3]=0;
	//D[1][4]=0;
	//D[1][5]=0;
	//D[2][2]=(2*m_gp + m_lamp - m_lamp*log(I3))/J;
	//D[2][3]=0;
	//D[2][4]=0;
	//D[2][5]=0;
	//D[3][3]=(m_gp - (m_lamp*log(I3))/2.)/J;
	//D[3][4]=0;
	//D[3][5]=0;
	//D[4][4]=(m_gp - (m_lamp*log(I3))/2.)/J;
	//D[4][5]=0;
	//D[5][5]=(m_gp - (m_lamp*log(I3))/2.)/J;


	//identity
	//D[0][0]=1;
	//D[0][1]=0;
	//D[0][2]=0;
	//D[0][3]=0;
	//D[0][4]=0;
	//D[0][5]=0;
	//D[1][1]=1;
	//D[1][2]=0;
	//D[1][3]=0;
	//D[1][4]=0;
	//D[1][5]=0;
	//D[2][2]=1;
	//D[2][3]=0;
	//D[2][4]=0;
	//D[2][5]=0;
	//D[3][3]=1;
	//D[3][4]=0;
	//D[3][5]=0;
	//D[4][4]=1;
	//D[4][5]=0;
	//D[5][5]=1;

	return tens4ds(D);
}
