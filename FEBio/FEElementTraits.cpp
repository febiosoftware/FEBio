// FEElementTraits.cpp: implementation of the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementTraits.h"
#include "FEElement.h"
#include "FEException.h"

//-----------------------------------------------------------------------------
//! Unpack solid element traits

void FESolidElementTraits::UnpackData(int nflag)
{
	double J[3][3], Ji[3][3];//, F[3][3];
	double det, deti;

	int i, n;

	double *Grn, *Gsn, *Gtn;
//	double GX, GY, GZ;

	for (n=0; n<nint; ++n)
	{
		Grn = Gr[n];
		Gsn = Gs[n];
		Gtn = Gt[n];

		// --- calculate jacobian Jt ---
		if (nflag & FE_UNPACK_JACT)
		{
			J[0][0] = J[0][1] = J[0][2] = 0.0;
			J[1][0] = J[1][1] = J[1][2] = 0.0;
			J[2][0] = J[2][1] = J[2][2] = 0.0;
			for (i=0; i<neln; ++i)
			{
				const double& Gri = Grn[i];
				const double& Gsi = Gsn[i];
				const double& Gti = Gtn[i];

				const double& x = rt[i].x;
				const double& y = rt[i].y;
				const double& z = rt[i].z;

				J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
				J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
				J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
			}

			// calculate the determinant
			det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

			// if det(J) > 0 calculate the inverse of the jacobian
			if (det > 0)
			{
				deti = 1.0 / det;

				Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
				Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
				Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
	
				Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
				Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
				Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);

				Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
				Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
				Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
			}
			else throw NegativeJacobian(m_pel->m_nID, n+1, det);

			matrix3_copy(m_Jt [n],  J);
			matrix3_copy(m_Jti[n],  Ji);
			m_detJt[n] = det;
		}

		// --- calculate jacobian J0 ---
		if ((nflag & FE_UNPACK_JAC0) || (nflag & FE_UNPACK_DEFGRAD))
		{
			J[0][0] = J[0][1] = J[0][2] = 0.0;
			J[1][0] = J[1][1] = J[1][2] = 0.0;
			J[2][0] = J[2][1] = J[2][2] = 0.0;
			for (i=0; i<neln; ++i)
			{
				const double& Gri = Grn[i];
				const double& Gsi = Gsn[i];
				const double& Gti = Gtn[i];

				const double& x = r0[i].x;
				const double& y = r0[i].y;
				const double& z = r0[i].z;

				J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
				J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
				J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
			}

			// calculate the determinant
			det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

			// if det(J) > 0 calculate the inverse of the jacobian
			if (det > 0)
			{
				deti = 1.0 / det;

				Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
				Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
				Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
	
				Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
				Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
				Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
	
				Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
				Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
				Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
			}
			else throw NegativeJacobian(m_pel->m_nID, n+1, det);

			if (nflag & FE_UNPACK_JAC0)
			{
				matrix3_copy(m_J0 [n],  J);
				matrix3_copy(m_J0i[n],  Ji);
				m_detJ0[n] = det;
			}
		}

		// --- calculate deformation gradient F ---
		// note that we need Ji0 for the deformation gradient!

/*		if (nflag & FE_UNPACK_DEFGRAD)
		{
			F[0][0] = F[0][1] = F[0][2] = 0;
			F[1][0] = F[1][1] = F[1][2] = 0;
			F[2][0] = F[2][1] = F[2][2] = 0;
			for (i=0; i<neln; ++i)
			{
				const double& Gri = Grn[i];
				const double& Gsi = Gsn[i];
				const double& Gti = Gtn[i];

				const double& x = rt[i].x;
				const double& y = rt[i].y;
				const double& z = rt[i].z;

				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !
				GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
				GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
				GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
	
				// calculate deformation gradient F
				F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
				F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
				F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
			}

			det =	(F[0][0]*(F[1][1]*F[2][2] - F[1][2]*F[2][1]) + 
					 F[0][1]*(F[1][2]*F[2][0] - F[2][2]*F[1][0]) + 
					 F[0][2]*(F[1][0]*F[2][1] - F[1][1]*F[2][0]));

			matrix3_copy(m_F[n],  F);
			m_detF[n] = det;
		}
*/
	}
}

//-----------------------------------------------------------------------------
//! Unpack shell element traits

void FEShellElementTraits::UnpackData(int nflag)
{
	double J[3][3], Ji[3][3];//, F[3][3];
	double det, deti;

	int i, n;

	double *Hrn, *Hsn, *Hn;
//	double NX, NY, NZ;
//	double MX, MY, MZ;

	FEShellElement* pe = dynamic_cast<FEShellElement*>(m_pel);
	double* h0 = pe->m_h0;
	double za;

	for (n=0; n<nint; ++n)
	{
		Hrn = Hr[n];
		Hsn = Hs[n];
		Hn = H[n];

		// --- calculate jacobian Jt ---
		if (nflag & FE_UNPACK_JACT)
		{
			J[0][0] = J[0][1] = J[0][2] = 0.0;
			J[1][0] = J[1][1] = J[1][2] = 0.0;
			J[2][0] = J[2][1] = J[2][2] = 0.0;
			for (i=0; i<neln; ++i)
			{
				const double& Hri = Hrn[i];
				const double& Hsi = Hsn[i];
				const double& Hi = Hn[i];

				const double& x = rt[i].x;
				const double& y = rt[i].y;
				const double& z = rt[i].z;

				const double& dx = Dt[i].x;
				const double& dy = Dt[i].y;
				const double& dz = Dt[i].z;

				za = 0.5*gt[n]*h0[i];

				J[0][0] += Hri*x + Hri*za*dx; J[0][1] += Hsi*x + Hsi*za*dx; J[0][2] += 0.5*h0[i]*Hi*dx;
				J[1][0] += Hri*y + Hri*za*dy; J[1][1] += Hsi*y + Hsi*za*dy; J[1][2] += 0.5*h0[i]*Hi*dy;
				J[2][0] += Hri*z + Hri*za*dz; J[2][1] += Hsi*z + Hsi*za*dz; J[2][2] += 0.5*h0[i]*Hi*dz;
			}

			// calculate the determinant
			det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

			// if det(J) > 0 calculate the inverse of the jacobian
			if (det > 0)
			{
				deti = 1.0 / det;

				Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
				Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
				Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
	
				Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
				Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
				Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);

				Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
				Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
				Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
			}
			else 
			{
				throw NegativeJacobian(m_pel->m_nID, n+1, det, m_pel);
			}

			matrix3_copy(m_Jt [n],  J);
			matrix3_copy(m_Jti[n],  Ji);
			m_detJt[n] = det;
		}

		// --- calculate jacobian J0 ---
		if ((nflag & FE_UNPACK_JAC0) || (nflag & FE_UNPACK_DEFGRAD))
		{
			J[0][0] = J[0][1] = J[0][2] = 0.0;
			J[1][0] = J[1][1] = J[1][2] = 0.0;
			J[2][0] = J[2][1] = J[2][2] = 0.0;
			for (i=0; i<neln; ++i)
			{
				const double& Hri = Hrn[i];
				const double& Hsi = Hsn[i];
				const double& Hi = Hn[i];

				const double& x = r0[i].x;
				const double& y = r0[i].y;
				const double& z = r0[i].z;

				const double& dx = D0[i].x;
				const double& dy = D0[i].y;
				const double& dz = D0[i].z;

				za = 0.5*gt[n]*h0[i];

				J[0][0] += Hri*x + Hri*za*dx; J[0][1] += Hsi*x + Hsi*za*dx; J[0][2] += 0.5*h0[i]*Hi*dx;
				J[1][0] += Hri*y + Hri*za*dy; J[1][1] += Hsi*y + Hsi*za*dy; J[1][2] += 0.5*h0[i]*Hi*dy;
				J[2][0] += Hri*z + Hri*za*dz; J[2][1] += Hsi*z + Hsi*za*dz; J[2][2] += 0.5*h0[i]*Hi*dz;
			}

			// calculate the determinant
			det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) 
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0]) 
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

			// if det(J) > 0 calculate the inverse of the jacobian
			if (det > 0)
			{
				deti = 1.0 / det;

				Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
				Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
				Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
	
				Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
				Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
				Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
	
				Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
				Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
				Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
			}
			else throw NegativeJacobian(m_pel->m_nID, n+1, det);

			if (nflag & FE_UNPACK_JAC0)
			{
				matrix3_copy(m_J0 [n],  J);
				matrix3_copy(m_J0i[n],  Ji);
				m_detJ0[n] = det;
			}
		}

		// --- calculate deformation gradient F ---
		// note that we need Ji0 for the deformation gradient!

/*		if (nflag & FE_UNPACK_DEFGRAD)
		{
			F[0][0] = F[0][1] = F[0][2] = 0;
			F[1][0] = F[1][1] = F[1][2] = 0;
			F[2][0] = F[2][1] = F[2][2] = 0;
			for (i=0; i<neln; ++i)
			{
				const double& Hri = Hrn[i];
				const double& Hsi = Hsn[i];
				const double& Hi  = Hn[i];

				const double& x = rt[i].x;
				const double& y = rt[i].y;
				const double& z = rt[i].z;

				const double& dx = Dt[i].x;
				const double& dy = Dt[i].y;
				const double& dz = Dt[i].z;

				za = 0.5*gt[n]*h0[i];

				// calculate global gradient of shape functions
				// note that we need the transposed of Ji, not Ji itself !
				NX = Ji[0][0]*Hri+Ji[1][0]*Hsi;
				NY = Ji[0][1]*Hri+Ji[1][1]*Hsi;
				NZ = Ji[0][2]*Hri+Ji[1][2]*Hsi;

				MX = za*Ji[0][0]*Hri + za*Ji[1][0]*Hsi + Ji[2][0]*0.5*h0[i]*Hi;
				MY = za*Ji[0][1]*Hri + za*Ji[1][1]*Hsi + Ji[2][1]*0.5*h0[i]*Hi;
				MZ = za*Ji[0][2]*Hri + za*Ji[1][2]*Hsi + Ji[2][2]*0.5*h0[i]*Hi;

				// calculate deformation gradient F
				F[0][0] += NX*x + MX*dx; F[0][1] += NY*x + MY*dx; F[0][2] += NZ*x + MZ*dx;
				F[1][0] += NX*y + MX*dy; F[1][1] += NY*y + MY*dy; F[1][2] += NZ*y + MZ*dy;
				F[2][0] += NX*z + MX*dz; F[2][1] += NY*z + MY*dz; F[2][2] += NZ*z + MZ*dz;
			}

			det =	(F[0][0]*(F[1][1]*F[2][2] - F[1][2]*F[2][1]) + 
					 F[0][1]*(F[1][2]*F[2][0] - F[2][2]*F[1][0]) + 
					 F[0][2]*(F[1][0]*F[2][1] - F[1][1]*F[2][0]));

			matrix3_copy(m_F[n],  F);
			m_detF[n] = det;
		}
*/
	}
}

//*****************************************************************************
//                          F E H E X E L E M E N T
//*****************************************************************************

void FEHexElementTraits::init()
{
	int n;

	// integration point coordinates
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gt[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gt[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gt[2] = -a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gt[3] = -a; gw[3] = 1;
	gr[4] = -a; gs[4] = -a; gt[4] =  a; gw[4] = 1;
	gr[5] =  a; gs[5] = -a; gt[5] =  a; gw[5] = 1;
	gr[6] =  a; gs[6] =  a; gt[6] =  a; gw[6] = 1;
	gr[7] = -a; gs[7] =  a; gt[7] =  a; gw[7] = 1;

	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][1] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][2] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][3] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][4] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][5] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][6] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 + gt[n]);
		H[n][7] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 + gt[n]);
	}

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][1] =  0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][2] =  0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][3] = -0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][4] = -0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][5] =  0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][6] =  0.125*(1 + gs[n])*(1 + gt[n]);
		Gr[n][7] = -0.125*(1 + gs[n])*(1 + gt[n]);

		Gs[n][0] = -0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][1] = -0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][2] =  0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][3] =  0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][4] = -0.125*(1 - gr[n])*(1 + gt[n]);
		Gs[n][5] = -0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][6] =  0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][7] =  0.125*(1 - gr[n])*(1 + gt[n]);

		Gt[n][0] = -0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][1] = -0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][2] = -0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][3] = -0.125*(1 - gr[n])*(1 + gs[n]);
		Gt[n][4] =  0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][5] =  0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][6] =  0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][7] =  0.125*(1 - gr[n])*(1 + gs[n]);
	}
}

//*****************************************************************************
//                          F E R I H E X E L E M E N T
//*****************************************************************************

void FERIHexElementTraits::init()
{
	int n;

	// This is for a six point integration rule
	// integration point coordinates
	const double a = 8.0 / 6.0;
	gr[0] = -1; gs[0] = 0; gt[0] = 0; gw[0] = a;
	gr[1] =  1; gs[1] = 0; gt[1] = 0; gw[1] = a;
	gr[2] =  0; gs[2] =-1; gt[2] = 0; gw[2] = a;
	gr[3] =  0; gs[3] = 1; gt[3] = 0; gw[3] = a;
	gr[4] =  0; gs[4] = 0; gt[4] =-1; gw[4] = a;
	gr[5] =  0; gs[5] = 0; gt[5] = 1; gw[5] = a;

	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][1] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][2] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][3] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][4] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][5] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][6] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 + gt[n]);
		H[n][7] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 + gt[n]);
	}

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][1] =  0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][2] =  0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][3] = -0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][4] = -0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][5] =  0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][6] =  0.125*(1 + gs[n])*(1 + gt[n]);
		Gr[n][7] = -0.125*(1 + gs[n])*(1 + gt[n]);

		Gs[n][0] = -0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][1] = -0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][2] =  0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][3] =  0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][4] = -0.125*(1 - gr[n])*(1 + gt[n]);
		Gs[n][5] = -0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][6] =  0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][7] =  0.125*(1 - gr[n])*(1 + gt[n]);

		Gt[n][0] = -0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][1] = -0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][2] = -0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][3] = -0.125*(1 - gr[n])*(1 + gs[n]);
		Gt[n][4] =  0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][5] =  0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][6] =  0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][7] =  0.125*(1 - gr[n])*(1 + gs[n]);
	}
}

//*****************************************************************************
//                          F E U D F H E X E L E M E N T
//*****************************************************************************

void FEUDFHexElementTraits::init()
{
	// single gauss-point integration rule
	gr[0] = 0; gs[0] = 0; gt[0] = 0; gw[0] = 8.0;

	// calculate shape function values at gauss points
	H[0][0] = 0.125*(1 - gr[0])*(1 - gs[0])*(1 - gt[0]);
	H[0][1] = 0.125*(1 + gr[0])*(1 - gs[0])*(1 - gt[0]);
	H[0][2] = 0.125*(1 + gr[0])*(1 + gs[0])*(1 - gt[0]);
	H[0][3] = 0.125*(1 - gr[0])*(1 + gs[0])*(1 - gt[0]);
	H[0][4] = 0.125*(1 - gr[0])*(1 - gs[0])*(1 + gt[0]);
	H[0][5] = 0.125*(1 + gr[0])*(1 - gs[0])*(1 + gt[0]);
	H[0][6] = 0.125*(1 + gr[0])*(1 + gs[0])*(1 + gt[0]);
	H[0][7] = 0.125*(1 - gr[0])*(1 + gs[0])*(1 + gt[0]);

	// calculate local derivatives of shape functions at gauss points
	Gr[0][0] = -0.125*(1 - gs[0])*(1 - gt[0]);
	Gr[0][1] =  0.125*(1 - gs[0])*(1 - gt[0]);
	Gr[0][2] =  0.125*(1 + gs[0])*(1 - gt[0]);
	Gr[0][3] = -0.125*(1 + gs[0])*(1 - gt[0]);
	Gr[0][4] = -0.125*(1 - gs[0])*(1 + gt[0]);
	Gr[0][5] =  0.125*(1 - gs[0])*(1 + gt[0]);
	Gr[0][6] =  0.125*(1 + gs[0])*(1 + gt[0]);
	Gr[0][7] = -0.125*(1 + gs[0])*(1 + gt[0]);

	Gs[0][0] = -0.125*(1 - gr[0])*(1 - gt[0]);
	Gs[0][1] = -0.125*(1 + gr[0])*(1 - gt[0]);
	Gs[0][2] =  0.125*(1 + gr[0])*(1 - gt[0]);
	Gs[0][3] =  0.125*(1 - gr[0])*(1 - gt[0]);
	Gs[0][4] = -0.125*(1 - gr[0])*(1 + gt[0]);
	Gs[0][5] = -0.125*(1 + gr[0])*(1 + gt[0]);
	Gs[0][6] =  0.125*(1 + gr[0])*(1 + gt[0]);
	Gs[0][7] =  0.125*(1 - gr[0])*(1 + gt[0]);

	Gt[0][0] = -0.125*(1 - gr[0])*(1 - gs[0]);
	Gt[0][1] = -0.125*(1 + gr[0])*(1 - gs[0]);
	Gt[0][2] = -0.125*(1 + gr[0])*(1 + gs[0]);
	Gt[0][3] = -0.125*(1 - gr[0])*(1 + gs[0]);
	Gt[0][4] =  0.125*(1 - gr[0])*(1 - gs[0]);
	Gt[0][5] =  0.125*(1 + gr[0])*(1 - gs[0]);
	Gt[0][6] =  0.125*(1 + gr[0])*(1 + gs[0]);
	Gt[0][7] =  0.125*(1 - gr[0])*(1 + gs[0]);
}

//*****************************************************************************
//                          F E T E T E L E M E N T
//*****************************************************************************

void FETetElementTraits::init()
{
	int n;

	// gaussian integration for tetrahedral elements
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 1.0 / 24.0;

	gr[0] = b; gs[0] = b; gt[0] = b; gw[0] = w;
	gr[1] = a; gs[1] = b; gt[1] = b; gw[1] = w;
	gr[2] = b; gs[2] = a; gt[2] = b; gw[2] = w;
	gr[3] = b; gs[3] = b; gt[3] = a; gw[3] = w;

	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1 - gr[n] - gs[n] - gt[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
		H[n][3] = gt[n];
	}

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;
		Gr[n][3] =  0;

		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
		Gs[n][3] =  0;

		Gt[n][0] = -1;
		Gt[n][1] =  0;
		Gt[n][2] =  0;
		Gt[n][3] =  1;
	}
}

//*****************************************************************************
//                          F E P E N T A E L E M E N T
//*****************************************************************************

void FEPentaElementTraits::init()
{
	int n;

	//gauss intergration points
	const double a = 1.0/6.0;
	const double b = 2.0/3.0;
	const double c = 1.0 / sqrt(3.0);
	const double w = 1.0 / 6.0;

	gr[0] = a; gs[0] = a; gt[0] = -c; gw[0] = w;
	gr[1] = b; gs[1] = a; gt[1] = -c; gw[1] = w;
	gr[2] = a; gs[2] = b; gt[2] = -c; gw[2] = w;
	gr[3] = a; gs[3] = a; gt[3] =  c; gw[3] = w;
	gr[4] = b; gs[4] = a; gt[4] =  c; gw[4] = w;
	gr[5] = a; gs[5] = b; gt[5] =  c; gw[5] = w;

	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.5*(1 - gt[n])*(1 - gr[n] - gs[n]);
		H[n][1] = 0.5*(1 - gt[n])*gr[n];
		H[n][2] = 0.5*(1 - gt[n])*gs[n];
		H[n][3] = 0.5*(1 + gt[n])*(1 - gr[n] - gs[n]);
		H[n][4] = 0.5*(1 + gt[n])*gr[n];
		H[n][5] = 0.5*(1 + gt[n])*gs[n];
	}

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.5*(1 - gt[n]);
		Gr[n][1] =  0.5*(1 - gt[n]);
		Gr[n][2] =  0.0;
		Gr[n][3] = -0.5*(1 + gt[n]);
		Gr[n][4] =  0.5*(1 + gt[n]);
		Gr[n][5] =  0.0;

		Gs[n][0] = -0.5*(1 - gt[n]);
		Gs[n][1] =  0.0;
		Gs[n][2] =  0.5*(1 - gt[n]);
		Gs[n][3] = -0.5*(1 + gt[n]);
		Gs[n][4] =  0.0;
		Gs[n][5] =  0.5*(1 + gt[n]);

		Gt[n][0] = -0.5*(1 - gr[n] - gs[n]);
		Gt[n][1] = -0.5*gr[n];
		Gt[n][2] = -0.5*gs[n];
		Gt[n][3] =  0.5*(1 - gr[n] - gs[n]);
		Gt[n][4] =  0.5*gr[n];
		Gt[n][5] =  0.5*gs[n];
	}
}

//*****************************************************************************
//                          F E Q U A D E L E M E N T
//*****************************************************************************

void FEQuadElementTraits::init()
{
	int n;

	const double a = 1.0 / sqrt(3.0);

	gr[0] = -a; gs[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gw[3] = 1;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}

	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.25*(1-gs[n]);
		Gr[n][1] =  0.25*(1-gs[n]);
		Gr[n][2] =  0.25*(1+gs[n]);
		Gr[n][3] = -0.25*(1+gs[n]);

		Gs[n][0] = -0.25*(1-gr[n]);
		Gs[n][1] = -0.25*(1+gr[n]);
		Gs[n][2] =  0.25*(1+gr[n]);
		Gs[n][3] =  0.25*(1-gr[n]);
	}
}

//*****************************************************************************
//                          F E N I Q U A D E L E M E N T
//*****************************************************************************

void FENIQuadElementTraits::init()
{
	int n;

	gr[0] = -1; gs[0] = -1; gw[0] = 1;
	gr[1] =  1; gs[1] = -1; gw[1] = 1;
	gr[2] =  1; gs[2] =  1; gw[2] = 1;
	gr[3] = -1; gs[3] =  1; gw[3] = 1;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}

	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.25*(1-gs[n]);
		Gr[n][1] =  0.25*(1-gs[n]);
		Gr[n][2] =  0.25*(1+gs[n]);
		Gr[n][3] = -0.25*(1+gs[n]);

		Gs[n][0] = -0.25*(1-gr[n]);
		Gs[n][1] = -0.25*(1+gr[n]);
		Gs[n][2] =  0.25*(1+gr[n]);
		Gs[n][3] =  0.25*(1-gr[n]);
	}
}

//*****************************************************************************
//                          F E T R I E L E M E N T
//*****************************************************************************

void FETriElementTraits::init()
{
	int n;

	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;

	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}

	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;

		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
	}
}

//*****************************************************************************
//                          F E N I T R I E L E M E N T
//*****************************************************************************

void FENITriElementTraits::init()
{
	int n;

	const double a = 1.0 / 6.0;

	gr[0] = 0; gs[0] = 0; gw[0] = a;
	gr[1] = 1; gs[1] = 0; gw[1] = a;
	gr[2] = 0; gs[2] = 1; gw[2] = a;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}

	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;

		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
	}
}

//*****************************************************************************
//                          F E S H E L L Q U A D E L E M E N T
//*****************************************************************************

void FEShellQuadElementTraits::init()
{
	int n;

	const double a = 1.0 / sqrt(3.0);
	const double b = sqrt(3.0/5.0);
	const double w = 5.0 / 9.0;

	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -b; gw[ 0] = w;
	gr[ 1] =  a; gs[ 1] = -a; gt[ 1] = -b; gw[ 1] = w;
	gr[ 2] =  a; gs[ 2] =  a; gt[ 2] = -b; gw[ 2] = w;
	gr[ 3] = -a; gs[ 3] =  a; gt[ 3] = -b; gw[ 3] = w;

	gr[ 4] = -a; gs[ 4] = -a; gt[ 4] =  0; gw[ 4] = w;
	gr[ 5] =  a; gs[ 5] = -a; gt[ 5] =  0; gw[ 5] = w;
	gr[ 6] =  a; gs[ 6] =  a; gt[ 6] =  0; gw[ 6] = w;
	gr[ 7] = -a; gs[ 7] =  a; gt[ 7] =  0; gw[ 7] = w;

	gr[ 8] = -a; gs[ 8] = -a; gt[ 8] =  b; gw[ 8] = w;
	gr[ 9] =  a; gs[ 9] = -a; gt[ 9] =  b; gw[ 9] = w;
	gr[10] =  a; gs[10] =  a; gt[10] =  b; gw[10] = w;
	gr[11] = -a; gs[11] =  a; gt[11] =  b; gw[11] = w;


	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -0.25*(1-gs[n]);
		Hr[n][1] =  0.25*(1-gs[n]);
		Hr[n][2] =  0.25*(1+gs[n]);
		Hr[n][3] = -0.25*(1+gs[n]);

		Hs[n][0] = -0.25*(1-gr[n]);
		Hs[n][1] = -0.25*(1+gr[n]);
		Hs[n][2] =  0.25*(1+gr[n]);
		Hs[n][3] =  0.25*(1-gr[n]);
	}
}

//*****************************************************************************
//                          F E S H E L L T R I E L E M E N T
//*****************************************************************************

void FEShellTriElementTraits::init()
{
	int n;

	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	const double c = sqrt(3.0/5.0);
	const double w = 5.0 / 9.0;

	gr[0] = a; gs[0] = a; gt[0] = -b; gw[0] = a*w;
	gr[1] = b; gs[1] = a; gt[1] = -b; gw[1] = a*w;
	gr[2] = a; gs[2] = b; gt[2] = -b; gw[2] = a*w;

	gr[3] = a; gs[3] = a; gt[3] =  0; gw[3] = a*w;
	gr[4] = b; gs[4] = a; gt[4] =  0; gw[4] = a*w;
	gr[5] = a; gs[5] = b; gt[5] =  0; gw[5] = a*w;

	gr[6] = a; gs[6] = a; gt[6] =  b; gw[6] = a*w;
	gr[7] = b; gs[7] = a; gt[7] =  b; gw[7] = a*w;
	gr[8] = a; gs[8] = b; gt[8] =  b; gw[8] = a*w;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -1;
		Hr[n][1] =  1;
		Hr[n][2] =  0;

		Hs[n][0] = -1;
		Hs[n][1] =  0;
		Hs[n][2] =  1;
	}
}
