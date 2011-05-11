#include "stdafx.h"
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
// Fiber material stress
//
mat3ds FEFiberMaterial::Stress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double J = pt.J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0/3.0);
	double twoJi = 2.0*Ji;

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = 0;
	if (lamd > 1)
	{
		double lamdi = 1.0/lamd;
		double Wl;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1-1))-1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
		}
		W4  = 0.5*lamdi*Wl;
	}
	else 
	{
		W4 = 0;
	}	

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// ---
	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);
	
	// calculate stress: 
	mat3ds s = T.dev()*twoJi;

	// --- active contraction contribution ---
	if (m_lcna >= 0)
	{
		double ctenslm = m_plc->Value();

		// current sarcomere length
		double strl = m_refl*lamd;

		// sarcomere length change
		double dl = strl - m_l0;

		if (dl >= 0)
		{
			double eca50i = (exp(m_beta*dl) - 1);

			// active fiber stress
			double saf = ctenslm*eca50i / ( eca50i + 1 );

			// add saf*(a x a)
			s += AxA*saf;
		}
	}
	//-----------------------------------------

	return s;
}

//-----------------------------------------------------------------------------
// Fiber material tangent
//
tens4ds FEFiberMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// get initial local material axis
	vec3d a0;
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	double W4, W44;
	if (lamd >= 1)
	{
		double lamdi = 1.0/lamd;
		double Wl, Wll;
		if (lamd < m_lam1)
		{
			Wl  = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
			Wll = m_c3*lamdi*(m_c4*exp(m_c4*(lamd - 1)) - lamdi*(exp(m_c4*(lamd-1))-1));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1-1))-1) - m_c5*m_lam1;
			Wl  = lamdi*(m_c5*lamd + c6);
			Wll = -c6*lamdi*lamdi;
		}
		W4  = 0.5*lamdi*Wl;
		W44 = 0.25*lamdi*lamdi*(Wll - lamdi*Wl);
	}
	else 
	{
		W4 = 0;
		W44 = 0;
	}

	// --- active contraction stiffness ---
/*	if (lcna >= 0)
	{
		double ctenslm = m_plc->Value();

		// current sarcomere length
		double strl = refl*lamd;

		// sarcomere length change
		double dl = strl - l0;

		if (dl >= 0) W44 += J*2*beta*refl*exp(-beta*dl);
	}
*/

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4  = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0/3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI/3.0)*(4.0/3.0*Ji*WC) + IxI*(4.0/9.0*Ji*CWWC) + cw;

	return c;
}
