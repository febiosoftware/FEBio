// FEMultiphasic.cpp: implementation of the FEMultiphasic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMultiphasic.h"
#include "FECore/FEModel.h"
#include <complex>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEMultiphasic, FEMaterial)
	ADD_PARAMETER(m_phi0   , FE_PARAM_DOUBLE, "phi0"                );
	ADD_PARAMETER(m_rhoTw  , FE_PARAM_DOUBLE, "fluid_density"       );
	ADD_PARAMETER(m_cFr    , FE_PARAM_DOUBLE, "fixed_charge_density");
	ADD_PARAMETER(m_penalty, FE_PARAM_DOUBLE, "penalty"             );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Polynomial root solver

// function whose roots needs to be evaluated
void fn(complex<double>& z, complex<double>& fz, vector<double> a)
{
	int n = a.size()-1;
	fz = a[0];
	complex<double> x(1,0);
	
	for (int i=1; i<=n; ++i) {
		x *= z;
		fz += a[i]*x;
	}
	return;
}

// deflation
bool dflate(complex<double> zero, const int i, int& kount,
			complex<double>& fzero, complex<double>& fzrdfl,
			complex<double>* zeros, vector<double> a)
{
	complex<double> den;
	++kount;
	fn(zero, fzero, a);
	fzrdfl = fzero;
	if (i < 1) return false;
	for (int j=0; j<i; ++j) {
		den = zero - zeros[j];
		if (abs(den) == 0) {
			zeros[i] = zero*1.001;
			return true;
		} else {
			fzrdfl = fzrdfl/den;
		}
	}
	return false;
}

// Muller's method for solving roots of a function
void muller(bool fnreal, complex<double>* zeros, const int n, const int nprev,
			const int maxit, const double ep1, const double ep2, vector<double> a)
{
	int kount;
	complex<double> dvdf1p, fzrprv, fzrdfl, divdf1, divdf2;
	complex<double> fzr, zero, c, den, sqr, z;
	
	// initialization
	double eps1 = (ep1 > 1e-12) ? ep1:1e-12;
	double eps2 = (ep2 > 1e-20) ? ep2:1e-20;
	
	for (int i=nprev; i<n; ++i) {
		kount = 0;
	eloop:
		zero = zeros[i];
		complex<double> h = 0.5;
		complex<double> hprev = -1.0;
		
		// compute first three estimates for zero as
		// zero+0.5, zero-0.5, zero
		z = zero + 0.5;
		if (dflate(z, i, kount, fzr, dvdf1p, zeros, a)) goto eloop;
		z = zero - 0.5;
		if (dflate(z, i, kount, fzr, fzrprv, zeros, a)) goto eloop;
		dvdf1p = (fzrprv - dvdf1p)/hprev;
		if (dflate(zero, i, kount, fzr, fzrdfl, zeros, a)) goto eloop;
		do {
			divdf1 = (fzrdfl - fzrprv)/h;
			divdf2 = (divdf1 - dvdf1p)/(h+hprev);
			hprev = h;
			dvdf1p = divdf1;
			c = divdf1 + h*divdf2;
			sqr = c*c - 4.*fzrdfl*divdf2;
			if (fnreal && (sqr.real() < 0)) sqr = 0;
			sqr = sqrt(sqr);
			if (c.real()*sqr.real()+c.imag()*sqr.imag() < 0) {
				den = c - sqr;
			} else {
				den = c + sqr;
			}
			if (abs(den) <= 0.) den = 1.;
			h = -2.*fzrdfl/den;
			fzrprv = fzrdfl;
			zero = zero + h;
		dloop:
			fn(zero,fzrdfl,a);
			// check for convergence
			if (abs(h) < eps1*abs(zero)) break;
			if (abs(fzrdfl) < eps2) break;
			// check for divergence
			if (abs(fzrdfl) >= 10.*abs(fzrprv)) {
				h /= 2.;
				zero -= h;
				goto dloop;
			}
		} while (kount < maxit);
		zeros[i] = zero;
	}
	return;
}

//-----------------------------------------------------------------------------
//! FEMultiphasic constructor

FEMultiphasic::FEMultiphasic()
{	m_pPerm = 0;
	m_pOsmC = 0;
	m_phi0 = 0;
	m_rhoTw = 0;
	m_cFr = 0;
	m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
	m_penalty = 1;

	AddComponent<FEElasticMaterial      >(&m_pSolid    , "solid"              );
	AddComponent<FEHydraulicPermeability>(&m_pPerm     , "permeability"       );
	AddComponent<FEOsmoticCoefficient   >(&m_pOsmC     , "osmotic_coefficient");
//	AddComponent<FESolute               >(&m_pSolute[0], "solute",           0);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::Init()
{
	FEMaterial::Init();
	m_pSolid->Init();
	m_pPerm->Init();
	m_pOsmC->Init();
	for (int i=0; i<(int)m_pSolute.size(); ++i) m_pSolute[i]->Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 <= phi0 <= 1");
	if (m_rhoTw < 0) throw MaterialError("fluid_density must be positive");
	if (m_penalty < 0) throw MaterialError("penalty must be positive");
	
	// Determine how to solve for the electric potential psi
	int isol;
	int zmin = 0, zmax = 0, z;
	for (isol=0; isol<(int)m_pSolute.size(); ++isol) {
		z = m_pSolute[isol]->ChargeNumber();
		if (z < zmin) zmin = z;
		if (z > zmax) zmax = z;
	}
	m_zmin = zmin;
	m_ndeg = zmax - zmin;	// polynomial degree

	m_Rgas = FEModel::GetGlobalConstant("R");
	m_Tabs = FEModel::GetGlobalConstant("T");
	m_Fc = FEModel::GetGlobalConstant("Fc");
	
	if (m_Rgas <= 0) throw MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) throw MaterialError("A positive absolute temperature T must be defined in Globals section");
	if ((zmin || zmax) && (m_Fc <= 0)) throw MaterialError("A positive Faraday constant Fc must be defined in Globals section");
	
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEMultiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.J;
	// porosity
	//	double phiw = 1 - m_phi0/J;
	double phi0 = pet.m_phi0;
	double phiw = 1 - phi0/J;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! Fixed charge density in current configuration
double FEMultiphasic::FixedChargeDensity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// relative volume
	double J = et.J;
	double phi0 = pet.m_phi0;
	double cF = m_cFr*(1-phi0)/(J-phi0);
	
	return cF;
}

//-----------------------------------------------------------------------------
//! Electric potential
double FEMultiphasic::ElectricPotential(FEMaterialPoint& pt, const bool eform)
{
	// check if solution is neutral
	if (m_ndeg == 0) {
		if (eform) return 1.0;
		else return 0.0;
	}
	
	int i, j;
	
	// if not neutral, solve electroneutrality polynomial for zeta
	FESolutesMaterialPoint& set = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = m_pSolute.size();
	double cF = FixedChargeDensity(pt);
	vector<double> c(nsol);		// effective concentration
	vector<double> khat(nsol);	// solubility
	vector<int> z(nsol);		// charge number
	for (i=0; i<nsol; ++i) {
		c[i] = set.m_c[i];
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
		z[i] = m_pSolute[i]->ChargeNumber();
	}
	double zeta, psi;
	
	// evaluate polynomial coefficients
	const int n = m_ndeg;
	vector<double> a(n+1,0);
	if (m_zmin < 0) {
		for (i=0; i<nsol; ++i) {
			j = z[i] - m_zmin;
			a[j] += z[i]*khat[i]*c[i];
		}
		a[-m_zmin] = cF;
	} else {
		for (i=0; i<nsol; ++i) {
			j = z[i];
			a[j] += z[i]*khat[i]*c[i];
		}
		a[0] = cF;
	}

	// solve polynomial
	zeta = 1.0;
	if (n == 1) {
		if (a[1]) {
			zeta = -a[0]/a[1];
		}	// linear
	} else if (n==2) {
		if (a[2]) {
			zeta = (-a[1]+sqrt(a[1]*a[1]-4*a[0]*a[2]))/(2*a[2]);	// quadratic
		} else if (a[1]) {
			zeta = -a[0]/a[1];			// linear
		}
	} else {
		// solve higher degree polynomial using Muller's method
		bool fnreal = true;
		vector< complex<double> > zeros(n);
		int maxit = 100;
		double ep1 = 1e-6;
		double ep2 = 1e-12;
		
		muller(fnreal, &zeros[0], n, 0, maxit, ep1, ep2, a);
		for (i=0; i<n; ++i) {
			if (zeros[i].real() > 0) {
				zeta = zeros[i].real();
				break;
			}
		}
	}
	
	// Return exponential (non-dimensional) form if desired
	if (eform) return zeta;
	
	// Otherwise return dimensional value of electric potential
	psi = -m_Rgas*m_Tabs/m_Fc*log(zeta);
	
	return psi;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEMultiphasic::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// effective concentration
	double c = spt.m_c[sol];
	// solubility
	double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
	int z = m_pSolute[sol]->ChargeNumber();
	double zeta = ElectricPotential(pt, true);
	double zz = pow(zeta, z);
	double kappa = zz*khat;
	
	// actual concentration
	double ca = kappa*c;
	
	return ca;
}

//-----------------------------------------------------------------------------
//! The stress of a triphasic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEMultiphasic::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = m_pSolid->Stress(mp);
	
	// fluid pressure
	double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEMultiphasic::Tangent(FEMaterialPoint& mp)
{
	int i;
	
	FEElasticMaterialPoint& ept = *mp.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	const int nsol = m_pSolute.size();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	
	// relative volume and solid volume fraction
	double J = ept.J;
	double phi0 = ppt.m_phi0;
	
	// get the charge density and its derivatives
	double cF = FixedChargeDensity(mp);
	double dcFdJ = -cF/(J - phi0);
	
	// fluid pressure and solute concentration
	double p = Pressure(mp);
	
	// get remaining variables
	double zeta = ElectricPotential(mp, true);
	vector<double> c(nsol);
	vector<int> z(nsol);
	vector<double> khat(nsol);
	vector<double> dkhdJ(nsol);
	vector<double> zz(nsol);
	vector<double> kappa(nsol);
	double den = 0;
	for (i=0; i<nsol; ++i) {
		c[i] = spt.m_c[i];
		z[i] = m_pSolute[i]->ChargeNumber();
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(mp);
		dkhdJ[i] = m_pSolute[i]->m_pSolub->Tangent_Solubility_Strain(mp);
		zz[i] = pow(zeta, z[i]);
		kappa[i] = zz[i]*khat[i];
		den += SQR(z[i])*kappa[i]*c[i];
	}
	
	// evaluate electric potential (nondimensional exponential form) and its derivatives
	// also evaluate partition coefficients and their derivatives
	double zidzdJ = 0;
	if (den > 0) {
		zidzdJ = dcFdJ;
		for (i=0; i<nsol; ++i)
			zidzdJ += z[i]*zz[i]*dkhdJ[i]*c[i];
		zidzdJ = -zidzdJ/den;
	}
	vector<double> dkdJ(nsol);
	for (i=0; i<nsol; ++i) dkdJ[i] = zz[i]*dkhdJ[i]+z[i]*kappa[i]*zidzdJ;
	
	// osmotic coefficient and its derivative w.r.t. strain
	double osmc = m_pOsmC->OsmoticCoefficient(mp);
	double dodJ = m_pOsmC->Tangent_OsmoticCoefficient_Strain(mp);
	
	double dp = 0;
	for (i=0; i<nsol; ++i) {
		dp += c[i]*(osmc*dkdJ[i]+dodJ*kappa[i]);
	}
	dp *= m_Rgas*m_Tabs*J;
	
	// adjust tangent for pressures
	double D[6][6] = {0};
	C.extract(D);
	
	D[0][0] -= -p + dp;
	D[1][1] -= -p + dp;
	D[2][2] -= -p + dp;
	
	D[0][1] -= p + dp; D[1][0] -= p + dp;
	D[1][2] -= p + dp; D[2][1] -= p + dp;
	D[0][2] -= p + dp; D[2][0] -= p + dp;
	
	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
	
	return tens4ds(D);
}

//-----------------------------------------------------------------------------
//! Calculate fluid flux

vec3d FEMultiphasic::FluidFlux(FEMaterialPoint& pt)
{
	int i;
	
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = m_pSolute.size();
	vector<double> c(nsol);
	vector<vec3d> gradc(nsol);
	vector<mat3ds> D(nsol);
	vector<double> D0(nsol);
	vector<double> khat(nsol);
	vector<int> z(nsol);
	vector<double> zz(nsol);
	vector<double> kappa(nsol);
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// electric potential
	double zeta = ElectricPotential(pt, true);
	
	
	for (i=0; i<nsol; ++i) {
		// concentration
		c[i] = spt.m_c[i];

		// concentration gradient
		gradc[i] = spt.m_gradc[i];
		
		// solute diffusivity in mixture
		D[i] = m_pSolute[i]->m_pDiff->Diffusivity(pt);
		
		// solute free diffusivity
		D0[i] = m_pSolute[i]->m_pDiff->Free_Diffusivity(pt);
		
		// solubility
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
		z[i] = m_pSolute[i]->ChargeNumber();
		zz[i] = pow(zeta, z[i]);
		kappa[i] = zz[i]*khat[i];
	}
	
	// identity matrix
	mat3dd I(1);
	
	// hydraulic permeability
	mat3ds kt = m_pPerm->Permeability(pt);
	
	// effective hydraulic permeability
	mat3ds ke;
	ke.zero();
	for (i=0; i<nsol; ++i)
		ke += (I-D[i]/D0[i])*(kappa[i]*c[i]/D0[i]);
	ke = (kt.inverse() + ke*(m_Rgas*m_Tabs/phiw)).inverse();
	
	// fluid flux w
	vec3d w(0);
	for (i=0; i<nsol; ++i)
		w += (D[i]*gradc[i])*(kappa[i]/D0[i]);
	w = -(ke*(gradp + w*m_Rgas*m_Tabs));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEMultiphasic::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// concentration
	double c = spt.m_c[sol];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[sol];
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute[sol]->m_pDiff->Diffusivity(pt);
	
	// solute free diffusivity
	double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
	
	// solubility
	double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
	int z = m_pSolute[sol]->ChargeNumber();
	double zeta = ElectricPotential(pt, true);
	double zz = pow(zeta, z);
	double kappa = zz*khat;
	
	// fluid flux w
	vec3d w = FluidFlux(pt);
	
	// solute flux j
	vec3d j = (D*(w*(c/D0) - gradc*phiw))*kappa;
	
	return j;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEMultiphasic::Pressure(FEMaterialPoint& pt)
{
	int i;
	
	FEBiphasicMaterialPoint& ppt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	const int nsol = m_pSolute.size();
	
	// effective pressure
	double p = ppt.m_p;
	
	// effective concentration
	vector<double> ca(nsol);
	for (i=0; i<nsol; ++i)
		ca[i] = Concentration(pt, i);
	
	// osmotic coefficient
	double osmc = m_pOsmC->OsmoticCoefficient(pt);
	
	// actual pressure
	double pa = 0;
	for (i=0; i<nsol; ++i) pa += ca[i];
	pa = p + m_Rgas*m_Tabs*osmc*pa;
	
	return pa;
}

//-----------------------------------------------------------------------------
//! Current density
vec3d FEMultiphasic::CurrentDensity(FEMaterialPoint& pt)
{
	int i;
	const int nsol = m_pSolute.size();
	
	vector<vec3d> j(nsol);
	vector<int> z(nsol);
	vec3d Ie(0);
	for (i=0; i<nsol; ++i) {
		j[i] = SoluteFlux(pt, i);
		z[i] = m_pSolute[i]->ChargeNumber();
		Ie += j[i]*z[i];
	}
	Ie *= m_Fc;
	
	return Ie;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FEMultiphasic::Serialize(DumpFile& ar)
{
	int i, nsol;
	
	FEMaterial::Serialize(ar);
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	if (ar.IsSaving())
	{
		ar << m_phi0 << m_rhoTw << m_cFr << m_Rgas << m_Tabs << m_Fc << m_pSolute.size();
		
		ar << febio.GetTypeStr<FEMaterial>(m_pSolid); m_pSolid->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pPerm ); m_pPerm ->Serialize(ar);
		ar << febio.GetTypeStr<FEMaterial>(m_pOsmC ); m_pOsmC ->Serialize(ar);
		for (i=0; i<(int)m_pSolute.size(); ++i) {
			ar << febio.GetTypeStr<FEMaterial>(m_pSolute[i]);
			m_pSolute[i] ->Serialize(ar);
		}
	}
	else
	{
		ar >> m_phi0 >> m_rhoTw >> m_cFr >> m_Rgas >> m_Tabs >> m_Fc >> nsol ;
		
		char sz[256] = {0};
		ar >> sz;
		m_pSolid = dynamic_cast<FEElasticMaterial*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pSolid); m_pSolid->Serialize(ar);
		m_pSolid->Init();
		
		ar >> sz;
		m_pPerm = dynamic_cast<FEHydraulicPermeability*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pPerm); m_pPerm->Serialize(ar);
		m_pPerm->Init();
		
		ar >> sz;
		m_pOsmC = dynamic_cast<FEOsmoticCoefficient*>(febio.Create<FEMaterial>(sz, ar.GetFEModel()));
		assert(m_pOsmC); m_pOsmC->Serialize(ar);
		m_pOsmC->Init();
		
		for (i=0; i<nsol; ++i) {
			ar >> sz;
			m_pSolute.push_back(dynamic_cast<FESolute*>(febio.Create<FEMaterial>(sz, ar.GetFEModel())));
			assert(m_pSolute[i]); m_pSolute[i]->Serialize(ar);
			m_pSolute[i]->Init();
		}
		
	}
}
