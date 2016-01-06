// FEMultiphasic.cpp: implementation of the FEMultiphasic class.
//
//////////////////////////////////////////////////////////////////////

#include "FEMultiphasic.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <complex>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEMultiphasic, FEMaterial)
	ADD_PARAMETER2(m_phi0   , FE_PARAM_DOUBLE, FE_RANGE_CLOSED     (0.0, 1.0), "phi0"         );
	ADD_PARAMETER2(m_rhoTw  , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "fluid_density");
	ADD_PARAMETER2(m_penalty, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "penalty"      );
	ADD_PARAMETER(m_cFr     , FE_PARAM_DOUBLE, "fixed_charge_density");
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

// Newton's method for finding nearest root of a polynomial
bool newton(double& zero, const int n, const int maxit, 
			const double ep1, const double ep2, vector<double> a)
{
	bool done = false;
	bool conv = false;
	int it = 0;
	double f, df, x, dx, xi;
	x = zero;
	
	while (!done) {
		// Evaluate function and its derivative
		xi = x;
		f = a[0] + a[1]*xi;
		df = a[1];
		for (int i=2; i<=n; ++i) {
			df += i*a[i]*xi;
			xi *= x;
			f += a[i]*xi;
		}
		if (df == 0) break;
		// check absolute convergence and don't update x if met
		if (abs(f) < ep2) {
			done = true;
			conv = true;
			zero = x;
			break;
		}
		// evaluate increment in x
		dx = -f/df;
		x += dx;
		++it;
		// check relative convergence
		if (abs(dx) < ep1*abs(x)) {
			done = true;
			conv = true;
			zero = x;
		}
		// check iteration count
		else if (it > maxit) {
			done = true;
			zero = x;
		}
	}
	return conv;
}

// linear
bool poly1(vector<double> a, double& x)
{
	if (a[1]) {
		x = -a[0]/a[1];
		return true;
	} else {
		return false;
	}
}

// quadratic
bool poly2(vector<double> a, double& x)
{
	if (a[2]) {
		x = (-a[1]+sqrt(SQR(a[1])-4*a[0]*a[2]))/(2*a[2]);
		return true;
	} else {
		return poly1(a,x);
	}
}

// higher order
bool polyn(int n, vector<double> a, double& x)
{
//	bool fnreal = true;
//	vector< complex<double> > zeros(n,complex<double>(1,0));
	int maxit = 100;
	double ep1 = 1e-6;
	double ep2 = 1e-12;
	
/*	muller(fnreal, &zeros[0], n, 0, maxit, ep1, ep2, a);
	for (int i=0; i<n; ++i) {
		if (zeros[i].real() > 0) {
			x = zeros[i].real();
			return true;
		}
	}*/
	return newton(x, n, maxit,ep1, ep2, a);
}

bool solvepoly(int n, vector<double> a, double& x)
{
	switch (n) {
		case 1:
			return poly1(a, x);
			break;
		case 2:
			return poly2(a, x);
		default:
			if (a[n]) {
				return polyn(n, a, x);
			} else {
				return solvepoly(n-1, a, x);
			}
			break;
	}
}

//=============================================================================
//   FEMultiphasic
//=============================================================================

//-----------------------------------------------------------------------------
//! FEMultiphasic constructor
FEMultiphasic::FEMultiphasic(FEModel* pfem) : FEMaterial(pfem)
{	
	m_phi0 = 0;
	m_rhoTw = 0;
	m_cFr = 0;
	m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
	m_penalty = 1;

	// define the material properties
	AddProperty(&m_pSolid , "solid"              );
	AddProperty(&m_pPerm  , "permeability"       );
	AddProperty(&m_pOsmC  , "osmotic_coefficient");
	AddProperty(&m_pSupp  , "solvent_supply"     , false);
	AddProperty(&m_pSolute, "solute"             , false);
	AddProperty(&m_pSBM   , "solid_bound"        , false);
	AddProperty(&m_pReact , "reaction"           , false);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddSolute(FESolute* psol)
{
	m_pSolute.SetProperty(psol);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddSolidBoundMolecule(FESolidBoundMolecule* psbm)
{
	m_pSBM.SetProperty(psbm);
}

//-----------------------------------------------------------------------------
void FEMultiphasic::AddChemicalReaction(FEChemicalReaction* pcr)
{
	m_pReact.SetProperty(pcr);
}

//-----------------------------------------------------------------------------
//! Returns the local ID of the SBM, given the global ID.
//! \param nid global ID (one - based)
//! \return the local ID (zero-based index) or -1 if not found.
int FEMultiphasic::FindLocalSBMID(int nid)
{
	int lsbm = -1;
	int nsbm = (int) SBMs();
	for (int isbm=0; isbm<nsbm; ++isbm) {
		if (m_pSBM[isbm]->GetSBMID() == nid - 1) {
			lsbm = isbm;
			break;
		}
	}
	return lsbm;
}

//-----------------------------------------------------------------------------
//! Returns the local solute index given the global ID
int FEMultiphasic::FindLocalSoluteID(int nid)
{
	int lsid = -1;
	for (int isol=0; isol<Solutes(); ++isol)
		if (GetSolute(isol)->GetSoluteID() == nid) {
                lsid = isol;
                break;
		}
	return lsid;
}

//-----------------------------------------------------------------------------
bool FEMultiphasic::InitializeReaction(FEChemicalReaction* m_pReact)
{
	int isol, isbm, itot;
	
	const int nsol = (int)m_pSolute.size();
	const int nsbm = (int)m_pSBM.size();
	const int ntot = nsol + nsbm;
	
	// initialize the stoichiometric coefficients to zero
	m_pReact->m_nsol = nsol;
	m_pReact->m_vR.assign(ntot,0);
	m_pReact->m_vP.assign(ntot,0);
	m_pReact->m_v.assign(ntot,0);
	
	// cycle through all the solutes in the mixture and determine
	// if they participate in this reaction
	itrmap it;
	intmap solR = m_pReact->m_solR;
	intmap solP = m_pReact->m_solP;
	for (isol=0; isol<nsol; ++isol) {
		int sid = m_pSolute[isol]->GetSoluteID();
		it = solR.find(sid);
		if (it != solR.end()) m_pReact->m_vR[isol] = it->second;
		it = solP.find(sid);
		if (it != solP.end()) m_pReact->m_vP[isol] = it->second;
	}
	
	// cycle through all the solid-bound molecules in the mixture
	// and determine if they participate in this reaction
	intmap sbmR = m_pReact->m_sbmR;
	intmap sbmP = m_pReact->m_sbmP;
	for (isbm=0; isbm<nsbm; ++isbm) {
		int sid = m_pSBM[isbm]->GetSBMID();
		it = sbmR.find(sid);
		if (it != sbmR.end()) m_pReact->m_vR[nsol+isbm] = it->second;
		it = sbmP.find(sid);
		if (it != sbmP.end()) m_pReact->m_vP[nsol+isbm] = it->second;
	}
	
	// evaluate the net stoichiometric coefficient
	for (itot=0; itot<ntot; ++itot) {
		m_pReact->m_v[itot] = m_pReact->m_vP[itot] - m_pReact->m_vR[itot];
	}
    
    // evaluate the weighted molar volume of reactants and products
    if (!m_pReact->m_Vovr) {
        m_pReact->m_Vbar = 0;
        for (isol=0; isol<nsol; ++isol)
            m_pReact->m_Vbar += m_pReact->m_v[isol]*m_pSolute[isol]->MolarMass()/m_pSolute[isol]->Density();
        for (isbm=0; isbm<nsbm; ++isbm)
            m_pReact->m_Vbar += m_pReact->m_v[nsol+isbm]*m_pSBM[isbm]->MolarMass()/m_pSBM[isbm]->Density();
    }
	
	// check that the chemical reaction satisfies electroneutrality
	int znet = 0;
	for (isol=0; isol<nsol; ++isol)
		znet += m_pReact->m_v[isol]*m_pSolute[isol]->ChargeNumber();
	for (isbm=0; isbm<nsbm; ++isbm)
		znet += m_pReact->m_v[nsol+isbm]*m_pSBM[isbm]->ChargeNumber();
	if (znet != 0) return MaterialError("chemical reaction must satisfy electroneutrality");
	
	// set pointer to this multiphasic material
	m_pReact->m_pMP = this;

	return true;
}

//-----------------------------------------------------------------------------
bool FEMultiphasic::Init()
{
	// we first have to set the parent material
	// TODO: This seems redundant since each material already has a pointer to its parent
	for (int i=0; i<Reactions(); ++i)
	{
		m_pReact[i]->m_pMP = this;
	}

	// call the base class.
	// This also initializes all properties
	if (FEMaterial::Init() == false) return false;

	// set the solute IDs first, since they are referenced in FESolute::Init()
	for (int i=0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }

	// initialize chemical reactions
	for (int i=0; i<Reactions(); ++i)
	{
		if (InitializeReaction(m_pReact[i]) == false) return false;
	}

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

	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	m_Fc   = GetFEModel()->GetGlobalConstant("Fc");
	
	if (m_Rgas <= 0) return MaterialError("A positive universal gas constant R must be defined in Globals section");
	if (m_Tabs <= 0) return MaterialError("A positive absolute temperature T must be defined in Globals section");
	if ((zmin || zmax) && (m_Fc <= 0)) return MaterialError("A positive Faraday constant Fc must be defined in Globals section");

	return true;
}

//-----------------------------------------------------------------------------
//! Solid referential apparent density
double FEMultiphasic::SolidReferentialApparentDensity(FEMaterialPoint& pt)
{
	FEBiphasicMaterialPoint& pet = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
		
	// evaluate referential apparent density of base solid
	double rhosr = pet.m_phi0*m_pSolid->Density();

	// add contribution from solid-bound molecules
	for (int isbm=0; isbm<(int)spt.m_sbmr.size(); ++isbm)
		rhosr += spt.m_sbmr[isbm];
	
	return rhosr;
}

//-----------------------------------------------------------------------------
//! Solid referential volume fraction
double FEMultiphasic::SolidReferentialVolumeFraction(FEMaterialPoint& pt)
{
	// get referential apparent density of base solid (assumed constant)
	double phisr = m_phi0;
    
	// add contribution from solid-bound molecules
	for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm)
		phisr += SBMReferentialVolumeFraction(pt, isbm);
	
	return phisr;
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEMultiphasic::Porosity(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	// solid referential volume fraction
	double phisr = bt.m_phi0;
    
	// relative volume
	double J = et.m_J;
	
	double phiw = 1 - phisr/J;
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
	FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// relative volume
	double J = et.m_J;
	double phi0 = bt.m_phi0;
	double ce = 0;

	// add contribution from charged solid-bound molecules
	for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm)
		ce += SBMChargeNumber(isbm)*spt.m_sbmr[isbm]/SBMMolarMass(isbm);
	
	double cF = (m_cFr*(1-phi0)+ce)/(J-phi0);

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
	const int nsol = (int)m_pSolute.size();
	double cF = FixedChargeDensity(pt);

	vector<double> c(nsol);		// effective concentration
	vector<double> khat(nsol);	// solubility
	vector<int> z(nsol);		// charge number
	for (i=0; i<nsol; ++i) {
		c[i] = set.m_c[i];
		khat[i] = m_pSolute[i]->m_pSolub->Solubility(pt);
		z[i] = m_pSolute[i]->ChargeNumber();
	}
	
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
	double psi = set.m_psi;		// use previous solution as initial guess
	double zeta = exp(-m_Fc*psi/m_Rgas/m_Tabs);
	if (!solvepoly(n, a, zeta)) {
		zeta = 1.0;
	}
	
	// Return exponential (non-dimensional) form if desired
	if (eform) return zeta;
	
	// Otherwise return dimensional value of electric potential
	psi = -m_Rgas*m_Tabs/m_Fc*log(zeta);
	
	return psi;
}

//-----------------------------------------------------------------------------
//! partition coefficient
double FEMultiphasic::PartitionCoefficient(FEMaterialPoint& pt, const int sol)
{
	
	// solubility
	double khat = m_pSolute[sol]->m_pSolub->Solubility(pt);
	// charge number
	int z = m_pSolute[sol]->ChargeNumber();
	// electric potential
	double zeta = ElectricPotential(pt, true);
	double zz = pow(zeta, z);
	// partition coefficient
	double kappa = zz*khat;
	
	return kappa;
}

//-----------------------------------------------------------------------------
//! partition coefficients and their derivatives
void FEMultiphasic::PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                                  vector<double>& dkdJ,
                                                  vector< vector<double> >& dkdc,
                                                  vector< vector<double> >& dkdr,
                                                  vector< vector<double> >& dkdJr,
                                                  vector< vector< vector<double> > >& dkdrc)
{
	int isol, jsol, ksol;
    int isbm;
	
	FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
	FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
	
	const int nsol = (int)m_pSolute.size();
    const int nsbm = (int)m_pSBM.size();
	
	vector<double> c(nsol);
	vector<int> z(nsol);
	vector<double> khat(nsol);
	vector<double> dkhdJ(nsol);
	vector<double> dkhdJJ(nsol);
	vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
	vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
	vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);	// use dkhdc to initialize only
	vector<double> zz(nsol);
	kappa.resize(nsol);

	double den = 0;
    double num = 0;
	double zeta = ElectricPotential(mp, true);

	for (isol=0; isol<nsol; ++isol) {
		// get the effective concentration, its gradient and its time derivative
		c[isol] = spt.m_c[isol];
		// get the charge number
		z[isol] = m_pSolute[isol]->ChargeNumber();
		// evaluate the solubility and its derivatives w.r.t. J and c
		khat[isol] = m_pSolute[isol]->m_pSolub->Solubility(mp);
		dkhdJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(mp);
		dkhdJJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Strain(mp);
		for (jsol=0; jsol<nsol; ++jsol) {
			dkhdc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(mp,jsol);
			dkhdJc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Concentration(mp,jsol);
			for (ksol=0; ksol<nsol; ++ksol) {
				dkhdcc[isol][jsol][ksol] = 
				m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration_Concentration(mp,jsol,ksol);
			}
		}
		zz[isol] = pow(zeta, z[isol]);
		kappa[isol] = zz[isol]*khat[isol];
		den += SQR(z[isol])*kappa[isol]*c[isol];
        num += pow((double)z[isol],3)*kappa[isol]*c[isol];
	}
	
	// get the charge density and its derivatives
	double J = ept.m_J;
	double phi0 = ppt.m_phi0;
	double cF = FixedChargeDensity(mp);
	double dcFdJ = -cF/(J - phi0);
	double dcFdJJ = 2*cF/SQR(J-phi0);
	
	// evaluate electric potential (nondimensional exponential form) and its derivatives
	// also evaluate partition coefficients and their derivatives
	double zidzdJ = 0;
	double zidzdJJ = 0, zidzdJJ1 = 0, zidzdJJ2 = 0;
	vector<double> zidzdc(nsol,0);
	vector<double> zidzdJc(nsol,0), zidzdJc1(nsol,0), zidzdJc2(nsol,0);
	vector< vector<double> > zidzdcc(nsol, vector<double>(nsol,0));
	vector< vector<double> > zidzdcc1(nsol, vector<double>(nsol,0));
	vector<double> zidzdcc2(nsol,0);
	double zidzdcc3 = 0;

	if (den > 0) {
		
		for (isol=0; isol<nsol; ++isol)
			zidzdJ += z[isol]*zz[isol]*dkhdJ[isol]*c[isol];
		zidzdJ = -(dcFdJ+zidzdJ)/den;
		
		for (isol=0; isol<nsol; ++isol) {
			for (jsol=0; jsol<nsol; ++jsol) {
				zidzdJJ1 += SQR(z[jsol])*c[jsol]*(z[jsol]*zidzdJ*kappa[jsol]+zz[jsol]*dkhdJ[jsol]);
				zidzdJJ2 += z[jsol]*zz[jsol]*c[jsol]*(zidzdJ*z[jsol]*dkhdJ[jsol]+dkhdJJ[jsol]);
				zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
			}
			zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
			zidzdcc3 += pow(double(z[isol]),3)*kappa[isol]*c[isol];
		}
		zidzdJJ = zidzdJ*(zidzdJ-zidzdJJ1/den)-(dcFdJJ+zidzdJJ2)/den;
		
		for (isol=0; isol<nsol; ++isol) {
			for (jsol=0; jsol<nsol; ++jsol) {
				zidzdJc1[isol] += SQR(z[jsol])*c[jsol]*(zidzdc[isol]*z[jsol]*kappa[jsol]+zz[jsol]*dkhdc[jsol][isol]);
				zidzdJc2[isol] += z[jsol]*zz[jsol]*c[jsol]*(zidzdc[isol]*z[jsol]*dkhdJ[jsol]+dkhdJc[jsol][isol]);
				zidzdcc2[isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
				for (ksol=0; ksol<nsol; ++ksol)
					zidzdcc1[isol][jsol] += z[ksol]*zz[ksol]*c[ksol]*dkhdcc[ksol][isol][jsol];
			}
			zidzdJc[isol] = zidzdJ*(zidzdc[isol]-(SQR(z[isol])*kappa[isol] + zidzdJc1[isol])/den)
			-(z[isol]*zz[isol]*dkhdJ[isol] + zidzdJc2[isol])/den;
		}
		
		for (isol=0; isol<nsol; ++isol) {
			for (jsol=0; jsol<nsol; ++jsol) {
				zidzdcc[isol][jsol] = zidzdc[isol]*zidzdc[jsol]*(1 - zidzdcc3/den)
				- zidzdcc1[isol][jsol]/den
				- z[isol]*(z[isol]*kappa[isol]*zidzdc[jsol]+zz[isol]*dkhdc[isol][jsol])/den
				- z[jsol]*(z[jsol]*kappa[jsol]*zidzdc[isol]+zz[jsol]*dkhdc[jsol][isol])/den
				- zidzdc[jsol]*zidzdcc2[isol]/den
				- zidzdc[isol]*zidzdcc2[jsol]/den;
			}
		}
	}
	
	dkdJ.resize(nsol);
	dkdc.resize(nsol, vector<double>(nsol,0));
	
	for (isol=0; isol<nsol; ++isol) {
		dkdJ[isol] = zz[isol]*dkhdJ[isol]+z[isol]*kappa[isol]*zidzdJ;
		for (jsol=0; jsol<nsol; ++jsol) {
			dkdc[isol][jsol] = zz[isol]*dkhdc[isol][jsol]+z[isol]*kappa[isol]*zidzdc[jsol];
		}
	}
    vector<double> zidzdr(nsbm,0);
    vector<double> zidzdJr(nsbm,0);
	vector< vector<double> > zidzdrc(nsbm, vector<double>(nsol,0));
	dkdr.resize(nsol, vector<double>(nsbm));
	dkdJr.resize(nsol, vector<double>(nsbm));
	dkdrc.resize(nsol, zidzdrc);	// use zidzdrc for initialization only
    
	if (den > 0) {
		
        for (isbm=0; isbm<nsbm; ++isbm) {
            zidzdr[isbm] = -(cF/SBMDensity(isbm) + SBMChargeNumber(isbm)/SBMMolarMass(isbm))/(J-phi0)/den;
            
            for (isol=0; isol<nsol; ++isol) {
                zidzdJr[isbm] += SQR(z[isol])*dkdJ[isol]*c[isol];
            }
            zidzdJr[isbm] = 1/(J-phi0) + zidzdJr[isbm]/den;
            zidzdJr[isbm] = (zidzdJr[isbm] + zidzdJ)*zidzdr[isbm];
            zidzdJr[isbm] += cF/SBMDensity(isbm)/SQR(J-phi0)/den;
            
            for (isol=0; isol<nsol; ++isol) {
                zidzdrc[isbm][isol] = SQR(z[isol])*kappa[isol];
                for (jsol=0; jsol<nsol; ++jsol)
                    zidzdrc[isbm][isol] += SQR(z[jsol])*zz[jsol]*c[jsol]*dkhdc[jsol][isol];
                zidzdrc[isbm][isol] = zidzdr[isbm]*(zidzdc[isol]*(1+num/den) - zidzdrc[isbm][isol]/den);
            }
        }
	}
    
    for (isbm=0; isbm<nsbm; ++isbm) {
        for (isol=0; isol<nsol; ++isol) {
            dkdr[isol][isbm] = z[isol]*kappa[isol]*zidzdr[isbm];
            dkdJr[isol][isbm] = z[isol]*((dkhdJ[isol]*zz[isol]+(z[isol]-1)*kappa[isol]*zidzdJ)*zidzdr[isbm]
                                         +kappa[isol]*zidzdJr[isbm]);
            for (jsol=0; jsol<nsol; ++jsol) {
                dkdrc[isol][isbm][jsol] = z[isol]*(zz[isol]*dkhdc[isol][jsol]*zidzdr[isbm]
                                                   +(z[isol]-1)*kappa[isol]*zidzdr[isbm]*zidzdc[jsol]
                                                   +kappa[isol]*zidzdrc[isbm][jsol]);
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEMultiphasic::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// effective concentration
	double c = spt.m_c[sol];
	
	// partition coefficient
	double kappa = PartitionCoefficient(pt, sol);
	
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
	FEBiphasicMaterialPoint& bpt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	const int nsol = (int)m_pSolute.size();
	
	// call solid tangent routine
	tens4ds C = m_pSolid->Tangent(mp);
	double D[6][6] = {0};
	C.extract(D);
	
	// relative volume and solid volume fraction
	double J = ept.m_J;
	double phi0 = bpt.m_phi0;
	
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
	const int nsol = (int)m_pSolute.size();
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
	vec3d w(0,0,0);
	for (i=0; i<nsol; ++i)
		w += (D[i]*gradc[i])*(kappa[i]/D0[i]);
	w = -(ke*(gradp + w*m_Rgas*m_Tabs));
	
	return w;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEMultiphasic::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
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
	vec3d w = bpt.m_w;
	
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
	const int nsol = (int)m_pSolute.size();
	
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
	const int nsol = (int)m_pSolute.size();
	
	vector<vec3d> j(nsol);
	vector<int> z(nsol);
	vec3d Ie(0,0,0);
	for (i=0; i<nsol; ++i) {
		j[i] = SoluteFlux(pt, i);
		z[i] = m_pSolute[i]->ChargeNumber();
		Ie += j[i]*z[i];
	}
	Ie *= m_Fc;
	
	return Ie;
}
