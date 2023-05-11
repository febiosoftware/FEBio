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
#include "FEFluidSolutes.h"
#include "FECore/FEModel.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>
#include <complex>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidSolutes, FEMaterial)
// material properties
// ADD_PARAMETER(m_penalty, FE_RANGE_GREATER_OR_EQUAL(0.0), "penalty"      );
ADD_PARAMETER(m_diffMtmSupp , "dms")->setLongName("include osmosis");
ADD_PROPERTY(m_pFluid, "fluid");
ADD_PROPERTY(m_pOsmC  , "osmotic_coefficient");
ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Polynomial root solver

// function whose roots needs to be evaluated
void fnFS(complex<double>& z, complex<double>& fz, vector<double> a)
{
    int n = (int)a.size()-1;
    fz = a[0];
    complex<double> x(1,0);
    
    for (int i=1; i<=n; ++i) {
        x *= z;
        fz += a[i]*x;
    }
    return;
}

// deflation
bool dflateFS(complex<double> zero, const int i, int& kount,
            complex<double>& fzero, complex<double>& fzrdfl,
            complex<double>* zeros, vector<double> a)
{
    complex<double> den;
    ++kount;
    fnFS(zero, fzero, a);
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
void mullerFS(bool fnreal, complex<double>* zeros, const int n, const int nprev,
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
        if (dflateFS(z, i, kount, fzr, dvdf1p, zeros, a)) goto eloop;
        z = zero - 0.5;
        if (dflateFS(z, i, kount, fzr, fzrprv, zeros, a)) goto eloop;
        dvdf1p = (fzrprv - dvdf1p)/hprev;
        if (dflateFS(zero, i, kount, fzr, fzrdfl, zeros, a)) goto eloop;
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
            fnFS(zero,fzrdfl,a);
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
bool newtonFS(double& zero, const int n, const int maxit,
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
bool poly1FS(vector<double> a, double& x)
{
    if (a[1]) {
        x = -a[0]/a[1];
        return true;
    } else {
        return false;
    }
}

// quadratic
bool poly2FS(vector<double> a, double& x)
{
    if (a[2]) {
        x = (-a[1]+sqrt(SQR(a[1])-4*a[0]*a[2]))/(2*a[2]);
        return true;
    } else {
        return poly1FS(a,x);
    }
}

// higher order
bool polynFS(int n, vector<double> a, double& x)
{
    //    bool fnreal = true;
    //    vector< complex<double> > zeros(n,complex<double>(1,0));
    int maxit = 100;
    double ep1 = 1e-6;
    double ep2 = 1e-12;
    
    /*    mullerFS(fnreal, &zeros[0], n, 0, maxit, ep1, ep2, a);
     for (int i=0; i<n; ++i) {
     if (zeros[i].real() > 0) {
     x = zeros[i].real();
     return true;
     }
     }*/
    return newtonFS(x, n, maxit,ep1, ep2, a);
}

bool solvepolyFS(int n, vector<double> a, double& x)
{
    switch (n) {
        case 1:
            return poly1FS(a, x);
            break;
        case 2:
            return poly2FS(a, x);
        default:
            if (a[n]) {
                return polynFS(n, a, x);
            } else {
                return solvepolyFS(n-1, a, x);
            }
            break;
    }
}

//============================================================================
// FEFluidSolutesMaterialPoint
//============================================================================
FEFluidSolutesMaterialPoint::FEFluidSolutesMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEFluidSolutesMaterialPoint::Copy()
{
    FEFluidSolutesMaterialPoint* pt = new FEFluidSolutesMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_nsol & m_psi & m_Ie & m_pe;
    ar & m_c & m_ca & m_gradc & m_j & m_cdot & m_k & m_dkdJ;
    ar & m_dkdc;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Init()
{
    m_nsol = 0;
    m_psi = 0;
    m_pe = 0;
    m_Ie = vec3d(0,0,0);
    m_c.clear();
    m_ca.clear();
    m_gradc.clear();
    m_j.clear();
    m_cdot.clear();
    m_k.clear();
    m_dkdJ.clear();
    m_dkdJJ.clear();
    m_dkdc.clear();
    m_dkdJc.clear();
    m_dkdcc.clear();

	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
double FEFluidSolutesMaterialPoint::Osmolarity() const
{
    double ew = 0.0;
    for (int isol = 0; isol < (int)m_ca.size(); ++isol)
    {
        ew += m_ca[isol];
    }
    return ew;
}

//============================================================================
// FEFluidSolutes
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidSolutes constructor

FEFluidSolutes::FEFluidSolutes(FEModel* pfem) : FEMaterial(pfem)
{
    m_pFluid = 0;
    m_Rgas = 0; m_Tabs = 0; m_Fc = 0;
    m_diffMtmSupp = false;
    m_penalty = 1;
    m_pOsmC = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEFluidSolutes::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    return new FEFluidSolutesMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
void FEFluidSolutes::AddChemicalReaction(FEChemicalReaction* pcr)
{
    m_pReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidSolutes::Init()
{
    // set the solute IDs first, since they are referenced in FESolute::Init()
    for (int i = 0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }
    
    // call the base class.
    // This also initializes all properties
    if (FEMaterial::Init() == false) return false;
    
    // Determine how to solve for the electric potential psi
    int isol;
    int zmin = 0, zmax = 0, z;
    for (isol=0; isol<(int)m_pSolute.size(); ++isol) {
        z = m_pSolute[isol]->ChargeNumber();
        if (z < zmin) zmin = z;
        if (z > zmax) zmax = z;
    }
    m_zmin = zmin;
    m_ndeg = zmax - zmin;    // polynomial degree
    
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    m_Fc   = GetFEModel()->GetGlobalConstant("Fc");
    
    if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
    if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");     return false; }
    if ((zmin || zmax) && (m_Fc <= 0)) {
        feLogError("A positive Faraday constant Fc must be defined in Globals section");
        return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutes::Serialize(DumpStream& ar)
{
    FEMaterial::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_Rgas & m_Tabs & m_Fc;
    ar & m_zmin & m_ndeg;
}

//-----------------------------------------------------------------------------
//! Electric potential
double FEFluidSolutes::ElectricPotential(const FEMaterialPoint& pt, const bool eform)
{
    // check if solution is neutral
    if (m_ndeg == 0) {
        if (eform) return 1.0;
        else return 0.0;
    }
    
    int i, j;
    
    // if not neutral, solve electroneutrality polynomial for zeta
    const FEFluidSolutesMaterialPoint& set = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    double cF = 0.0;
    
    FEMaterialPoint mp = pt;
    vector<double> c(nsol);        // effective concentration
    vector<double> khat(nsol);    // solubility
    vector<int> z(nsol);        // charge number
    for (i=0; i<nsol; ++i) {
        c[i] = set.m_c[i];
        khat[i] = m_pSolute[i]->m_pSolub->Solubility(mp);
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
    double psi = set.m_psi;        // use previous solution as initial guess
    double zeta = exp(-m_Fc*psi/m_Rgas/m_Tabs);
    if (!solvepolyFS(n, a, zeta)) {
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
double FEFluidSolutes::PartitionCoefficient(const FEMaterialPoint& pt, const int sol)
{
    FEMaterialPoint mp = pt;
    // solubility
    double khat = m_pSolute[sol]->m_pSolub->Solubility(mp);
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
void FEFluidSolutes::PartitionCoefficientFunctions(const FEMaterialPoint& mp, vector<double>& kappa,
                                                  vector<double>& dkdJ,
                                                  vector< vector<double> >& dkdc)
{
    //TODO: Include dkdcc and dkdJc
    
    int isol, jsol, ksol;
    
    const FEFluidMaterialPoint& fpt = *(mp.ExtractData<FEFluidMaterialPoint>());
    const FEFluidSolutesMaterialPoint& spt = *(mp.ExtractData<FEFluidSolutesMaterialPoint>());
    FEMaterialPoint pt = mp;
    
    const int nsol = (int)m_pSolute.size();
    
    vector<double> c(nsol);
    vector<int> z(nsol);
    vector<double> khat(nsol);
    vector<double> dkhdJ(nsol);
    vector<double> dkhdJJ(nsol);
    vector< vector<double> > dkhdc(nsol, vector<double>(nsol));
    vector< vector<double> > dkhdJc(nsol, vector<double>(nsol));
    vector< vector< vector<double> > > dkhdcc(nsol, dkhdc);    // use dkhdc to initialize only
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
        khat[isol] = m_pSolute[isol]->m_pSolub->Solubility(pt);
        dkhdJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain(pt);
        dkhdJJ[isol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Strain(pt);
        for (jsol=0; jsol<nsol; ++jsol) {
            dkhdc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration(pt,jsol);
            dkhdJc[isol][jsol] = m_pSolute[isol]->m_pSolub->Tangent_Solubility_Strain_Concentration(pt,jsol);
            for (ksol=0; ksol<nsol; ++ksol) {
                dkhdcc[isol][jsol][ksol] =
                m_pSolute[isol]->m_pSolub->Tangent_Solubility_Concentration_Concentration(pt,jsol,ksol);
            }
        }
        zz[isol] = pow(zeta, z[isol]);
        kappa[isol] = zz[isol]*khat[isol];
        den += SQR(z[isol])*kappa[isol]*c[isol];
        num += pow((double)z[isol],3)*kappa[isol]*c[isol];
    }

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
        zidzdJ = -(zidzdJ)/den;
        
        for (isol=0; isol<nsol; ++isol) {
            for (jsol=0; jsol<nsol; ++jsol) {
                zidzdJJ1 += SQR(z[jsol])*c[jsol]*(z[jsol]*zidzdJ*kappa[jsol]+zz[jsol]*dkhdJ[jsol]);
                zidzdJJ2 += z[jsol]*zz[jsol]*c[jsol]*(zidzdJ*z[jsol]*dkhdJ[jsol]+dkhdJJ[jsol]);
                zidzdc[isol] += z[jsol]*zz[jsol]*dkhdc[jsol][isol]*c[jsol];
            }
            zidzdc[isol] = -(z[isol]*kappa[isol]+zidzdc[isol])/den;
            zidzdcc3 += pow(double(z[isol]),3)*kappa[isol]*c[isol];
        }
        zidzdJJ = zidzdJ*(zidzdJ-zidzdJJ1/den)-(zidzdJJ2)/den;
        
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
}

//-----------------------------------------------------------------------------
//! Current density
vec3d FEFluidSolutes::CurrentDensity(const FEMaterialPoint& pt)
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

//-----------------------------------------------------------------------------
//! actual concentration
double FEFluidSolutes::ConcentrationActual(const FEMaterialPoint& pt, const int sol)
{
    const FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    
    // effective concentration
    double c = spt.m_c[sol];
    
    // partition coefficient
    double kappa = PartitionCoefficient(pt, sol);
    
    double ca = kappa*c;
    
    return ca;
}

//-----------------------------------------------------------------------------
//! actual fluid pressure
double FEFluidSolutes::PressureActual(const FEMaterialPoint& pt)
{
    int i;
    
    const FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    const FEFluidSolutesMaterialPoint& spt = *(pt.ExtractData<FEFluidSolutesMaterialPoint>());
    const int nsol = (int)m_pSolute.size();
    
    // effective pressure
    double p = spt.m_pe;
    
    // actual concentration
    vector<double> c(nsol);
    for (i=0; i<nsol; ++i)
        c[i] = ConcentrationActual(pt, i);
    
    // osmotic coefficient
    FEMaterialPoint mp = pt;
    double osmc = m_pOsmC->OsmoticCoefficient(mp);
    
    // actual pressure
    double pa = 0;
    for (i=0; i<nsol; ++i) pa += c[i];
    pa = p + m_Rgas*m_Tabs*osmc*pa;
    
    return pa;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEFluidSolutes::SoluteFlux(const FEMaterialPoint& pt, const int sol)
{
    const FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    const FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    FEMaterialPoint mp = pt;
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(mp);
    double kappa = PartitionCoefficient(pt, sol);
    
    double c = spt.m_c[sol];
    vec3d v = fpt.m_vft;
    
    // solute flux j
    vec3d j = -gradc*D0*kappa + v*c*kappa;
    
    return j;
}

//-----------------------------------------------------------------------------
//! Calculate diffusive solute molar flux

vec3d FEFluidSolutes::SoluteDiffusiveFlux(const FEMaterialPoint& pt, const int sol)
{
    const FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    const FEFluidMaterialPoint& fpt = *pt.ExtractData<FEFluidMaterialPoint>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    FEMaterialPoint mp = pt;
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(mp);
    double kappa = PartitionCoefficient(pt, sol);
    
    // diffusive solute flux j
    vec3d jd = -gradc*D0*kappa;
    
    return jd;
}

//-----------------------------------------------------------------------------
// solute interface functions

double FEFluidSolutes::GetEffectiveSoluteConcentration(FEMaterialPoint& mp, int soluteIndex)
{
    FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    return spt.m_c[soluteIndex];
}

double FEFluidSolutes::GetFreeDiffusivity(FEMaterialPoint& mp, int soluteIndex)
{
    return m_pSolute[soluteIndex]->m_pDiff->Free_Diffusivity(mp);
}

double FEFluidSolutes::GetPartitionCoefficient(FEMaterialPoint& mp, int soluteIndex)
{
    FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    return spt.m_k[soluteIndex];
}

double FEFluidSolutes::GetOsmolarity(const FEMaterialPoint& mp)
{
    double osm = 0;
    const FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    const int nsol = (int)m_pSolute.size();
    for (int i=0; i<nsol; ++i) osm += spt.m_ca[i];
    return osm;
}

double FEFluidSolutes::dkdc(const FEMaterialPoint& mp, int i, int j)
{
    const FEFluidSolutesMaterialPoint& spt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    return spt.m_dkdc[i][j];
}
