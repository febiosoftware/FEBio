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
#include <cmath>
#include "FEReactivePlasticReturnMap.h"
#include "FEElasticMaterial.h"
#include "FEDamageCriterion.h"

//-----------------------------------------------------------------------------
FEReactivePlasticReturnMap::FEReactivePlasticReturnMap(FEElasticMaterial* pBase, FEDamageCriterion* pCrit, bool isochoric, double rtol)
{
    m_pBase = pBase;
    m_pCrit = pCrit;
    m_isochrc = isochoric;
    m_rtol = rtol;
}

//-----------------------------------------------------------------------------
// get the tensor M such that dPhi = M : L, where L = dF*F^-1 (Eq. 3.27)
mat3d FEReactivePlasticReturnMap::YieldStressTangent(FEMaterialPoint& mp)
{
    mat3ds s = m_pBase->Stress(mp);
    tens4ds c = m_pBase->Tangent(mp);
    mat3ds dPhi = m_pCrit->CriterionStressTangent(mp);
    return dPhi*s*2 - mat3dd((dPhi*s).trace()) + c.dot(dPhi);
}

//-----------------------------------------------------------------------------
// get the yield surface normal dPhi/dU (Eq. 3.26)
mat3ds FEReactivePlasticReturnMap::YieldSurfaceNormal(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d M = YieldStressTangent(mp);
    mat3ds Ui = pe.RightStretchInverse();
    mat3d R = pe.m_F*Ui;
    mat3ds N = (R.transpose()*M*R*Ui).sym();
    return N;
}

//-----------------------------------------------------------------------------
// original formulation: Uv = xi*sym(Ue*(I - lam*Nhat)), Fv = R*Uv (Eqs. 3.20-3.22)
FEReactivePlasticReturnMap::Result FEReactivePlasticReturnMap::Original(FEMaterialPoint& pt, const mat3d& Fs, const mat3d& Fusi, double Ky)
{
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();

    // store safe copy of total deformation gradient
    mat3d Ftmp = pe.m_F;
    double Jtmp = pe.m_J;
    mat3d R = pe.m_F*pe.RightStretchInverse();

    // elastic predictor
    mat3d Fe = Fs*Fusi;
    pe.m_F = Fe; pe.m_J = Fe.det();
    mat3ds Ue = pe.RightStretch();
    double Kv = m_pCrit->DamageCriterion(pt);
    double phi = Kv - Ky;

    // find Fv
    bool conv = false;
    bool triedLroot = false;
    int maxit = 20;
    int iter = 0;
    double lam = 0;
    double tol = m_rtol*Ky;     // convergence tolerance on phi (same units as Ky)
    mat3d Fv = Fe;
    mat3ds Uv = Ue;
    mat3ds Nv = YieldSurfaceNormal(pt);
    double Nvmag = Nv.norm();
    mat3dd I(1);
    double beta = 1;
    mat3ds ImN = I;
    double phi0=0, phi1=0, phi2=0, lam0 = 0, lam1=0, lam2=0, a=0, b=0, c=0, d=0;
    bool mroot = false;
    double lroot = 0;
    double lbest = 0, phibest = phi;    // iterate closest to the yield surface
    // evaluate Fv for the current value of lam
    auto updateFv = [&]() {
        ImN = I - Nv*(lam/Nvmag);
        if (m_isochrc) beta = pow((Fusi*ImN).det(), -1./3.);
        Uv = (Ue*ImN).sym()*beta;
        Fv = R*Uv;
    };
    while (true) {
        // evaluate the actual yield function at the current Fv
        pe.m_F = Fv; pe.m_J = Fv.det();
        Kv = m_pCrit->DamageCriterion(pt);
        phi = Kv - Ky;    // phi = 0 => stay on yield surface
        if (fabs(phi) <= tol) { conv = true; break; }
        if (fabs(phi) < fabs(phibest)) { phibest = phi; lbest = lam; }
        if (iter > maxit) {
            // if a quadratic fit produced two roots of the same sign, try the larger one
            if (mroot && !triedLroot) {
                triedLroot = true;
                mroot = false;
                iter = 0;
                lam = lroot;
                updateFv();
                continue;
            }
            break;
        }
        int it = iter % 3;
        if (it == 0) {
            phi0 = phi;
            lam0 = lam;
        }
        else if (it == 1) {
            phi1 = phi;
            lam1 = lam;
        }
        else if (it == 2) {
            phi2 = phi;
            lam2 = lam;
        }
        mat3d dUvdlam = -Ue*Nv*(beta/Nvmag);
        if (m_isochrc)
            dUvdlam += Ue*ImN*((ImN.inverse()*Nv/Nvmag).trace()*beta/3.);
        double dlam = -phi/(Nv*dUvdlam.transpose()).trace();
        lam += dlam;
        // every third iteration, replace the Newton update with a quadratic fit
        // of phi(lam) through the last three samples (keep the Newton update if
        // the samples are not distinct)
        if (it == 2) {
            d = (lam0-lam1)*(lam0-lam2)*(lam1-lam2);
            if (d != 0) {
                a = (lam2*(phi1-phi0)+lam1*(phi0-phi2)+lam0*(phi2-phi1))/d;
                b = (pow(lam2,2)*(phi0-phi1)+pow(lam0,2)*(phi1-phi2)+pow(lam1,2)*(phi2-phi0))/d;
                c = (lam0*lam2*(lam2-lam0)*phi1+pow(lam1,2)*(lam2*phi0-lam0*phi2)+lam1*(phi2*pow(lam0,2)-phi0*pow(lam2,2)))/d;
                double disc = b*b - 4*a*c;
                if (disc >= 0) {
                    if (a != 0) {
                        // numerically stable evaluation of both roots
                        double q = -0.5*(b + copysign(sqrt(disc), b));
                        double r1 = q/a;
                        double r2 = (q != 0) ? c/q : r1;
                        if (r1*r2 < 0) lam = (r1 > r2) ? r1 : r2;
                        else {
                            mroot = true;
                            lam = (fabs(r1) < fabs(r2)) ? r1 : r2;
                            lroot = (fabs(r1) >= fabs(r2)) ? r1 : r2;
                        }
                    }
                    else if (b != 0) lam = -c/b;
                    else lam = 0;
                }
                else if (a != 0) {
                    // first try a least squares linear fit
                    double sx = lam0 + lam1 + lam2;
                    double sy = phi0 + phi1 + phi2;
                    double sx2 = pow(lam0,2) + pow(lam1,2) + pow(lam2,2);
                    double sxy = lam0*phi0 + lam1*phi1 + lam2*phi2;
                    double slope = (3*sxy - sx*sy)/(3*sx2 - pow(sx,2));
                    double intercept = (sy - slope*sx)/3;
                    lam = -intercept/slope;
                    // if that produces a negative root, use the minimum of the parabola
                    if (lam < 0) lam = -b/(2*a);
                }
                else
                    lam = (fabs(b) > 0) ? -c/b : 0;
            }
        }
        updateFv();
        ++iter;
    }
    if (!conv && (lam != lbest)) {
        // the yield surface may not be reachable along the fixed flow direction
        // (e.g., pressure-dependent yield criteria far from equilibrium), so fall
        // back to the iterate closest to the yield surface rather than the last one
        lam = lbest;
        updateFv();
        pe.m_F = Fv; pe.m_J = Fv.det();
        Kv = m_pCrit->DamageCriterion(pt);
        phi = Kv - Ky;
    }

    // restore the total deformation gradient
    pe.m_F = Ftmp; pe.m_J = Jtmp;

    Result r;
    r.Fvsi = Fs.inverse()*Fv;
    r.Kv = Kv;
    r.phi = phi;
    r.conv = conv;
    return r;
}

//-----------------------------------------------------------------------------
// Symmetric formulation based on (Fvs)^-1 = xi*(Fus)^-1*G (Eq. 3.18), with
// G = I - lam*Nhat symmetric (zero plastic spin). The yielded elastic deformation
// gradient Fv = xi*Fe*G is evaluated directly, without assuming that Fv shares
// the rotation of Fs, so no symmetrization is needed and Eq. 3.18 holds exactly.
FEReactivePlasticReturnMap::Result FEReactivePlasticReturnMap::Symmetric(FEMaterialPoint& pt, const mat3d& Fs, const mat3d& Fusi, double Ky, bool mandel)
{
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();

    // store safe copy of total deformation gradient
    mat3d Ftmp = pe.m_F; double Jtmp = pe.m_J;

    const mat3dd I(1);
    const double tol = m_rtol*Ky;
    const mat3d Fe = Fs*Fusi;     // elastic predictor

    mat3ds Nh; Nh.zero();
    mat3ds G = I;
    double xi = 1;
    mat3d Fv = Fe;
    double Kv = 0;

    // evaluate G, xi, Fv and the yield function for a given lam
    auto evaluate = [&](double lam) {
        G = I - Nh*lam;
        // xi enforces det(Fvs) = 1 when plastic deformation is isochoric (Eq. 3.19)
        xi = m_isochrc ? pow((Fusi*G).det(), -1./3.) : 1.0;
        Fv = Fe*G*xi;
        pe.m_F = Fv; pe.m_J = Fv.det();
        Kv = m_pCrit->DamageCriterion(pt);
        return Kv - Ky;
    };

    // flow direction, evaluated at the elastic predictor
    double phi = evaluate(0);
    if (mandel) {
        // classical associated flow: normal to the yield surface in Mandel stress space
        mat3ds dPhi = m_pCrit->CriterionStressTangent(pt);
        Nh = (Fv.inverse()*dPhi*Fv).sym();
    }
    else {
        // normal to the yield surface with respect to G, i.e., dPhi = Nh : dG
        mat3d M = YieldStressTangent(pt);
        Nh = (Fv.transpose()*M*Fv.transpose().inverse()).sym();
    }
    double Nmag = Nh.norm();
    bool conv = false;
    double lam = 0;
    if (Nmag > 0) {
        Nh = Nh/Nmag;

        // G must remain positive definite
        double ev[3];
        Nh.eigen2(ev);
        double emax = ev[0];
        if (ev[1] > emax) emax = ev[1];
        if (ev[2] > emax) emax = ev[2];

        // safeguarded Newton iterations on phi(lam) along the fixed flow direction,
        // bracketing the root in [lo, hi] and falling back to bisection
        double lo = 0;
        bool hasHi = (emax > 0);
        double hi = hasHi ? (1 - 1e-6)/emax : 0;
        double lbest = 0, phibest = phi;
        const int maxit = 50;
        for (int iter = 0; ; ++iter) {
            if (fabs(phi) <= tol) { conv = true; break; }
            if (iter >= maxit) break;

            // dphi/dlam = M : (dFv/dlam*Fv^-1), with
            // dFv/dlam*Fv^-1 = (dxi/dlam/xi)*I - Fe*Nh*G^-1*Fe^-1
            mat3d M = YieldStressTangent(pt);
            mat3ds Gi = G.inverse();
            double dxi = m_isochrc ? (Gi*Nh).trace()/3. : 0;
            mat3d L = mat3dd(dxi) - Fe*Nh*Gi*Fe.inverse();
            double slope = (M*L.transpose()).trace();

            // update the bracket
            if ((phi > 0) && (slope < 0)) lo = lam;
            else { hi = lam; hasHi = true; }   // root (if any) lies below lam

            double lnew = (slope != 0) ? lam - phi/slope : lam;
            if ((lnew <= lo) || (hasHi && (lnew >= hi)) || (lnew == lam))
                lnew = hasHi ? 0.5*(lo + hi) : 2*lam + 1e-6;
            lam = lnew;
            phi = evaluate(lam);
            if (fabs(phi) < fabs(phibest)) { phibest = phi; lbest = lam; }
        }

        if (!conv && (lam != lbest)) {
            // the yield surface may not be reachable along this flow direction;
            // use the iterate closest to the yield surface
            lam = lbest;
            phi = evaluate(lam);
        }
    }

    // restore the total deformation gradient
    pe.m_F = Ftmp; pe.m_J = Jtmp;

    Result r;
    r.Fvsi = Fusi*G*xi;
    r.Kv = Kv;
    r.phi = phi;
    r.conv = conv;
    return r;
}
