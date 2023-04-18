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
#include "FE3FieldElasticShellDomain.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FE3FieldElasticShellDomain, FEElasticShellDomain)
	ADD_PARAMETER(m_blaugon, "laugon");
	ADD_PARAMETER(m_augtol , "atol");
	ADD_PARAMETER(m_naugmin, "minaug");
	ADD_PARAMETER(m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
void FE3FieldElasticShellDomain::ELEM_DATA::Serialize(DumpStream& ar)
{
	ar & eJ;
	ar & ep;
	ar & Lk;
}

//-----------------------------------------------------------------------------
FE3FieldElasticShellDomain::FE3FieldElasticShellDomain(FEModel* pfem) : FEElasticShellDomain(pfem)
{
	m_blaugon = false;
	m_augtol = 0.01;
	m_naugmin = 0;
	m_naugmax = 0;
}

//-----------------------------------------------------------------------------
FE3FieldElasticShellDomain& FE3FieldElasticShellDomain::operator = (FE3FieldElasticShellDomain& d) 
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
bool FE3FieldElasticShellDomain::DoAugmentations() const
{
	return m_blaugon;
}

//-----------------------------------------------------------------------------
//! Initialize the 3-field domain data
bool FE3FieldElasticShellDomain::Init()
{
    // make sure the domain material uses an uncoupled formulation
    if (dynamic_cast<FEUncoupledMaterial*>(m_pMat) == 0) return false;
	if (FEElasticShellDomain::Init() == false) return false;
    
    // allocate element data
    int NE = Elements();
    m_Data.resize(NE);
    
    // initialize element data
    for (int i=0; i<NE; ++i)
    {
        ELEM_DATA& d = m_Data[i];
        d.eJ = 1.0;
        d.ep = 0.0;
        d.Lk = 0.0;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FE3FieldElasticShellDomain::Reset()
{
    FEElasticShellDomain::Reset();
    // initialize element data
    int NE = (int)m_Data.size();
    for (int i=0; i<NE; ++i)
    {
        ELEM_DATA& d = m_Data[i];
        d.eJ = 1.0;
        d.ep = 0.0;
        d.Lk = 0.0;
    }
}

//-----------------------------------------------------------------------------
//! Stiffness matrix for three-field domain
void FE3FieldElasticShellDomain::StiffnessMatrix(FELinearSystem& LS)
{
    FEModel& fem = *GetFEModel();
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FEShellElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);
        
        // create the element's stiffness matrix
        int ndof = 6*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate material and geometrical stiffness (i.e. constitutive component)
        ElementStiffness(iel, ke);

        // Calculate dilatational stiffness
        ElementDilatationalStiffness(fem, iel, ke);
        
        // get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates dilatational element stiffness component for element iel

void FE3FieldElasticShellDomain::ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke)
{
    int i, j, i6, j6, n;
    
    FEShellElement& elem = Element(iel);
    ELEM_DATA& ed = m_Data[iel];
    
    const int nint = elem.GaussPoints();
    const int neln = elem.Nodes();
    
    // get the material
    FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_pMat);
    assert(pmi);
    
    // average global derivatives
    vector<vec3d> gradMu(neln,vec3d(0,0,0));
    vector<vec3d> gradMd(neln,vec3d(0,0,0));
    
    // initial element volume
    double Ve = 0;
    
    // global derivatives of shape functions
    const double *gw = elem.GaussWeights();
    double eta;
    vec3d gcnt[3];

    // jacobian
    double Jt, J0;
    
    const double* Mr, *Ms, *M;

    // repeat over gauss-points
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        J0 = detJ0(elem, n);
        Jt = detJ(elem, n);
        
        Jt *= gw[n];
        
        Ve += J0*gw[n];
        
        eta = elem.gt(n);
        Mr = elem.Hr(n);
        Ms = elem.Hs(n);
        M  = elem.H(n);

        ContraBaseVectors(elem, n, gcnt);

        for (i=0; i<neln; ++i)
        {
            vec3d gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] += ((gradM*(1+eta) + gcnt[2]*M[i])/2)*Jt;
            gradMd[i] += ((gradM*(1-eta) - gcnt[2]*M[i])/2)*Jt;
        }
    }
    
    // get effective modulus
    double k = pmi->UJJ(ed.eJ);
    
    // next, we add the Lagrangian contribution
    // note that this term will always be zero if the material does not
    // use the augmented lagrangian
    k += ed.Lk*pmi->hpp(ed.eJ);
    
    // divide by initial volume
    k /= Ve;
    
    // calculate dilatational stiffness component
    for (i=0, i6=0; i<neln; ++i, i6 += 6)
    {
        for (j=0, j6 = 0; j<neln; ++j, j6 += 6)
        {
            mat3d Kuu = (gradMu[i] & gradMu[j])*k;
            mat3d Kud = (gradMu[i] & gradMd[j])*k;
            mat3d Kdu = (gradMd[i] & gradMu[j])*k;
            mat3d Kdd = (gradMd[i] & gradMd[j])*k;
            
            ke[i6  ][j6  ] += Kuu(0,0); ke[i6  ][j6+1] += Kuu(0,1); ke[i6  ][j6+2] += Kuu(0,2);
            ke[i6+1][j6  ] += Kuu(1,0); ke[i6+1][j6+1] += Kuu(1,1); ke[i6+1][j6+2] += Kuu(1,2);
            ke[i6+2][j6  ] += Kuu(2,0); ke[i6+2][j6+1] += Kuu(2,1); ke[i6+2][j6+2] += Kuu(2,2);
            
            ke[i6  ][j6+3] += Kud(0,0); ke[i6  ][j6+4] += Kud(0,1); ke[i6  ][j6+5] += Kud(0,2);
            ke[i6+1][j6+3] += Kud(1,0); ke[i6+1][j6+4] += Kud(1,1); ke[i6+1][j6+5] += Kud(1,2);
            ke[i6+2][j6+3] += Kud(2,0); ke[i6+2][j6+4] += Kud(2,1); ke[i6+2][j6+5] += Kud(2,2);
            
            ke[i6+3][j6  ] += Kdu(0,0); ke[i6+3][j6+1] += Kdu(0,1); ke[i6+3][j6+2] += Kdu(0,2);
            ke[i6+4][j6  ] += Kdu(1,0); ke[i6+4][j6+1] += Kdu(1,1); ke[i6+4][j6+2] += Kdu(1,2);
            ke[i6+5][j6  ] += Kdu(2,0); ke[i6+5][j6+1] += Kdu(2,1); ke[i6+5][j6+2] += Kdu(2,2);
            
            ke[i6+3][j6+3] += Kdd(0,0); ke[i6+3][j6+4] += Kdd(0,1); ke[i6+3][j6+5] += Kdd(0,2);
            ke[i6+4][j6+3] += Kdd(1,0); ke[i6+4][j6+4] += Kdd(1,1); ke[i6+4][j6+5] += Kdd(1,2);
            ke[i6+5][j6+3] += Kdd(2,0); ke[i6+5][j6+4] += Kdd(2,1); ke[i6+5][j6+5] += Kdd(2,2);
        }
    }
}

//-----------------------------------------------------------------------------
//! Calculates the shell element material and geometrical stiffness matrix

void FE3FieldElasticShellDomain::ElementStiffness(int iel, matrix& ke)
{
    FEShellElement& el = Element(iel);
    ELEM_DATA& ed = m_Data[iel];

    // get the material
    FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_pMat);
    assert(pmi);
    
    int i, i6, j, j6, n;
    
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    const double* Mr, *Ms, *M;
    vec3d gradMu[FEElement::MAX_NODES], gradMd[FEElement::MAX_NODES];
    
    // jacobian matrix determinant
    double detJt;
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    double eta;
    
    vec3d gcnt[3];
    
    // get the material
    FEUncoupledMaterial& mat = dynamic_cast<FEUncoupledMaterial&>(*m_pMat);
    
    // we need the following tensors for the dilational stiffness
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds Cp = IxI - I4*2;
    
    // calculate element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        // get the material's tangent
        // Note that we are only grabbing the deviatoric tangent.
        // The other tangent terms depend on the pressure p
        // which we seperately
        tens4ds C = Cp*ed.ep + mat.DevTangent(mp);
        
        // get the stress
        mat3ds s = pt.m_s;
        
        eta = el.gt(n);
        
        Mr = el.Hr(n);
        Ms = el.Hs(n);
        M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);
        
        // ------------ constitutive component --------------
        
        // setup the material point
        
        for (i=0; i<neln; ++i)
        {
            vec3d gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            gradMu[i] = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            gradMd[i] = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        }
        
        for (i=0, i6=0; i<neln; ++i, i6 += 6)
        {
            for (j=0, j6 = 0; j<neln; ++j, j6 += 6)
            {
                mat3d Kuu = vdotTdotv(gradMu[i], C, gradMu[j])*detJt;
                mat3d Kud = vdotTdotv(gradMu[i], C, gradMd[j])*detJt;
                mat3d Kdu = vdotTdotv(gradMd[i], C, gradMu[j])*detJt;
                mat3d Kdd = vdotTdotv(gradMd[i], C, gradMd[j])*detJt;
                
                ke[i6  ][j6  ] += Kuu(0,0); ke[i6  ][j6+1] += Kuu(0,1); ke[i6  ][j6+2] += Kuu(0,2);
                ke[i6+1][j6  ] += Kuu(1,0); ke[i6+1][j6+1] += Kuu(1,1); ke[i6+1][j6+2] += Kuu(1,2);
                ke[i6+2][j6  ] += Kuu(2,0); ke[i6+2][j6+1] += Kuu(2,1); ke[i6+2][j6+2] += Kuu(2,2);
                
                ke[i6  ][j6+3] += Kud(0,0); ke[i6  ][j6+4] += Kud(0,1); ke[i6  ][j6+5] += Kud(0,2);
                ke[i6+1][j6+3] += Kud(1,0); ke[i6+1][j6+4] += Kud(1,1); ke[i6+1][j6+5] += Kud(1,2);
                ke[i6+2][j6+3] += Kud(2,0); ke[i6+2][j6+4] += Kud(2,1); ke[i6+2][j6+5] += Kud(2,2);
                
                ke[i6+3][j6  ] += Kdu(0,0); ke[i6+3][j6+1] += Kdu(0,1); ke[i6+3][j6+2] += Kdu(0,2);
                ke[i6+4][j6  ] += Kdu(1,0); ke[i6+4][j6+1] += Kdu(1,1); ke[i6+4][j6+2] += Kdu(1,2);
                ke[i6+5][j6  ] += Kdu(2,0); ke[i6+5][j6+1] += Kdu(2,1); ke[i6+5][j6+2] += Kdu(2,2);
                
                ke[i6+3][j6+3] += Kdd(0,0); ke[i6+3][j6+4] += Kdd(0,1); ke[i6+3][j6+5] += Kdd(0,2);
                ke[i6+4][j6+3] += Kdd(1,0); ke[i6+4][j6+4] += Kdd(1,1); ke[i6+4][j6+5] += Kdd(1,2);
                ke[i6+5][j6+3] += Kdd(2,0); ke[i6+5][j6+4] += Kdd(2,1); ke[i6+5][j6+5] += Kdd(2,2);
            }
        }
        
        // ------------ initial stress component --------------
        
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                double Kuu = gradMu[i]*(s*gradMu[j])*detJt;
                double Kud = gradMu[i]*(s*gradMd[j])*detJt;
                double Kdu = gradMd[i]*(s*gradMu[j])*detJt;
                double Kdd = gradMd[i]*(s*gradMd[j])*detJt;
                
                // the u-u component
                ke[6*i  ][6*j  ] += Kuu;
                ke[6*i+1][6*j+1] += Kuu;
                ke[6*i+2][6*j+2] += Kuu;
                
                // the u-d component
                ke[6*i  ][6*j+3] += Kud;
                ke[6*i+1][6*j+4] += Kud;
                ke[6*i+2][6*j+5] += Kud;
                
                // the d-u component
                ke[6*i+3][6*j  ] += Kdu;
                ke[6*i+4][6*j+1] += Kdu;
                ke[6*i+5][6*j+2] += Kdu;
                
                // the d-d component
                ke[6*i+3][6*j+3] += Kdd;
                ke[6*i+4][6*j+4] += Kdd;
                ke[6*i+5][6*j+5] += Kdd;
            }
        
    } // end loop over gauss-points
    
}

//-----------------------------------------------------------------------------
//! This function loops over all elements and updates the stress
void FE3FieldElasticShellDomain::Update(const FETimeInfo& tp)
{
	FESSIShellDomain::Update(tp);

    bool berr = false;
    int NE = (int) m_Elem.size();
#pragma omp parallel for shared(NE, berr)
    for (int i=0; i<NE; ++i)
    {
        try
        {
            UpdateElementStress(i);
        }
        catch (NegativeJacobian e)
        {
#pragma omp critical
            {
                berr = true;
                if (e.DoOutput()) feLogError(e.what());
            }
        }
    }
    
    if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
//! This function updates the stresses for elements using the three-field formulation.
//! For such elements, the stress is a sum of a deviatoric stress, calculate by the
//! material and a dilatational term.
void FE3FieldElasticShellDomain::UpdateElementStress(int iel)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // get the material
    FEUncoupledMaterial& mat = *(dynamic_cast<FEUncoupledMaterial*>(m_pMat));
    
    // get the solid element
    FEShellElement& el = m_Elem[iel];
    ELEM_DATA& ed = m_Data[iel];
    
    // get the number of integration points
    int nint = el.GaussPoints();
    
    // get the integration weights
    double* gw = el.GaussWeights();
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    const int NELN = FEElement::MAX_NODES;
    vec3d r0[NELN], s0[NELN], r[NELN], s[NELN];
    vec3d v[NELN], w[NELN];
    vec3d a[NELN], b[NELN];
    // nodal coordinates
    GetCurrentNodalCoordinates(el, r, m_alphaf, false);
    GetCurrentNodalCoordinates(el, s, m_alphaf, true);
    GetReferenceNodalCoordinates(el, r0, false);
    GetReferenceNodalCoordinates(el, s0, true);

    // update dynamic quantities
    if (m_update_dynamic)
    {
        for (int j=0; j<neln; ++j)
        {
            FENode& node = m_pMesh->Node(el.m_node[j]);
            v[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*m_alphaf + node.m_vp*(1-m_alphaf);
            w[j] = node.get_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2])*m_alphaf + node.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2])*(1-m_alphaf);
            a[j] = node.m_at*m_alpham + node.m_ap*(1-m_alpham);
            b[j] = node.get_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2])*m_alpham + node.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2])*(1-m_alpham);
        }
    }

    // calculate the average dilatation and pressure
    double vt = 0, V = 0;
    for (int n=0; n<nint; ++n)
    {
        vt += detJ(el, n)*gw[n];
        V += detJ0(el, n)*gw[n];
    }
    
    // calculate volume ratio
    ed.eJ = vt / V;
    
    // Calculate pressure. This is a sum of a Lagrangian term and a penalty term
    //      <--- Lag. mult. -->  <-- penalty -->
    ed.ep = ed.Lk*mat.hp(ed.eJ) + mat.UJ(ed.eJ);
    //    ed.ep = mat.UJ(ed.eJ);
    
    // loop over the integration points and calculate
    // the stress at the integration point
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        pt.m_p = ed.ep;
        
        // material point coordinates
        // TODO: I'm not entirly happy with this solution
        //         since the material point coordinates are not used by most materials.
        mp.m_r0 = evaluate(el, r0, s0, n);
        mp.m_rt = evaluate(el, r, s, n);

        // get the deformation gradient and determinant at intermediate time
        double Jt, Jp;
        mat3d Ft, Fp;
        Jt = defgrad(el, Ft, n);
        Jp = defgradp(el, Fp, n);
        pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
        pt.m_J = pt.m_F.det();
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (Ft - Fp)*Fi/dt;
        if (m_update_dynamic)
        {
            pt.m_v = evaluate(el, v, w, n);
            pt.m_a = evaluate(el, a, b, n);
        }
        
        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the stress at this material point
        // Note that we don't call the material's Stress member function.
        // The reason is that we need to use the averaged pressure for the element
        // and the Stress function uses the pointwise pressure.
        // Therefore we call the DevStress function and add the pressure term
        // seperately.
        pt.m_s = mat3dd(ed.ep) + mat.DevStress(mp);
        
        
        // adjust stress for strain energy conservation
        if (m_alphaf == 0.5)
        {
            // evaluate strain energy at current time
			mat3d Ftmp = pt.m_F;
			double Jtmp = pt.m_J;
			pt.m_F = Ft;
            pt.m_J = Jt;
            FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(m_pMat);
            pt.m_Wt = pme->StrainEnergyDensity(mp);
			pt.m_F = Ftmp;
			pt.m_J = Jtmp;

            mat3ds D = pt.m_L.sym();
            double D2 = D.dotdot(D);
            if (D2 > 0)
                pt.m_s += D*(((pt.m_Wt-pt.m_Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
        }
    }
}

//-----------------------------------------------------------------------------
//! Do augmentation
bool FE3FieldElasticShellDomain::Augment(int naug)
{
    FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_pMat);
    assert(pmi);
    
    // make sure Augmented Lagrangian flag is on
    if (m_blaugon == false) return true;
    
    // do the augmentation
    int n;
    double normL0 = 0, normL1 = 0, L0, L1;
    double k = pmi->m_K;
    int NE = Elements();
    
    for (n=0; n<NE; ++n)
    {
        ELEM_DATA& ed = m_Data[n];
        
        L0 = ed.Lk;
        normL0 += L0*L0;
        
        L1 = L0 + k*pmi->h(ed.eJ);
        normL1 += L1*L1;
    }
    
    normL0 = sqrt(normL0);
    normL1 = sqrt(normL1);
    
    // check convergence
    double pctn = 0;
    if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
    
    feLog(" material %d\n", pmi->GetID());
    feLog("                        CURRENT         CHANGE        REQUIRED\n");
    feLog("   pressure norm : %15le%15le%15le\n", normL1, pctn, m_augtol);
    
    // check convergence
    bool bconv = true;
    if (pctn >= m_augtol) bconv = false;
    if (m_naugmin > naug) bconv = false;
    if ((m_naugmax > 0) && (m_naugmax <= naug)) bconv = true;
    
    // do the augmentation only if we have not yet converged
    if (bconv == false)
    {
        for (n=0; n<NE; ++n)
        {
            ELEM_DATA& ed = m_Data[n];
            
            ed.Lk += k*pmi->h(ed.eJ);
            ed.ep = ed.Lk*pmi->hp(ed.eJ) + k*log(ed.eJ)/ed.eJ;
        }
    }
    
    return bconv;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FE3FieldElasticShellDomain::Serialize(DumpStream &ar)
{
    FEElasticShellDomain::Serialize(ar);
	ar & m_Data;
}
