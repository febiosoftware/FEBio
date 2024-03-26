/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2024 University of Utah, The Trustees of Columbia University in
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
#include "FEElasticSoluteSolidDomain.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>
#include <FECore/sys.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FEElasticSoluteSolidDomain::FEElasticSoluteSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticSoluteDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dof(pfem)
{
    m_pMat = nullptr;
	
    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
    }
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
    m_pMat = dynamic_cast<FEElasticSolute*>(pmat);
    assert(m_pMat);
}

//-----------------------------------------------------------------------------
// get total dof list
const FEDofList& FEElasticSoluteSolidDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data.
void FEElasticSoluteSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
    // get nodal DOFS
    const int n_sol = m_pMat->Solutes();
    
    int n_m = el.Nodes();
    int n_dpn = 3 + n_sol;
    size_t sz_ln = n_m * (n_dpn + 3);
    lm.resize(sz_ln);
    
    for (int i_a = 0; i_a < n_m; ++i_a)
    {
        int a = el.m_node[i_a];
        FENode& node = m_pMesh->Node(a);
        
        vector<int>& id = node.m_ID;
        
        int lm_u_a = n_dpn * i_a;
        // first the displacement dofs
        lm[lm_u_a  ] = id[m_dofU[0]];
        lm[lm_u_a+1] = id[m_dofU[1]];
        lm[lm_u_a+2] = id[m_dofU[2]];

        int lm_c_a = n_dpn * i_a + 3;
        // concentration dofs
        for (int isol = 0; isol < n_sol; ++isol)
            lm[lm_c_a+isol] = id[m_dofC + m_pMat->GetSolute(isol)->GetSoluteDOF()];
        
        // rigid rotational dofs
        // TODO: Do we really need this?
        int lm_r_a = n_dpn * i_a + 3 * i_a;
        lm[lm_r_a  ] = id[m_dofR[0]];
        lm[lm_r_a+1] = id[m_dofR[1]];
        lm[lm_r_a+2] = id[m_dofR[2]];
    }
}

//-----------------------------------------------------------------------------
bool FEElasticSoluteSolidDomain::Init()
{
    // initialize base class
	if (FESolidDomain::Init() == false) return false;
    
    // extract the initial concentrations of the solid-bound molecules
    const int n_sol = m_pMat->Solutes();
    int n_e = m_Elem.size();

    for (int i_e = 0; i_e < n_e; ++i_e)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i_e];
        
        // get the number of integration points
        int n_int = el.GaussPoints();
        
        // loop over the integration points
        for (int i_k = 0; i_k < n_int; ++i_k)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            // initialize elastic solute solutes
            ps.m_nsol = n_sol;
            ps.m_c.assign(n_sol,0);
            ps.m_ca.assign(n_sol,0);
            ps.m_crp.assign(n_sol, 0);
            ps.m_gradc.assign(n_sol,vec3d(0,0,0));
            ps.m_j.assign(n_sol,vec3d(0,0,0));
            
            if (m_pMat->m_phi0(mp) > 1.0) {
                feLogError("Referential solid volume fraction of multiphasic material cannot exceed unity!\nCheck ratios of sbm apparent and true densities.");
                return false;
            }
        }
    }

	// set the active degrees of freedom list
	FEDofList dofs(GetFEModel());
    for (int i_sol = 0; i_sol < n_sol; ++i_sol)
	{
		int i_sol_dof = m_pMat->GetSolute(i_sol)->GetSoluteDOF();
		dofs.AddDof(m_dofC + i_sol_dof);
	}
	m_dof = dofs;

    return true;
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::Activate()
{
    for (int i = 0; i < Nodes(); ++i)
    {
        FENode& node = Node(i);
        if (node.HasFlags(FENode::EXCLUDE) == false)
        {
            if (node.m_rid < 0)
            {
                node.set_active(m_dofU[0]);
                node.set_active(m_dofU[1]);
                node.set_active(m_dofU[2]);
            }
        }
    }
    
    const int n_sol = m_pMat->Solutes();

    // Activate dof_P and dof_C, except when a solid element is connected to the
    // back of a shell element, in which case activate dof_Q and dof_D for those nodes.
    FEMesh& mesh = *GetMesh();
    for (int i_e = 0; i_e < Elements(); ++i_e) {
        FESolidElement& el = m_Elem[i_e];
        int n_m = el.Nodes();
        for (int i_a = 0; i_a < n_m; ++i_a)
        {
            FENode& node = mesh.Node(el.m_node[i_a]);
            for (int i_sol = 0; i_sol < n_sol; ++i_sol)
                node.set_active(m_dofC + m_pMat->GetSolute(i_sol)->GetSoluteDOF());
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::InitMaterialPoints()
{
    const int n_sol = m_pMat->Solutes();
    FEMesh& mesh = *GetMesh();
   
    const int NE = FEElement::MAX_NODES;
    vector< vector<double> > c0(n_sol, vector<double>(NE));
    vector<int> sid(n_sol);
    for (int i_sol = 0; i_sol < n_sol; ++i_sol)
        sid[i_sol] = m_pMat->GetSolute(i_sol)->GetSoluteDOF();
    
    int n_e = m_Elem.size();
    for (int i_e = 0; i_e < n_e; ++i_e)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i_e];
        
        // get the number of nodes
        int n_m = el.Nodes();
        // get initial values of solute concentration
        for (int i_a = 0; i_a < n_m; ++i_a) {
            FENode& na = mesh.Node(el.m_node[i_a]);
            for (int i_sol = 0; i_sol < n_sol; ++i_sol)
                c0[i_sol][i_a] = na.get(m_dofC + sid[i_sol]);
        }
        
        // get the number of integration points
        int n_int = el.GaussPoints();
        
        // loop over the integration points
        for (int i_k = 0; i_k < n_int; ++i_k)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
            FEElasticMaterialPoint& pm = *(mp.ExtractData<FEElasticMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            // initialize multiphasic solutes
            ps.m_nsol = n_sol;
            
            // initialize effective solute concentrations
            for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
                ps.m_c[i_sol] = el.Evaluate(c0[i_sol], i_k);
                ps.m_gradc[i_sol] = gradient(el, c0[i_sol], i_k);
            }
            
            for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
                ps.m_j[i_sol] = m_pMat->SoluteFlux(mp, i_sol);
                ps.m_crp[i_sol] = m_pMat->Porosity(mp) * ps.m_ca[i_sol];
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::Reset()
{
    // reset base class
    FESolidDomain::Reset();
    
    const int n_sol = m_pMat->Solutes();
    
    int n_e = m_Elem.size();
    for (int i_e = 0; i_e < n_e; ++i_e)
    {
        // get the solid element
        FESolidElement& el = m_Elem[i_e];
        
        // get the number of integration points
        int n_int = el.GaussPoints();
        
        // loop over the integration points
        for (int i_k = 0; i_k < n_int; ++i_k)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
                        
            // initialize multiphasic solutes
            ps.m_nsol = n_sol;
            ps.m_c.assign(n_sol,0);
            ps.m_ca.assign(n_sol,0);
            ps.m_crp.assign(n_sol, 0);
            ps.m_gradc.assign(n_sol,vec3d(0,0,0));
            ps.m_j.assign(n_sol,vec3d(0,0,0));

            // reset chemical reaction element data
            ps.m_cri.clear();
            ps.m_crd.clear();
            int n_r = m_pMat->Reactions();
            for (int i_r = 0; i_r < n_r; ++i_r)
                m_pMat->GetReaction(i_r)->ResetElementData(mp);
        }
    }
    m_breset = true;
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    FESolidDomain::PreSolveUpdate(timeInfo);
    
    const int max_eln = FEElement::MAX_NODES;
    vec3d x0[max_eln], xt[max_eln], r0, rt;
    FEMesh& mesh = *GetMesh();
    size_t n_m = m_Elem.size();
    for (size_t i_e = 0; i_e < n_m; ++i_e)
    {
        FESolidElement& el = m_Elem[i_e];
        int tot_eln = el.Nodes();
        for (int i_a = 0; i_a < n_m; ++i_a)
        {
            x0[i_a] = mesh.Node(el.m_node[i_a]).m_r0;
            xt[i_a] = mesh.Node(el.m_node[i_a]).m_rt;
        }
        
        int n_int = el.GaussPoints();
        for (int i_k = 0; i_k < n_int; ++i_k)
        {
            r0 = el.Evaluate(x0, i_k);
            rt = el.Evaluate(xt, i_k);
            
            FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
            FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            mp.m_r0 = r0;
            mp.m_rt = rt;
            
            pe.m_J = defgrad(el, pe.m_F, i_k);
            
            // reset referential actual solute concentration at previous time
            int n_sol = m_pMat->Solutes();
            for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
                ps.m_crp[i_sol] = m_pMat->Porosity(mp) * ps.m_ca[i_sol];
            }
                        
            // reset chemical reaction element data
            int n_r = m_pMat->Reactions();
            for (int i_r = 0; i_r < n_r; ++i_r)
                m_pMat->GetReaction(i_r)->InitializeElementData(mp);
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::InternalForces(FEGlobalVector& R)
{
    size_t n_e = m_Elem.size();
    
    // get nodal DOFS
    int n_sol = m_pMat->Solutes();
    int n_dpn = 3 + n_sol;
    
#pragma omp parallel for
    for (int i_e = 0; i_e < n_e; ++i_e)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i_e];
        
        // get the element force vector and initialize it to zero
        int ndof = n_dpn * el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEElasticSoluteSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradNa;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *N;
    
    int n_int = el.GaussPoints();
    int n_m = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    const int n_sol = m_pMat->Solutes();
    int n_dpn = 3 + n_sol;
    
    const int n_r = m_pMat->Reactions();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // repeat for all integration points
    for (int i_k = 0; i_k < n_int; ++i_k)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
        
        // calculate the jacobian
        detJt = invjact(el, Ji, i_k);
        
        Gr = el.Gr(i_k);
        Gs = el.Gs(i_k);
        Gt = el.Gt(i_k);
        
        N = el.H(i_k);
        
        // get the stress for this integration point
        s = pt.m_s;
        
        vector<vec3d> j(spt.m_j);
        
        // evaluate the porosity, its derivative w.r.t. J, and its gradient
        double phiw = m_pMat->Porosity(mp);
        vector<double> chat(n_sol,0);
        
        // get the solvent supply
        double phiwhat = 0;
        
        // chemical reactions
        for (int i_r = 0; i_r < n_r; ++i_r) {
            FEChemicalReaction* pri = m_pMat->GetReaction(i_r);
            double zhat = pri->ReactionSupply(mp);
            for (int i_sol = 0; i_sol < n_sol; ++i_sol)
                chat[i_sol] += phiw * zhat * pri->m_v[i_sol];
        }
        
        for (int i_a = 0; i_a < n_m; ++i_a)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradNa = vec3d(Ji[0][0]*Gr[i_a]+Ji[1][0]*Gs[i_a]+Ji[2][0]*Gt[i_a],
                          Ji[0][1]*Gr[i_a]+Ji[1][1]*Gs[i_a]+Ji[2][1]*Gt[i_a],
                          Ji[0][2]*Gr[i_a]+Ji[1][2]*Gs[i_a]+Ji[2][2]*Gt[i_a]);
            
            // calculate internal force
            vec3d fu = s * gradNa;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            int i_fe_sx = n_dpn * i_a;
            int i_fe_sy = n_dpn * i_a + 1;
            int i_fe_sz = n_dpn * i_a + 2;
            int i_fe_c = n_dpn * i_a + 3;
            fe[i_fe_sx] -= fu.x * detJt * gw[i_k];
            fe[i_fe_sy] -= fu.y * detJt * gw[i_k];
            fe[i_fe_sz] -= fu.z * detJt * gw[i_k];
            for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
                double fc = 0.0;
                fc += gradNa * j[i_sol] + N[i_a] * (chat[i_sol] - (phiw * spt.m_ca[i_sol] - spt.m_crp[i_sol]) / dt);
                fc *= dt;
                fe[i_fe_c+i_sol] -= fc * detJt * gw[i_k];
            }
                
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
    const int n_sol = m_pMat->Solutes();
    int n_dpn = 3 + n_sol;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    
#pragma omp parallel for
    for (int iel=0; iel<NE; ++iel)
    {
		FESolidElement& el = m_Elem[iel];

        // element stiffness matrix
        FEElementMatrix ke(el);

        // allocate stiffness matrix
        int neln = el.Nodes();
        int ndof = neln*n_dpn;
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementElasticSoluteStiffness(el, ke, bsymm);

		// get the lm vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEElasticSoluteSolidDomain::ElementElasticSoluteStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
    int n_int = el.GaussPoints();
    int n_m = el.Nodes();
    
    double *Gr, *Gs, *Gt, *N;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(n_m);
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    const int n_sol = m_pMat->Solutes();
    int n_dpn = 3 + n_sol;
    
    const int n_r = m_pMat->Reactions();
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (int i_k = 0; i_k < n_int; ++i_k)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, i_k);
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        Gr = el.Gr(i_k);
        Gs = el.Gs(i_k);
        Gt = el.Gt(i_k);
        
        N = el.H(i_k);
        
        // calculate global gradient of shape functions
        for (int i_a = 0; i_a < n_m; ++i_a)
            gradN[i_a] = g1*Gr[i_a] + g2*Gs[i_a] + g3*Gt[i_a];
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // next we get the determinant
        double J = ept.m_J;
                
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        
        // evaluate the porosity
        double phiw = m_pMat->Porosity(mp);
        double phi0 = m_pMat->m_phi0(mp);
        double phis = 1.0 - phiw;
        
        vector<mat3ds> D(n_sol);
        vector<mat3ds> ImD(n_sol);
        mat3dd I(1);
        
        // evaluate the solvent supply and its derivatives
        vector<double> Phic(n_sol,0);
        vector<mat3ds> dchatde(n_sol);
        
        // chemical reactions
		vector<double> reactionSupply(n_r, 0.0);
		vector< vector<double> > tangentReactionSupplyConcentration(n_r, vector<double>(n_sol));
        for (int i_r = 0; i_r < n_r; ++i_r)
		{
			FEChemicalReaction* reacti = m_pMat->GetReaction(i_r);

			reactionSupply[i_r] = reacti->ReactionSupply(mp);

			for (int i_sol = 0; i_sol < n_sol; ++i_sol)
				tangentReactionSupplyConcentration[i_r][i_sol] = reacti->Tangent_ReactionSupply_Concentration(mp, i_sol);
		}
        
        for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
			FESolute* soli = m_pMat->GetSolute(i_sol);

            // evaluate the diffusivity tensor
            D[i_sol] = soli->m_pDiff->Diffusivity(mp);
                                    
            // chemical reactions
            dchatde[i_sol].zero();
            //SL: this will create a dependence on the strain. Maybe the correct tangent reaction supply would be one that zeros this out?
            //SL: for now zero this.
            //SL: Check thisif (!m_pMat->m_bool_refC) {
                for (int i_r = 0; i_r < n_r; ++i_r) {
                    FEChemicalReaction* reacti = m_pMat->GetReaction(i_r);
                    dchatde[i_sol] += reacti->m_v[i_sol] * (I * reactionSupply[i_r]);
                    Phic[i_sol] += phiw* reacti->m_Vbar*tangentReactionSupplyConcentration[i_r][i_sol];
                }
            //}
        }
        
        // Miscellaneous constants
        double R = m_pMat->m_Rgas;
        double T = m_pMat->m_Tabs;
        
        // calculate all the matrices
        vec3d vtmp,gp,qpu;
        vector<vec3d> gc(n_sol),qcu(n_sol),wc(n_sol),jce(n_sol);
        vector< vector<vec3d> > jc(n_sol, vector<vec3d>(n_sol));
        vector<mat3d> ju(n_sol);
        vector< vector<double> > qcc(n_sol, vector<double>(n_sol));
        vector< vector<double> > dchatdc(n_sol, vector<double>(n_sol));
        for (int i_a = 0; i_a < n_m; ++i_a)
        {
            for (int i_b = 0; i_b < n_m; ++i_b)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i_a] * (s * gradN[i_b])) + vdotTdotv(gradN[i_a], C, gradN[i_b])) * detJ * gw[i_k];
                ke.add(n_dpn*i_a, n_dpn*i_b, Kuu);

                // do nothing for the kcu matrix?
                //for (int isol = 0; isol < n_sol; ++isol) {
                //    ke[n_dpn*a+3+isol][n_dpn*b] += 0.0;
                //    ke[n_dpn*a+3+isol][n_dpn*b+1] += 0.0;
                //    ke[n_dpn*a+3+isol][n_dpn*b+2] += 0.0;
                    
                    // do nothing for the kuc matrix?
                    //ke[n_dpn*a][n_dpn*b+3+isol] += 0.0;
                    //ke[n_dpn*a+1][n_dpn*b+3+isol] += 0.0;
                    //ke[n_dpn*a+2][n_dpn*b+3+isol] += 0.0;                    
                //}
                
                // calculate data for the kcc matrix
                for (int isol = 0; isol < n_sol; ++isol) {
                    for (int jsol = 0; jsol < n_sol; ++jsol) {
                        if (jsol = isol)
                            jc[isol][jsol] = -(D[isol] * (gradN[i_b] * (phiw)));
                    }
                }
                
                // calculate the kcc matrix
                int ke_c_a = n_dpn * i_a + 3;
                int ke_c_b = n_dpn * i_b + 3;
                for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
                    for (int j_sol = 0; j_sol < n_sol; ++j_sol) {
                        ke[ke_c_a+i_sol][ke_c_b+j_sol] += (gradN[i_a] * jc[i_sol][j_sol] + N[i_a] * (qcc[i_sol][j_sol] + N[i_b] * phiw * dchatdc[i_sol][j_sol])) * (detJ * gw[i_k] * dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    double tmp;
    if (bsymm) {
        for (int i = 0; i < n_dpn * n_m; ++i)
            for (int j = i + 1; j < n_dpn * n_m; ++j) {
                tmp = 0.5 * (ke[i][j] + ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FEElasticSoluteSolidDomain::Update(const FETimeInfo& tp)
{
    FEModel& fem = *GetFEModel();
    bool berr = false;
    int n_e = (int) m_Elem.size();
    double dt = fem.GetTime().timeIncrement;
#pragma omp parallel for shared(n_e, berr)
    for (int i_e = 0; i_e < n_e; ++i_e)
    {
        try
        {
            UpdateElementStress(i_e, dt);
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
void FEElasticSoluteSolidDomain::UpdateElementStress(int i_e, double dt)
{
    double* gw;
    vec3d r0[FEElement::MAX_NODES];
    vec3d rt[FEElement::MAX_NODES];
    
    // get the mesh
    FEMesh& mesh = *m_pMesh;
    
    // get the elastic solute material and dependencies
    FEElasticSolute* pmb = m_pMat;
    // get the number of solutes
    const int n_sol = (int)pmb->Solutes();
    // get the number of reactions
    int n_r = m_pMat->Reactions();

    // get the solid element and dependencies
    FESolidElement& el = m_Elem[i_e];
    // get the number of integration points
    int n_int = el.GaussPoints();
    // get the number of nodes
    int n_m = el.Nodes();
    // get the integration weights
    gw = el.GaussWeights();

    // get the solute ids
    vector< vector<double> > ct(n_sol, vector<double>(FEElement::MAX_NODES));
    vector<int> sid(n_sol);
    for (int i_sol = 0; i_sol < n_sol; ++i_sol)
        sid[i_sol] = pmb->GetSolute(i_sol)->GetSoluteDOF();
    
    // get the nodal solute data at the current time
    for (int i_a = 0; i_a < n_m; ++i_a)
    {
        FENode& node = mesh.Node(el.m_node[i_a]);
        r0[i_a] = node.m_r0;
        rt[i_a] = node.m_rt;
        for (int i_sol = 0; i_sol < n_sol; ++i_sol)
            ct[i_sol][i_a] = node.get(m_dofC + sid[i_sol]);
    }
    
    // loop over the integration points and calculate
    // the stress at the integration point
    for (int i_k = 0; i_k < n_int; ++i_k)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(i_k);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        // material point coordinates
        // TODO: I'm not entirly happy with this solution
        //		 since the material point coordinates are used by most materials.
        mp.m_r0 = el.Evaluate(r0, i_k);
        mp.m_rt = el.Evaluate(rt, i_k);
        
        // get the deformation gradient and determinant
        pt.m_J = defgrad(el, pt.m_F, i_k);
        mat3d Fp;
        defgradp(el, Fp, i_k);
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (pt.m_F - Fp) * Fi / dt;

        // project nodal data to the elasitc solute material points
        FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
               
        for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
            // evaluate effective solute concentrations at gauss-point
            spt.m_c[i_sol] = el.Evaluate(&ct[i_sol][0], i_k);
            // calculate the gradient of c at gauss-point
            spt.m_gradc[i_sol] = gradient(el, &ct[i_sol][0], i_k);
        }
        
        // update the solute fluxes and evaluate the actual solute concentration
        for (int i_sol = 0; i_sol < n_sol; ++i_sol) {
            spt.m_ca[i_sol] = pmb->Concentration(mp,i_sol);
            spt.m_j[i_sol] = pmb->SoluteFlux(mp,i_sol);
        }

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // update the stress
        pt.m_s = pmb->Stress(mp);
                
        // update chemical reaction element data
        for (int i_r = 0; i_r < n_r; ++i_r)
            pmb->GetReaction(i_r)->UpdateElementData(mp);
        
    }
    if (m_breset) m_breset = false;
}
