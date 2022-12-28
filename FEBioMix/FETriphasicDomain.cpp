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
#include "FETriphasicDomain.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FETriphasicDomain::FETriphasicDomain(FEModel* pfem) : FESolidDomain(pfem), FEElasticDomain(pfem), m_dofU(pfem), m_dofR(pfem), m_dof(pfem)
{
	m_pMat = nullptr;

    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
        m_dofP = pfem->GetDOFIndex("p");
        m_dofC = pfem->GetDOFIndex("concentration", 0);
    }
}

//-----------------------------------------------------------------------------
//! get the material (overridden from FEDomain)
FEMaterial* FETriphasicDomain::GetMaterial()
{
	return m_pMat;
}

//-----------------------------------------------------------------------------
//! get the total dof
const FEDofList& FETriphasicDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FETriphasic*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FETriphasicDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int dofc0 = m_dofC + m_pMat->m_pSolute[0]->GetSoluteDOF();
	int dofc1 = m_dofC + m_pMat->m_pSolute[1]->GetSoluteDOF();

	int N = el.Nodes();
	lm.resize(N*9);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);

		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[6*i  ] = id[m_dofU[0]];
		lm[6*i+1] = id[m_dofU[1]];
		lm[6*i+2] = id[m_dofU[2]];

		// now the pressure dofs
		lm[6*i+3] = id[m_dofP];
        
        // concentration dofs
        lm[6*i + 4] = id[dofc0];
        lm[6*i + 5] = id[dofc1];
        
        // rigid rotational dofs
		// TODO: Do I really need these?
		lm[6*N + 3*i  ] = id[m_dofR[0]];
		lm[6*N + 3*i+1] = id[m_dofR[1]];
		lm[6*N + 3*i+2] = id[m_dofR[2]];
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Activate()
{
	int dofc0 = m_dofC + m_pMat->m_pSolute[0]->GetSoluteDOF();
	int dofc1 = m_dofC + m_pMat->m_pSolute[1]->GetSoluteDOF();

	for (int i=0; i<Nodes(); ++i)
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

			node.set_active(m_dofP);
			node.set_active(dofc0 );
			node.set_active(dofc1 );
		}
	}

	// get the triphasic material
	FETriphasic* pmb = m_pMat;
	const int nsol = 2;
	const int nsbm = 0;

	const int NE = FEElement::MAX_NODES;
	double p0[NE];
	vector< vector<double> > c0(nsol, vector<double>(NE));
	FEMesh& m = *GetMesh();

	int id[2] = { m_pMat->m_pSolute[0]->GetSoluteID() - 1, m_pMat->m_pSolute[1]->GetSoluteID() - 1 };

	for (int i = 0; i<(int)m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the number of nodes
		int neln = el.Nodes();
		// get initial values of fluid pressure and solute concentrations
		for (int i = 0; i<neln; ++i)
		{
			p0[i] = m.Node(el.m_node[i]).get(m_dofP);
			for (int isol = 0; isol<nsol; ++isol)
				c0[isol][i] = m.Node(el.m_node[i]).get(m_dofC + id[isol]);
		}

		// get the number of integration points
		int nint = el.GaussPoints();

		// loop over the integration points
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEElasticMaterialPoint& pm = *(mp.ExtractData<FEElasticMaterialPoint>());
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());

			// initialize referential solid volume fraction
			pt.m_phi0t = pmb->m_phi0(mp);

			// initialize effective fluid pressure, its gradient, and fluid flux
			pt.m_p = el.Evaluate(p0, n);
			pt.m_gradp = gradient(el, p0, n);
			pt.m_w = pmb->FluidFlux(mp);


			// initialize multiphasic solutes
			ps.m_nsol = nsol;
			ps.m_nsbm = nsbm;

			// initialize effective solute concentrations
			for (int isol = 0; isol<nsol; ++isol) {
				ps.m_c[isol] = el.Evaluate(c0[isol], n);
				ps.m_gradc[isol] = gradient(el, c0[isol], n);
			}

			ps.m_psi = pmb->ElectricPotential(mp);
			for (int isol = 0; isol<nsol; ++isol) {
				ps.m_ca[isol] = pmb->Concentration(mp, isol);
				ps.m_j[isol] = pmb->SoluteFlux(mp, isol);
				ps.m_crp[isol] = pm.m_J*m_pMat->Porosity(mp)*ps.m_ca[isol];
			}
			pt.m_pa = pmb->Pressure(mp);
			ps.m_cF = pmb->FixedChargeDensity(mp);
			ps.m_Ie = pmb->CurrentDensity(mp);

			pm.m_s = pmb->Stress(mp);

		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Reset()
{
	// reset base class
	FESolidDomain::Reset();
	
	// get the multiphasic material
	FETriphasic* pmb = m_pMat;
	const int nsol = 2;
	const int nsbm = 0;
	
	// loop over all material points
	ForEachMaterialPoint([=](FEMaterialPoint& mp) 
	{
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
			
		// initialize referential solid volume fraction
		pt.m_phi0 = pt.m_phi0t = pmb->m_phi0(mp);
			
		// initialize multiphasic solutes
		ps.m_nsol = nsol;
		ps.m_c.assign(nsol,0);
		ps.m_ca.assign(nsol,0);
        ps.m_crp.assign(nsol, 0);
		ps.m_gradc.assign(nsol,vec3d(0,0,0));
        ps.m_bsb.assign(nsol, false);
		ps.m_k.assign(nsol, 0);
		ps.m_dkdJ.assign(nsol, 0);
		ps.m_dkdc.resize(nsol, vector<double>(nsol,0));
		ps.m_j.assign(nsol,vec3d(0,0,0));
		ps.m_nsbm = nsbm;
	});
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	FESolidDomain::PreSolveUpdate(timeInfo);

	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
    double pn[NE], p;
	FEMesh& m = *GetMesh();
	for (size_t iel=0; iel<m_Elem.size(); ++iel)
	{
		FESolidElement& el = m_Elem[iel];
		int neln = el.Nodes();
		for (int i=0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
			xt[i] = m.Node(el.m_node[i]).m_rt;
            pn[i] = m.Node(el.m_node[i]).get(m_dofP);
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);
            p = el.Evaluate(pn, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(mp.ExtractData<FESolutesMaterialPoint>());
            
            mp.m_r0 = r0;
			mp.m_rt = rt;

			pe.m_J = defgrad(el, pe.m_F, j);

            // reset referential solid volume fraction at previous time
            pt.m_phi0p = pt.m_phi0t;
            
            // reset determinant of solid deformation gradient at previous time
            pt.m_Jp = pe.m_J;
            
            pt.m_p = p;
            pt.m_gradp = gradient(el, pn, j);
            
            // reset referential actual solute concentration at previous time
            for (int j=0; j<2; ++j) {
                ps.m_crp[j] = pe.m_J*m_pMat->Porosity(mp)*ps.m_ca[j];
            }
			mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalForces(FEGlobalVector& R)
{
	size_t NE = m_Elem.size();
	#pragma omp parallel for shared (NE)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> lm;
		
		// get the element
		FESolidElement& el = m_Elem[i];

		// get the element force vector and initialize it to zero
		int ndof = 6*el.Nodes();
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

void FETriphasicDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3], detJt;

    vec3d gradN;
    mat3ds s;

	const double* Gr, *Gs, *Gt, *H;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(el.GetMaterialPoint(n)->ExtractData<FESolutesMaterialPoint>());
        
		// calculate the jacobian
		detJt = invjact(el, Ji, n);

		detJt *= gw[n];

		Gr = el.Gr(n);
		Gs = el.Gs(n);
		Gt = el.Gt(n);

        H = el.H(n);
        
        // next we get the determinant
        double Jp = bpt.m_Jp;
        double J = pt.m_J;
        
        // and then finally
        double divv = ((J-Jp)/dt)/J;
        
        // get the stress for this integration point
        s = pt.m_s;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solute flux
        vec3d j[2] = {spt.m_j[0],spt.m_j[1]};
        // get the charge number
        int z[2] = {m_pMat->m_pSolute[0]->ChargeNumber(), m_pMat->m_pSolute[1]->ChargeNumber()};
        
        vec3d je = j[0]*z[0] + j[1]*z[1];
        
        // evaluate the porosity, its derivative w.r.t. J, and its gradient
        double phiw = m_pMat->Porosity(mp);
        
        for (i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
            gradN = vec3d(Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i],
                          Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i],
                          Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i]);
            
            // calculate internal force
            vec3d fu = s*gradN;

            // the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
            fe[6*i  ] -= fu.x*detJt;
            fe[6*i+1] -= fu.y*detJt;
            fe[6*i+2] -= fu.z*detJt;
            fe[6*i+3] -= dt*(w*gradN - divv*H[i])*detJt;
            fe[6*i+4] -= dt*(gradN*(j[0]+je*m_pMat->m_penalty)
                         - H[i]*((phiw*spt.m_ca[0] - spt.m_crp[0]/J)/dt)
                         )*detJt;
            fe[6*i+5] -= dt*(gradN*(j[1]+je*m_pMat->m_penalty)
                             - H[i]*((phiw*spt.m_ca[1] - spt.m_crp[1]/J)/dt)
                             )*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::InternalForcesSS(FEGlobalVector& R)
{
    size_t NE = m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInternalForceSS(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements (steady-state)

void FETriphasicDomain::ElementInternalForceSS(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN;
    mat3ds s;
    
    const double* Gr, *Gs, *Gt, *H;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double*	gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // repeat for all integration points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint& spt = *(el.GetMaterialPoint(n)->ExtractData<FESolutesMaterialPoint>());
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
        detJt *= gw[n];
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // get the stress for this integration point
        s = pt.m_s;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solute flux
        vec3d j[2] = {spt.m_j[0],spt.m_j[1]};
        // get the charge number
        int z[2] = {m_pMat->m_pSolute[0]->ChargeNumber(), m_pMat->m_pSolute[1]->ChargeNumber()};
        
        vec3d je = j[0]*z[0] + j[1]*z[1];
        
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN = vec3d(Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i],
                          Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i],
                          Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i]);
            
            // calculate internal force
            vec3d fu = s*gradN;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[6*i  ] -= fu.x*detJt;
            fe[6*i+1] -= fu.y*detJt;
            fe[6*i+2] -= fu.z*detJt;
            fe[6*i+3] -= dt*(w*gradN)*detJt;
            fe[6*i+4] -= dt*(gradN*(j[0]+je*m_pMat->m_penalty))*detJt;
            fe[6*i+5] -= dt*(gradN*(j[1]+je*m_pMat->m_penalty))*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
	// repeat over all solid elements
	size_t NE = m_Elem.size();
    
	#pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);

		// get the lm vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 6;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementTriphasicStiffness(el, ke, bsymm);
		
		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
	// repeat over all solid elements
	size_t NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);

		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndpn = 6;
		int ndof = neln*ndpn;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementTriphasicStiffnessSS(el, ke, bsymm);

		//  get the lm vector
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
bool FETriphasicDomain::ElementTriphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, isol, jsol, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    FETriphasic* pm = m_pMat;
    const int nsol = 2;
    int ndpn = 4+nsol;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        vec3d w = ppt.m_w;
        vec3d gradp = ppt.m_gradp;
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (isol=0; isol<nsol; ++isol)
            z[isol] = pm->m_pSolute[isol]->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        vector< vector<double> > dkdr(spt.m_dkdr);
        vector< vector<double> > dkdJr(spt.m_dkdJr);
        vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
        
        // evaluate the porosity and its derivative
        double phiw = pm->Porosity(mp);
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the osmotic coefficient
        double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
        
        // evaluate the permeability
        mat3ds K = pm->m_pPerm->Permeability(mp);
        tens4dmm dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
        
        vector<mat3ds> dKdc(nsol);
        vector<mat3ds> D(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
        vector<double> D0(nsol);
        vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
        vector<double> dodc(nsol);
        vector<mat3ds> dTdc(nsol);
        vector<mat3ds> ImD(nsol);
        mat3dd I(1);
        
        for (isol=0; isol<nsol; ++isol) {
            // evaluate the permeability derivatives
            dKdc[isol] = pm->m_pPerm->Tangent_Permeability_Concentration(mp,isol);
            
            // evaluate the diffusivity tensor and its derivatives
            D[isol] = pm->m_pSolute[isol]->m_pDiff->Diffusivity(mp);
            dDdE[isol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(mp);
            
            // evaluate the solute free diffusivity
            D0[isol] = pm->m_pSolute[isol]->m_pDiff->Free_Diffusivity(mp);
            
            // evaluate the derivative of the osmotic coefficient
            dodc[isol] = pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            
            // evaluate the stress tangent with concentration
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
            dTdc[isol] = mat3ds(0,0,0,0,0,0);
            
            ImD[isol] = I-D[isol]/D0[isol];
            
            for (jsol=0; jsol<nsol; ++jsol) {
                dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
                dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
            }
        }
        
        // Miscellaneous constants
        double R = pm->m_Rgas;
        double T = pm->m_Tabs;
        double penalty = pm->m_penalty;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds Ke(0,0,0,0,0,0);
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE);
        vector<mat3ds> Gc(nsol);
        vector<mat3ds> dKedc(nsol);
        for (isol=0; isol<nsol; ++isol) {
            Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
            G += dyad1(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
            +(dyad1(I,I) - dyad2(I,I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
            Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
            for (jsol=0; jsol<nsol; ++jsol) {
                Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
                -(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
            }
            Gc[isol] *= R*T/phiw;
        }
        Ke = (Ki + Ke*(R*T/phiw)).inverse();
        tens4d dKedE = (dyad1(Ke,I) - 2*dyad4(Ke,I))*2 - ddot(dyad2(Ke,Ke),G);
        for (isol=0; isol<nsol; ++isol)
            dKedc[isol] = -(Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke).sym();
        
        // calculate all the matrices
        vec3d vtmp,gp,qpu;
        vector<vec3d> gc(nsol),qcu(nsol),wc(nsol),jce(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        mat3d wu, jue;
        vector<mat3d> ju(nsol);
        vector< vector<double> > qcc(nsol, vector<double>(nsol));
        double sum;
        mat3ds De;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                ke[6*i  ][6*j  ] += Kuu[0][0]; ke[6*i  ][6*j+1] += Kuu[0][1]; ke[6*i  ][6*j+2] += Kuu[0][2];
                ke[6*i+1][6*j  ] += Kuu[1][0]; ke[6*i+1][6*j+1] += Kuu[1][1]; ke[6*i+1][6*j+2] += Kuu[1][2];
                ke[6*i+2][6*j  ] += Kuu[2][0]; ke[6*i+2][6*j+1] += Kuu[2][1]; ke[6*i+2][6*j+2] += Kuu[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradN[j]);
                for (isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                qpu = -gradN[j]*(1.0/dt);
                vtmp = (wu.transpose()*gradN[i] + qpu*H[i])*(detJ*dt);
                ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
                ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
                ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
                ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
                ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+3][ndpn*j+3] += (- gradN[i]*(Ke*gradN[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero();
                De.zero();
                for (isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
                       +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                    qcu[isol] = qpu*(c[isol]*(kappa[isol]+J*phiw*dkdJ[isol]));
                }
                
                for (isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i]
                            + qcu[isol]*H[i])*(detJ*dt);
                    ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
                    ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
                    ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
                                                              (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                               +De*penalty)
                                                              *(Ke*gradN[j])
                                                              ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*detJ;
                    ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
                    ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
                    ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-H[j])
                    -Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
                    ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*c[isol]*dkdc[isol][jsol];
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                            
                            qcc[isol][jsol] = -H[j]*phiw/dt*(c[isol]*dkdc[isol][jsol] + kappa[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty)
                                                             + H[i]*(qcc[isol][jsol]))*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    if (bsymm) {
        for (i=0; i<ndpn*neln; ++i)
            for (j=i+1; j<ndpn*neln; ++j) {
                tmp = 0.5*(ke[i][j]+ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for steady-state response (zero solid velocity, zero time derivative of
//! solute concentration)
//!
bool FETriphasicDomain::ElementTriphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, isol, jsol, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Gradient of shape functions
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    FETriphasic* pm = m_pMat;
    const int nsol = 2;
    int ndpn = 4+nsol;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint&  ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        FESolutesMaterialPoint&  spt = *(mp.ExtractData<FESolutesMaterialPoint >());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n)*gw[n];
        
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // calculate global gradient of shape functions
        for (i=0; i<neln; ++i)
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds C = m_pMat->Tangent(mp);
        
        // next we get the determinant
        double J = ept.m_J;
        
        // get the fluid flux and pressure gradient
        vec3d w = ppt.m_w;
        vec3d gradp = ppt.m_gradp;
        
        vector<double> c(spt.m_c);
        vector<vec3d> gradc(spt.m_gradc);
        vector<int> z(nsol);
        
        vector<double> kappa(spt.m_k);
        
        // get the charge number
        for (isol=0; isol<nsol; ++isol)
            z[isol] = pm->m_pSolute[isol]->ChargeNumber();
        
        vector<double> dkdJ(spt.m_dkdJ);
        vector< vector<double> > dkdc(spt.m_dkdc);
        vector< vector<double> > dkdr(spt.m_dkdr);
        vector< vector<double> > dkdJr(spt.m_dkdJr);
        vector< vector< vector<double> > > dkdrc(spt.m_dkdrc);
        
        // evaluate the porosity and its derivative
        double phiw = pm->Porosity(mp);
        double phis = 1. - phiw;
        double dpdJ = phis/J;
        
        // evaluate the osmotic coefficient
        double osmc = pm->m_pOsmC->OsmoticCoefficient(mp);
        
        // evaluate the permeability
        mat3ds K = pm->m_pPerm->Permeability(mp);
        tens4dmm dKdE = pm->m_pPerm->Tangent_Permeability_Strain(mp);
        
        vector<mat3ds> dKdc(nsol);
        vector<mat3ds> D(nsol);
        vector<tens4dmm> dDdE(nsol);
        vector< vector<mat3ds> > dDdc(nsol, vector<mat3ds>(nsol));
        vector<double> D0(nsol);
        vector< vector<double> > dD0dc(nsol, vector<double>(nsol));
        vector<double> dodc(nsol);
        vector<mat3ds> dTdc(nsol);
        vector<mat3ds> ImD(nsol);
        mat3dd I(1);
        
        for (isol=0; isol<nsol; ++isol) {
            // evaluate the permeability derivatives
            dKdc[isol] = pm->m_pPerm->Tangent_Permeability_Concentration(mp,isol);
            
            // evaluate the diffusivity tensor and its derivatives
            D[isol] = pm->m_pSolute[isol]->m_pDiff->Diffusivity(mp);
            dDdE[isol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Strain(mp);
            
            // evaluate the solute free diffusivity
            D0[isol] = pm->m_pSolute[isol]->m_pDiff->Free_Diffusivity(mp);
            
            // evaluate the derivative of the osmotic coefficient
            dodc[isol] = pm->m_pOsmC->Tangent_OsmoticCoefficient_Concentration(mp,isol);
            
            // evaluate the stress tangent with concentration
            //			dTdc[isol] = pm->GetSolid()->Tangent_Concentration(mp,isol);
            dTdc[isol] = mat3ds(0,0,0,0,0,0);
            
            ImD[isol] = I-D[isol]/D0[isol];
            
            for (jsol=0; jsol<nsol; ++jsol) {
                dDdc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Diffusivity_Concentration(mp,jsol);
                dD0dc[isol][jsol] = pm->m_pSolute[isol]->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp,jsol);
            }
        }
        
        // Miscellaneous constants
        double R = pm->m_Rgas;
        double T = pm->m_Tabs;
        double penalty = pm->m_penalty;
        
        // evaluate the effective permeability and its derivatives
        mat3ds Ki = K.inverse();
        mat3ds Ke(0,0,0,0,0,0);
        tens4d G = (dyad1(Ki,I) - dyad4(Ki,I)*2)*2 - ddot(dyad2(Ki,Ki),dKdE);
        vector<mat3ds> Gc(nsol);
        vector<mat3ds> dKedc(nsol);
        for (isol=0; isol<nsol; ++isol) {
            Ke += ImD[isol]*(kappa[isol]*c[isol]/D0[isol]);
            G += dyad1(ImD[isol],I)*(R*T*c[isol]*J/D0[isol]/phiw*(dkdJ[isol]-kappa[isol]/phiw*dpdJ))
            +(dyad1(I,I) - dyad2(I,I)*2 - dDdE[isol]/D0[isol])*(R*T*kappa[isol]*c[isol]/phiw/D0[isol]);
            Gc[isol] = ImD[isol]*(kappa[isol]/D0[isol]);
            for (jsol=0; jsol<nsol; ++jsol) {
                Gc[isol] += ImD[jsol]*(c[jsol]/D0[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol]))
                -(dDdc[jsol][isol]-D[jsol]*(dD0dc[jsol][isol]/D0[jsol])*(kappa[jsol]*c[jsol]/SQR(D0[jsol])));
            }
            Gc[isol] *= R*T/phiw;
        }
        Ke = (Ki + Ke*(R*T/phiw)).inverse();
        tens4d dKedE = (dyad1(Ke,I) - 2*dyad4(Ke,I))*2 - ddot(dyad2(Ke,Ke),G);
        for (isol=0; isol<nsol; ++isol)
            dKedc[isol] = -(Ke*(-Ki*dKdc[isol]*Ki + Gc[isol])*Ke).sym();
        
        // calculate all the matrices
        vec3d vtmp,gp;
        vector<vec3d> gc(nsol),qcu(nsol),wc(nsol),jce(nsol);
        vector< vector<vec3d> > jc(nsol, vector<vec3d>(nsol));
        mat3d wu, jue;
        vector<mat3d> ju(nsol);
        double sum;
        mat3ds De;
        for (i=0; i<neln; ++i)
        {
            for (j=0; j<neln; ++j)
            {
                // Kuu matrix
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], C, gradN[j]))*detJ;
                ke[6*i  ][6*j  ] += Kuu[0][0]; ke[6*i  ][6*j+1] += Kuu[0][1]; ke[6*i  ][6*j+2] += Kuu[0][2];
                ke[6*i+1][6*j  ] += Kuu[1][0]; ke[6*i+1][6*j+1] += Kuu[1][1]; ke[6*i+1][6*j+2] += Kuu[1][2];
                ke[6*i+2][6*j  ] += Kuu[2][0]; ke[6*i+2][6*j+1] += Kuu[2][1]; ke[6*i+2][6*j+2] += Kuu[2][2];
                
                // calculate the kpu matrix
                gp = vec3d(0,0,0);
                for (isol=0; isol<nsol; ++isol) gp += (D[isol]*gradc[isol])*(kappa[isol]/D0[isol]);
                gp = gradp+gp*(R*T);
                wu = vdotTdotv(-gp, dKedE, gradN[j]);
                for (isol=0; isol<nsol; ++isol) {
                    wu += (((Ke*(D[isol]*gradc[isol])) & gradN[j])*(J*dkdJ[isol] - kappa[isol])
                           +Ke*(2*kappa[isol]*(gradN[j]*(D[isol]*gradc[isol]))))*(-R*T/D0[isol])
                    + (Ke*vdotTdotv(gradc[isol], dDdE[isol], gradN[j]))*(-kappa[isol]*R*T/D0[isol]);
                }
                vtmp = (wu.transpose()*gradN[i])*(detJ*dt);
                ke[ndpn*i+3][ndpn*j  ] += vtmp.x;
                ke[ndpn*i+3][ndpn*j+1] += vtmp.y;
                ke[ndpn*i+3][ndpn*j+2] += vtmp.z;
                
                // calculate the kup matrix
                vtmp = -gradN[i]*H[j]*detJ;
                ke[ndpn*i  ][ndpn*j+3] += vtmp.x;
                ke[ndpn*i+1][ndpn*j+3] += vtmp.y;
                ke[ndpn*i+2][ndpn*j+3] += vtmp.z;
                
                // calculate the kpp matrix
                ke[ndpn*i+3][ndpn*j+3] += (- gradN[i]*(Ke*gradN[j]))*(detJ*dt);
                
                // calculate kcu matrix data
                jue.zero();
                De.zero();
                for (isol=0; isol<nsol; ++isol) {
                    gc[isol] = -gradc[isol]*phiw + w*c[isol]/D0[isol];
                    ju[isol] = ((D[isol]*gc[isol]) & gradN[j])*(J*dkdJ[isol])
                    + vdotTdotv(gc[isol], dDdE[isol], gradN[j])*kappa[isol]
                    + (((D[isol]*gradc[isol]) & gradN[j])*(-phis)
                       +(D[isol]*((gradN[j]*w)*2) - ((D[isol]*w) & gradN[j]))*c[isol]/D0[isol]
                       )*kappa[isol]
                    +D[isol]*wu*(kappa[isol]*c[isol]/D0[isol]);
                    jue += ju[isol]*z[isol];
                    De += D[isol]*(z[isol]*kappa[isol]*c[isol]/D0[isol]);
                }
                
                for (isol=0; isol<nsol; ++isol) {
                    
                    // calculate the kcu matrix
                    vtmp = ((ju[isol]+jue*penalty).transpose()*gradN[i]
                            + qcu[isol]*H[i])*(detJ*dt);
                    ke[ndpn*i+4+isol][ndpn*j  ] += vtmp.x;
                    ke[ndpn*i+4+isol][ndpn*j+1] += vtmp.y;
                    ke[ndpn*i+4+isol][ndpn*j+2] += vtmp.z;
                    
                    // calculate the kcp matrix
                    ke[ndpn*i+4+isol][ndpn*j+3] -= (gradN[i]*(
                                                              (D[isol]*(kappa[isol]*c[isol]/D0[isol])
                                                               +De*penalty)
                                                              *(Ke*gradN[j])
                                                              ))*(detJ*dt);
                    
                    // calculate the kuc matrix
                    sum = 0;
                    for (jsol=0; jsol<nsol; ++jsol)
                        sum += c[jsol]*(dodc[isol]*kappa[jsol]+osmc*dkdc[jsol][isol]);
                    vtmp = (dTdc[isol]*gradN[i] - gradN[i]*(R*T*(osmc*kappa[isol]+sum)))*H[j]*detJ;
                    ke[ndpn*i  ][ndpn*j+4+isol] += vtmp.x;
                    ke[ndpn*i+1][ndpn*j+4+isol] += vtmp.y;
                    ke[ndpn*i+2][ndpn*j+4+isol] += vtmp.z;
                    
                    // calculate the kpc matrix
                    vtmp = vec3d(0,0,0);
                    for (jsol=0; jsol<nsol; ++jsol)
                        vtmp += (D[jsol]*(dkdc[jsol][isol]-kappa[jsol]/D0[jsol]*dD0dc[jsol][isol])
                                 +dDdc[jsol][isol]*kappa[jsol])/D0[jsol]*gradc[jsol];
                    wc[isol] = (dKedc[isol]*gp)*(-H[j])
                    -Ke*((D[isol]*gradN[j])*(kappa[isol]/D0[isol])+vtmp*H[j])*(R*T);
                    ke[ndpn*i+3][ndpn*j+4+isol] += (gradN[i]*wc[isol])*(detJ*dt);
                    
                }
                
                // calculate data for the kcc matrix
                jce.assign(nsol, vec3d(0,0,0));
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        if (jsol != isol) {
                            jc[isol][jsol] =
                            ((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        else {
                            jc[isol][jsol] = (D[isol]*(gradN[j]*(-phiw)+w*(H[j]/D0[isol])))*kappa[isol]
                            +((D[isol]*dkdc[isol][jsol]+dDdc[isol][jsol]*kappa[isol])*gc[isol])*H[j]
                            +(D[isol]*(w*(-H[j]*dD0dc[isol][jsol]/D0[isol])+wc[jsol]))*(kappa[isol]*c[isol]/D0[isol]);
                        }
                        jce[jsol] += jc[isol][jsol]*z[isol];
                    }
                }
                
                // calculate the kcc matrix
                for (isol=0; isol<nsol; ++isol) {
                    for (jsol=0; jsol<nsol; ++jsol) {
                        ke[ndpn*i+4+isol][ndpn*j+4+jsol] += (gradN[i]*(jc[isol][jsol]+jce[jsol]*penalty))*(detJ*dt);
                    }
                }
            }
        }
    }
    
    // Enforce symmetry by averaging top-right and bottom-left corners of stiffness matrix
    if (bsymm) {
        for (i=0; i<ndpn*neln; ++i)
            for (j=i+1; j<ndpn*neln; ++j) {
                tmp = 0.5*(ke[i][j]+ke[j][i]);
                ke[i][j] = ke[j][i] = tmp;
            }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::Update(const FETimeInfo& tp)
{
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

	// if we encountered an error, throw an exception
	if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
void FETriphasicDomain::UpdateElementStress(int iel)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	// get the number of nodes
	int neln = el.Nodes();
		
	// get the biphasic-solute material
	int id0 = m_dofC + m_pMat->m_pSolute[0]->GetSoluteDOF();
	int id1 = m_dofC + m_pMat->m_pSolute[1]->GetSoluteDOF();
		
	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES], ct[2][FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
		pn[j] = mesh.Node(el.m_node[j]).get(m_dofP);
		ct[0][j] = mesh.Node(el.m_node[j]).get(id0);
		ct[1][j] = mesh.Node(el.m_node[j]).get(id1);
	}
		
	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			
		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		mp.m_r0 = el.Evaluate(r0, n);
		mp.m_rt = el.Evaluate(rt, n);
			
		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
			
		// solute-poroelastic data
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
		FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
			
		// evaluate fluid pressure at gauss-point
		ppt.m_p = el.Evaluate(pn, n);
			
		// calculate the gradient of p at gauss-point
		ppt.m_gradp = gradient(el, pn, n);
			
		// evaluate effective solute concentration at gauss-point
		spt.m_c[0] = el.Evaluate(ct[0], n);
		spt.m_c[1] = el.Evaluate(ct[1], n);
			
		// calculate the gradient of c at gauss-point
		spt.m_gradc[0] = gradient(el, ct[0], n);
		spt.m_gradc[1] = gradient(el, ct[1], n);
			
		// for biphasic-solute materials also update the porosity, fluid and solute fluxes
		// and evaluate the actual fluid pressure and solute concentration
		ppt.m_w = m_pMat->FluidFlux(mp);
		spt.m_psi = m_pMat->ElectricPotential(mp);
		spt.m_ca[0] = m_pMat->Concentration(mp,0);
		spt.m_ca[1] = m_pMat->Concentration(mp,1);
		ppt.m_pa = m_pMat->Pressure(mp);
		spt.m_j[0] = m_pMat->SoluteFlux(mp,0);
		spt.m_j[1] = m_pMat->SoluteFlux(mp,1);
		spt.m_cF = m_pMat->FixedChargeDensity(mp);
		spt.m_Ie = m_pMat->CurrentDensity(mp);
        m_pMat->PartitionCoefficientFunctions(mp, spt.m_k, spt.m_dkdJ, spt.m_dkdc);
			
        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the solid stress at this material point
        ppt.m_ss = m_pMat->GetElasticMaterial()->Stress(mp);
        
		pt.m_s = m_pMat->Stress(mp);
	}
}
