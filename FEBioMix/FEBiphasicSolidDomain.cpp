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
#include "FEBiphasicSolidDomain.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include <FECore/FEDataExport.h>
#include <FECore/FEModel.h>
#include <FEBioMech/FEBioMech.h>
#include <FECore/FELinearSystem.h>
#include "FEBioMix.h"

//-----------------------------------------------------------------------------
FEBiphasicSolidDomain::FEBiphasicSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEBiphasicDomain(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofR(pfem), m_dof(pfem)
{
	EXPORT_DATA(PLT_FLOAT, FMT_NODE, &m_nodePressure, "NPR fluid pressure");

	if (pfem)
	{
		m_varU = pfem->GetDOFS().GetVariableIndex(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT)); assert(m_varU >= 0);
		m_varP = pfem->GetDOFS().GetVariableIndex(FEBioMix::GetVariableName(FEBioMix::FLUID_PRESSURE)); assert(m_varP >= 0);

		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
		m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
	}
}

//-----------------------------------------------------------------------------
//! get the total dof list
const FEDofList& FEBiphasicSolidDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Serialize(DumpStream& ar)
{
    FESolidDomain::Serialize(ar);

    if (ar.IsShallow() == false)
    {
        ar & m_nodePressure;
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FEBiphasic*>(pmat);
	assert(m_pMat);
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEBiphasicSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
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
            FENode& node = m.Node(el.m_node[i]);
			x0[i] = node.m_r0;
			xt[i] = node.m_rt;
            if (el.m_bitfc.size()>0 && el.m_bitfc[i] && node.m_ID[m_dofQ] != -1)
                pn[i] = node.get(m_dofQ);
            else
                pn[i] = node.get(m_dofP);
        }

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);
            p = el.Evaluate(pn, j);
            
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint& pb = *mp.ExtractData<FEBiphasicMaterialPoint>();
			mp.m_r0 = r0;
			mp.m_rt = rt;

			pt.m_J = defgrad(el, pt.m_F, j);

            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, j);
            pb.m_gradpp = pb.m_gradp;
            pb.m_phi0p = pb.m_phi0t;
            
            mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
bool FEBiphasicSolidDomain::Init()
{
	// initialize base class
	if (FESolidDomain::Init() == false) return false;
    
    // initialize body forces
	FEModel& fem = *GetFEModel();
	m_pMat->m_bf.clear();
    for (int j=0; j<fem.ModelLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.ModelLoad(j));
        if (pbf) m_pMat->m_bf.push_back(pbf);
    }

	// allocate nodal pressures
	m_nodePressure.resize(Nodes(), 0.0);
    
	return true;
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Activate()
{
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
		}
	}

    // Activate dof_P, except when a biphasic solid is connected to the
    // back of a shell element, in which case activate dof_Q for those nodes.
    FEMesh& m = *GetMesh();
    for (int i=0; i<Elements(); ++i) {
        FESolidElement& el = m_Elem[i];
        int neln = el.Nodes();
        for (int j=0; j<neln; ++j)
        {
            FENode& node = m.Node(el.m_node[j]);
            if (el.m_bitfc.size()>0 && el.m_bitfc[j])
                node.set_active(m_dofQ);
            else
                node.set_active(m_dofP);
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEBiphasicSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

	// number of nodes for velocity interpolation
	int neln_d = el.ShapeFunctions(degree_d);

	// number of nodes for pressure interpolation
	int neln_p = el.ShapeFunctions(degree_p);

	// allocate lm
	lm.resize(neln_d*4 + 3*neln_d);

	// displacement dofs
	for (int i=0; i<neln_d; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

        // first the displacement dofs
        lm[4*i  ] = id[m_dofU[0]];
        lm[4*i+1] = id[m_dofU[1]];
        lm[4*i+2] = id[m_dofU[2]];

		// now the pressure dofs
		lm[4*i + 3] = id[m_dofP];

		// rigid rotational dofs
		lm[4 * neln_d + 3 * i    ] = id[m_dofR[0]];
		lm[4 * neln_d + 3 * i + 1] = id[m_dofR[1]];
		lm[4 * neln_d + 3 * i + 2] = id[m_dofR[2]];
	}

    // substitute interface dofs for solid-shell interfaces
	FESolidElement& sel = static_cast<FESolidElement&>(el);
	for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the back-face displacement dofs
            lm[4*i  ] = id[m_dofSU[0]];
            lm[4*i+1] = id[m_dofSU[1]];
            lm[4*i+2] = id[m_dofSU[2]];
            
            // now the pressure dof (if the shell has it)
            if (id[m_dofQ] > -1) lm[4*i + 3] = id[m_dofQ];
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Reset()
{
	// reset base class data
	FESolidDomain::Reset();

	// initialize all element data
	ForEachMaterialPoint([=](FEMaterialPoint& mp) {
		FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());

		// initialize referential solid volume fraction
		pt.m_phi0 = pt.m_phi0t = m_pMat->m_phi0(mp);
	});
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::InternalForces(FEGlobalVector& R)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

	int NE = (int)m_Elem.size();
	#pragma omp parallel for shared (NE)
	for (int i=0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> lm;
		
		// get the element
		FESolidElement& el = m_Elem[i];

		int nel_d = el.ShapeFunctions(degree_d);
		int nel_p = el.ShapeFunctions(degree_p);

		// get the element force vector and initialize it to zero
		int ndof = 4*nel_d;
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

void FEBiphasicSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3];

	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    int nint = el.GaussPoints();
	int nel_d = el.ShapeFunctions(degree_d);
	int nel_p = el.ShapeFunctions(degree_p);

    double*	gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // repeat for all integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
		// calculate the jacobian
		double Jw = invjact(el, Ji, n)*gw[n];

        // get the stress vector for this integration point
        mat3ds s = pt.m_s;
        
        double* Gr = el.Gr(n);
        double* Gs = el.Gs(n);
        double* Gt = el.Gt(n);
        
        double* H = el.H(n);

		// --- stress contribution
        for (int i=0; i<nel_d; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            vec3d gradN(Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i],
                        Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i],
                        Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i]);

            // calculate internal force
            vec3d fu = s*gradN;
            
            // the '-' sign is so that the internal forces get subtracted
            // from the global residual vector
            fe[4*i  ] -= fu.x*Jw;
            fe[4*i+1] -= fu.y*Jw;
            fe[4*i+2] -= fu.z*Jw;
        }

		// --- pressure contribution

		// next we get the determinant
		double Jp = bpt.m_Jp;
		double J = pt.m_J;

		// and then finally
		double divv = ((J - Jp) / dt) / J;

		// get the flux
		vec3d& w = bpt.m_w;

		// get the solvent supply
		double phiwhat = m_pMat->SolventSupply(mp);

		// pressure shape functions
		double* Hp  = el.H(degree_p, n);
		double* Gpr = el.Gr(degree_p, n);
		double* Gps = el.Gs(degree_p, n);
		double* Gpt = el.Gt(degree_p, n);

		for (int i = 0; i<nel_p; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			vec3d gradHp = vec3d(Ji[0][0] * Gpr[i] + Ji[1][0] * Gps[i] + Ji[2][0] * Gpt[i],
				                 Ji[0][1] * Gpr[i] + Ji[1][1] * Gps[i] + Ji[2][1] * Gpt[i],
				                 Ji[0][2] * Gpr[i] + Ji[1][2] * Gps[i] + Ji[2][2] * Gpt[i]);

			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[4*i + 3] -= dt*(w*gradHp + (phiwhat - divv)*Hp[i])*Jw;
		}
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::InternalForcesSS(FEGlobalVector& R)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 4*el.Nodes();
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

void FEBiphasicSolidDomain::ElementInternalForceSS(FESolidElement& el, vector<double>& fe)
{
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN, GradN;
   
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    int nint = el.GaussPoints();
	int nel_d = el.ShapeFunctions(degree_d);
	int nel_p = el.ShapeFunctions(degree_p);

    double*	gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // repeat for all integration points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        FEBiphasicMaterialPoint& bpt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        mat3d s = pt.m_s;
        
        double* Gr = el.Gr(n);
        double* Gs = el.Gs(n);
        double* Gt = el.Gt(n);
        
        double* H = el.H(n);
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solvent supply
        double phiwhat = m_pMat->SolventSupply(mp);
        
        for (int i=0; i<nel_d; ++i)
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
            fe[4*i  ] -= fu.x*detJt;
            fe[4*i+1] -= fu.y*detJt;
            fe[4*i+2] -= fu.z*detJt;
        }

		// --- pressure contribution

		double* Gpr = el.Gr(degree_p, n);
		double* Gps = el.Gs(degree_p, n);
		double* Gpt = el.Gt(degree_p, n);

		double* Hp = el.H(degree_p, n);

		for (int i = 0; i<nel_p; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			vec3d gradH(Ji[0][0] * Gpr[i] + Ji[1][0] * Gps[i] + Ji[2][0] * Gpt[i],
						Ji[0][1] * Gpr[i] + Ji[1][1] * Gps[i] + Ji[2][1] * Gpt[i],
						Ji[0][2] * Gpr[i] + Ji[1][2] * Gps[i] + Ji[2][2] * Gpt[i]);
			
			fe[4*i + 3] -= dt*(w*gradH + phiwhat*Hp[i])*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrix(FELinearSystem& LS, bool bsymm)
{
	// repeat over all solid elements
	int NE = (int)m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);
		int ndof = el.Nodes()*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffness(el, ke, bsymm);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrixSS(FELinearSystem& LS, bool bsymm)
{
	// repeat over all solid elements
	int NE = (int)m_Elem.size();

	#pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);
		int ndof = el.Nodes()*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffnessSS(el, ke, bsymm);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
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
bool FEBiphasicSolidDomain::ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    int nint = el.GaussPoints();
	int nel_d = el.ShapeFunctions(degree_d);
	int nel_p = el.ShapeFunctions(degree_p);

    // jacobian
    double Ji[3][3];
    
    // Bp-matrix
    vector<vec3d> gradNu(FEElement::MAX_NODES), gradNp(FEElement::MAX_NODES);
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    double tau = m_pMat->m_tau;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate jacobian
        double detJ = invjact(el, Ji, n);
        
        // contravariant basis vectors in spatial frame
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
		// displacement shape functions
        double* Hu = el.H(n);
        double* Gur = el.Gr(n);
        double* Gus = el.Gs(n);
        double* Gut = el.Gt(n);

		// pressure shape functions
		double* Hp  = el.H (degree_p, n);
		double* Gpr = el.Gr(degree_p, n);
		double* Gps = el.Gs(degree_p, n);
		double* Gpt = el.Gt(degree_p, n);

        for (int i=0; i<nel_d; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradNu[i] = g1*Gur[i] + g2*Gus[i] + g3*Gut[i];
        }

		for (int i = 0; i<nel_p; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			gradNp[i] = g1*Gpr[i] + g2*Gps[i] + g3*Gpt[i];
		}

        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp + (pt.m_gradp - pt.m_gradpp)*(tau/dt);
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        // evaluate the solvent supply and its derivatives
        double phiwhat = 0;
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        if (m_pMat->GetSolventSupply()) {
            phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // Miscellaneous constants
        mat3dd I(1);
        
        // Kuu matrix
        double Jw = detJ*gw[n];
        for (int i=0; i<nel_d; ++i)
            for (int j=0; j<nel_d; ++j)
            {
                mat3d Kuu = (mat3dd(gradNu[i]*(s*gradNu[j])) + vdotTdotv(gradNu[i], c, gradNu[j]))*Jw;
                
				ke.add(4 * i, 4 * j, Kuu);
            }
        
        // calculate the kpp matrix
        for (int i=0; i<nel_p; ++i)
            for (int j=0; j<nel_p; ++j)
            {
                ke[4*i+3][4*j+3] += (Hp[i]*Hp[j]*Phip - gradNp[i]*(K*gradNp[j])*(1+tau/dt))*(dt*Jw);
            }
        
        if (!bsymm) {
            // calculate the kup matrix
            for (int i=0; i<nel_d; ++i) {
                for (int j=0; j<nel_p; ++j)
                {
                    ke[4*i  ][4*j+3] -= Jw*gradNu[i].x*Hp[j];
                    ke[4*i+1][4*j+3] -= Jw*gradNu[i].y*Hp[j];
                    ke[4*i+2][4*j+3] -= Jw*gradNu[i].z*Hp[j];
                }
            }
            
            // calculate the kpu matrix
            mat3ds Q = Phie*ept.m_J + mat3dd(phiwhat - 1./dt);
            for (int i=0; i<nel_p; ++i) {
                for (int j=0; j<nel_d; ++j)
                {
                    vec3d vt = ((vdotTdotv(-gradNp[i], dKdE, gradNu[j])*gradp)
                                +(Q*gradNu[j])*Hp[i])*(Jw*dt);
                    ke[4*i+3][4*j  ] += vt.x;
                    ke[4*i+3][4*j+1] += vt.y;
                    ke[4*i+3][4*j+2] += vt.z;
                }
            }
            
        } else {
            // calculate the kup matrix and let kpu be its symmetric part
            for (int i=0; i<nel_d; ++i) {
                for (int j=0; j<nel_p; ++j)
                {
                    ke[4*i  ][4*j+3] -= Jw*gradNu[i].x*Hp[j];
                    ke[4*i+1][4*j+3] -= Jw*gradNu[i].y*Hp[j];
                    ke[4*i+2][4*j+3] -= Jw*gradNu[i].z*Hp[j];
                    
                    ke[4*j+3][4*i  ] -= Jw*gradNu[i].x*Hp[j];
                    ke[4*j+3][4*i+1] -= Jw*gradNu[i].y*Hp[j];
                    ke[4*j+3][4*i+2] -= Jw*gradNu[i].z*Hp[j];
                }
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//! for the steady-state response (zero solid velocity)
//!
bool FEBiphasicSolidDomain::ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm)
{
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

    int nint = el.GaussPoints();
	int nel_d = el.ShapeFunctions(degree_d);
	int nel_p = el.ShapeFunctions(degree_p);

    // jacobian
    double Ji[3][3];
    
    // Bp-matrix
    vector<vec3d> gradNu(nel_d), gradNp(nel_p);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate jacobian
        double detJ = invjact(el, Ji, n);
        
        // contravariant basis vectors in spatial frame
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        double* Hu = el.H(n);
		double* Hp = el.H(degree_p, n);

        double* Gur = el.Gr(n);
        double* Gus = el.Gs(n);
        double* Gut = el.Gt(n);

		double* Gpr = el.Gr(degree_p, n);
		double* Gps = el.Gs(degree_p, n);
		double* Gpt = el.Gt(degree_p, n);

        for (int i=0; i<nel_d; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradNu[i] = g1*Gur[i] + g2*Gus[i] + g3*Gut[i];
        }

		for (int i = 0; i<nel_p; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			gradNp[i] = g1*Gpr[i] + g2*Gps[i] + g3*Gpt[i];
		}

        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        // evaluate the solvent supply and its derivatives
        double phiwhat = 0;
        mat3ds Phie; Phie.zero();
        double Phip = 0;
        if (m_pMat->GetSolventSupply()) {
            phiwhat = m_pMat->GetSolventSupply()->Supply(mp);
            Phie = m_pMat->GetSolventSupply()->Tangent_Supply_Strain(mp);
            Phip = m_pMat->GetSolventSupply()->Tangent_Supply_Pressure(mp);
        }
        
        // Miscellaneous constants
        mat3dd I(1);
        
        // Kuu matrix
        tmp = detJ*gw[n];
        for (int i=0; i<nel_d; ++i)
            for (int j=0; j<nel_d; ++j)
            {
                mat3d Kuu = (mat3dd(gradNu[i]*(s*gradNu[j])) + vdotTdotv(gradNu[i], c, gradNu[j]))*tmp;
				ke.add(4 * i, 4 * j, Kuu);
            }
        
        // calculate the kpp matrix
        tmp = detJ*gw[n]*dt;
        for (int i=0; i<nel_p; ++i)
            for (int j=0; j<nel_p; ++j)
            {
                ke[4*i+3][4*j+3] += (Hp[i]*Hp[j]*Phip - gradNp[i]*(K*gradNp[j]))*tmp;
            }
        
        if (!bsymm) {
            // calculate the kup matrix
            for (int i=0; i<nel_d; ++i) {
                for (int j=0; j<nel_p; ++j)
                {
                    tmp = detJ*gw[n]*Hp[j];
                    ke[4*i    ][4*j+3] -= tmp*gradNu[i].x;
                    ke[4*i + 1][4*j+3] -= tmp*gradNu[i].y;
                    ke[4*i + 2][4*j+3] -= tmp*gradNu[i].z;
                }
            }
            
            // calculate the kpu matrix
            //			tmp = detJ*gw[n];
            tmp = detJ*gw[n]*dt;
            for (int i=0; i<nel_p; ++i) {
                for (int j=0; j<nel_d; ++j)
                {
                    vec3d vt = ((vdotTdotv(-gradp, dKdE, gradNu[j])*(gradNp[i]))
                                +(mat3dd(phiwhat) + Phie*ept.m_J)*gradNu[j]*Hp[i])*tmp;
                    ke[4*i+3][4*j  ] += vt.x;
                    ke[4*i+3][4*j+1] += vt.y;
                    ke[4*i+3][4*j+2] += vt.z;
                }
            }
            
        } else {
            // calculate the kup matrix and let kpu be its symmetric part
            tmp = detJ*gw[n];
            for (int i=0; i<nel_d; ++i) {
                for (int j=0; j<nel_p; ++j)
                {
                    ke[4*i  ][4*j+3] -= tmp*Hu[j]*gradNp[i].x;
                    ke[4*i+1][4*j+3] -= tmp*Hu[j]*gradNp[i].y;
                    ke[4*i+2][4*j+3] -= tmp*Hu[j]*gradNp[i].z;
                    
                    ke[4*j+3][4*i  ] -= tmp*Hu[j]*gradNp[i].x;
                    ke[4*j+3][4*i+1] -= tmp*Hu[j]*gradNp[i].y;
                    ke[4*j+3][4*i+2] -= tmp*Hu[j]*gradNp[i].z;
                }
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Update(const FETimeInfo& tp)
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

    if (berr) throw NegativeJacobianDetected();

	// also update the nodal pressures
	UpdateNodalPressures();
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::UpdateElementStress(int iel)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
   // extract the elastic component
    FEElasticMaterial* pme = m_pMat->GetElasticMaterial();

	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	DOFS& dofs = GetFEModel()->GetDOFS();
	int degree_d = dofs.GetVariableInterpolationOrder(m_varU);
	int degree_p = dofs.GetVariableInterpolationOrder(m_varP);

	// get the number of nodes
	int nel_d = el.ShapeFunctions(degree_d);
	int nel_p = el.ShapeFunctions(degree_p);

	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES];
	for (int j=0; j<nel_d; ++j)
	{
        FENode& node = mesh.Node(el.m_node[j]);
		r0[j] = node.m_r0;
		rt[j] = node.m_rt;
	}

	for (int j = 0; j<nel_p; ++j)
	{
		FENode& node = mesh.Node(el.m_node[j]);
		if (el.m_bitfc.size()>0 && el.m_bitfc[j] && node.m_ID[m_dofQ] != -1)
			pn[j] = node.get(m_dofQ);
		else
			pn[j] = node.get(m_dofP);
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
        mat3d Fp;
        defgradp(el, Fp, n);
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (pt.m_F - Fp)*Fi / dt;

		// poroelasticity data
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
		// evaluate fluid pressure at gauss-point
		ppt.m_p = el.Evaluate(degree_p, pn, n);
			
		// calculate the gradient of p at gauss-point
		ppt.m_gradp = gradient(el, degree_p, pn, n);
			
		// for biphasic materials also update the fluid flux
		ppt.m_w = FluidFlux(mp);
		ppt.m_pa = m_pMat->Pressure(mp);
			
        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, GetFEModel()->GetTime());
        
        // calculate the solid stress at this material point
        ppt.m_ss = m_pMat->GetElasticMaterial()->Stress(mp);
        
		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
	FEBodyForce* bf = &BF;
	LoadVector(R, m_dofU, [=](FEMaterialPoint& mp, int node_a, vector<double>& fa) {

		// get true solid and fluid densities
		double rhoTs = m_pMat->SolidDensity(mp);
		double rhoTw = m_pMat->FluidDensity();

		// Jacobian
		double detJ = mp.m_Jt;

		// get the force
		vec3d b = bf->force(mp);

		// evaluate apparent solid and fluid densities and mixture density
		double phiw = m_pMat->Porosity(mp);
		double rhos = (1 - phiw)*rhoTs;
		double rhow = phiw*rhoTw;
		double rho = rhos + rhow;

		double* H = mp.m_shape;

		fa[0] = -H[node_a] * rho*b.x*detJ;
		fa[1] = -H[node_a] * rho*b.y*detJ;
		fa[2] = -H[node_a] * rho*b.z*detJ;
	});
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(GetMaterial()); assert(pmb);
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);
        int neln = el.Nodes();
        int ndof = 4*neln;
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // TODO: the problem here is that the LM array that is returned by the UnpackLM
        // function does not give the equation numbers in the right order. For this reason we
        // have to create a new lm array and place the equation numbers in the right order.
        // What we really ought to do is fix the UnpackLM function so that it returns
        // the LM vector in the right order for poroelastic elements.
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);
        
        // assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEBiphasicSolidDomain::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    int neln = el.Nodes();
    
    // get true solid and fluid densities
    double rhoTw = m_pMat->FluidDensity();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // jacobian
    double detJt, Ji[3][3];
    double *N;
    double* gw = el.GaussWeights();
    vec3d gradN[FEElement::MAX_NODES];
    double *Grn, *Gsn, *Gtn;
    double Gr, Gs, Gt;
    
    vec3d b, kpu;
    mat3ds gradb;
    mat3d Kw, Kuu;
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        
        // get the body force
        b = BF.force(mp);
        
        // get the body force stiffness
        gradb = BF.stiffness(mp);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*m_pMat->SolidDensity(mp);
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4dmm dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
        N = el.H(n);
        
        // calculate jacobian
        detJt = invjact(el, Ji, n)*gw[n];
        
        Grn = el.Gr(n);
        Gsn = el.Gs(n);
        Gtn = el.Gt(n);
        
        for (int i=0; i<neln; ++i)
        {
            Gr = Grn[i];
            Gs = Gsn[i];
            Gt = Gtn[i];
            
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN[i] = vec3d(Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt,
                             Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt,
                             Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt);
        }
        
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                Kw = b & gradN[j];
                Kuu = (gradb*(N[j]*rho) + Kw*rhoTw)*(N[i]*detJt);
                ke[4*i  ][4*j  ] += Kuu(0,0); ke[4*i  ][4*j+1] += Kuu(0,1); ke[4*i  ][4*j+2] += Kuu(0,2);
                ke[4*i+1][4*j  ] += Kuu(1,0); ke[4*i+1][4*j+1] += Kuu(1,1); ke[4*i+1][4*j+2] += Kuu(1,2);
                ke[4*i+2][4*j  ] += Kuu(2,0); ke[4*i+2][4*j+1] += Kuu(2,1); ke[4*i+2][4*j+2] += Kuu(2,2);
                
                kpu = (vdotTdotv(gradN[i], dKdE, gradN[j])*b
                       + (Kw + gradb*N[j])*K*gradN[i])*(rhoTw*detJt*dt);
                ke[4*i+3][4*j  ] -= kpu.x; ke[4*i+3][4*j+1] -= kpu.y; ke[4*i+3][4*j+2] -= kpu.z;
            }
    }
}

//-----------------------------------------------------------------------------
vec3d FEBiphasicSolidDomain::FluidFlux(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	
	// pressure gradient
	vec3d gradp = ppt.m_gradp;
	
	// fluid flux w = -k*grad(p)
	mat3ds kt = m_pMat->Permeability(mp);
    
    vec3d w = -(kt*gradp);
    
    double tau = m_pMat->m_tau;
    if (tau > 0) {
        double dt = GetFEModel()->GetTime().timeIncrement;
        w -= kt*(gradp - ppt.m_gradpp)*(tau/dt);
    }
    
	// get true fluid density
	double rhoTw = m_pMat->FluidDensity();
    
    // body force contribution
	FEModel& fem = *m_pMat->GetFEModel();
    int nbf = fem.ModelLoads();
    if (nbf) {
        vec3d b(0,0,0);
        for (int i=0; i<nbf; ++i)
		{
			FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.ModelLoad(i));
			if (pbf && pbf->IsActive())
			{
				// negate b because body forces are defined with a negative sign in FEBio
				b -= pbf->force(mp);
			}
		}
		w += (kt*b)*(rhoTw);
    }
    
    // active momentum supply contribution
	FEActiveMomentumSupply* pAmom = m_pMat->GetActiveMomentumSupply();
    if (pAmom) {
        vec3d pw = pAmom->ActiveSupply(mp);
        w += kt*pw;
    }
    
    return w;
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::UpdateNodalPressures()
{
	vector<double> pi(FEElement::MAX_INTPOINTS);
	vector<double> pn(FEElement::MAX_NODES);

	int NN = Nodes();
	vector<int> tag(NN, 0);
	m_nodePressure.assign(NN, 0.0);

	for (int i = 0; i<Elements(); ++i)
	{
		FESolidElement& el = Element(i);

		int nint = el.GaussPoints();
		int neln = el.Nodes();

		// get integration point pressures
		double pavg = 0.0;
		int c = 0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());

			if (pt) { pavg += pt->m_pa; c++; }
		}
		if (c > 0) pavg /= (double) c;

		// store the nodal values
		for (int j=0; j<neln; ++j)
		{
			int m = el.m_lnode[j];
			m_nodePressure[m] += pavg;
			tag[m]++;
		}
	}

	for (int i=0; i<NN; ++i)
		if (tag[i] > 0) m_nodePressure[i] /= (double) tag[i];
}

//-----------------------------------------------------------------------------
// Note that the data vector stores the values for all of the nodes of the mesh, not just the domain nodes.
// The values will be set to zero for nodes that don't belong to this domain.
void FEBiphasicSolidDomain::GetNodalPressures(vector<double>& data)
{
	FEMesh& mesh = *GetMesh();	
	data.resize(mesh.Nodes(), 0.0);

	int NN = Nodes();
	for (int i=0; i<NN; ++i)
	{
		data[NodeIndex(i)] = m_nodePressure[i];
	}
}
