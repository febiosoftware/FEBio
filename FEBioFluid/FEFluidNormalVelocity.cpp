/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEFluidNormalVelocity.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/log.h"
#include "FECore/LinearSolver.h"
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidNormalVelocity, FESurfaceLoad)
	ADD_PARAMETER(m_velocity, "velocity");
	ADD_PARAMETER(m_VC      , "value"   );
	ADD_PARAMETER(m_bpv     , "prescribe_nodal_velocities");
	ADD_PARAMETER(m_bpar    , "parabolic"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalVelocity::FEFluidNormalVelocity(FEModel* pfem) : FESurfaceLoad(pfem), m_VC(FE_DOUBLE), m_dofWE(pfem)
{
    m_velocity = 0.0;
    m_bpv = true;
    m_bpar = false;
    
	m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
	m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNormalVelocity::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_VC.Create(ps, 1.0);
}

//-----------------------------------------------------------------------------
double FEFluidNormalVelocity::NormalVelocity(FESurfaceMaterialPoint& mp)
{
	FESurfaceElement& el = *mp.SurfaceElement();
	double vn = 0;
	double* N = mp.m_shape;
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
	{
		vn += N[i] * m_VN[el.m_lnode[i]];
	}
	return vn;
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal velocity
void FEFluidNormalVelocity::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEDofList dofE(GetFEModel());
	dofE.AddDof(m_dofWE[3]);
	m_psurf->LoadVector(R, dofE, false, [=](FESurfaceMaterialPoint& mp, int node_a, vector<double>& fa) {

		FESurfaceElement& el = *mp.SurfaceElement();

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);

		// calculate the tangent vectors
		vec3d dxr = el.eval_deriv1(rt, mp.m_index);
		vec3d dxs = el.eval_deriv2(rt, mp.m_index);

		// normal velocity
		double vn = NormalVelocity(mp);
		vn *= m_velocity;
		double da = (dxr ^ dxs).norm();

		double* N = mp.m_shape;
		fa[0] = N[node_a] * vn * da;
	});
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidNormalVelocity::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    FEMesh* mesh = ps->GetMesh();

    // evaluate surface normals
    vector<vec3d> sn(ps->Elements(),vec3d(0,0,0));
    m_nu.resize(ps->Nodes(),vec3d(0,0,0));
    m_VN.resize(ps->Nodes(),0);
    vector<int> nf(ps->Nodes(),0);
    vec3d r0[FEElement::MAX_NODES];
    for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            r0[i] = mesh->Node(el.m_node[i]).m_r0;
            m_VN[el.m_lnode[i]] += m_VC.value<double>(iel, i);
            ++nf[el.m_lnode[i]];
        }
        
        double* Nr, *Ns;
        
        vec3d dxr, dxs;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the tangent vectors at integration point
            dxr = dxs = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                dxr += r0[i]*Nr[i];
                dxs += r0[i]*Ns[i];
            }
            
            sn[iel] = dxr ^ dxs;
            sn[iel].unit();
        }
        
        // evaluate nodal normals by averaging surface normals
        for (int i=0; i<neln; ++i)
            m_nu[el.m_lnode[i]] += sn[iel];
    }
    
    for (int i=0; i<ps->Nodes(); ++i) {
        m_nu[i].unit();
        m_VN[i] /= nf[i];
    }
    
    // Set parabolic velocity profile if requested.
    // This will override velocity boundary cards in m_VC
    // and nodal cards in m_VN
    if (m_bpar) return SetParabolicVelocity();
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidNormalVelocity::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofWE[0], DOF_PRESCRIBED);
        node.set_bc(m_dofWE[1], DOF_PRESCRIBED);
        node.set_bc(m_dofWE[2], DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the velocities
void FEFluidNormalVelocity::Update()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        // evaluate the velocity
        vec3d v = m_nu[i]*(m_velocity*m_VN[i]);
        FENode& node = ps->Node(i);
        if (node.m_ID[m_dofWE[0]] < -1) node.set(m_dofWE[0], v.x);
        if (node.m_ID[m_dofWE[1]] < -1) node.set(m_dofWE[1], v.y);
        if (node.m_ID[m_dofWE[2]] < -1) node.set(m_dofWE[2], v.z);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate normal velocities by solving Poiseuille flow across the surface
bool FEFluidNormalVelocity::SetParabolicVelocity()
{
    // find surface boundary nodes
    FEElemElemList EEL;
    FESurface* ps = &GetSurface();
    EEL.Create(ps);

    vector<bool> boundary(ps->Nodes(), false);
    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        for (int j=0; j<el.facet_edges(); ++j) {
            FEElement* nel = EEL.Neighbor(i, j);
            if (nel == nullptr) {
                int en[3] = {-1,-1,-1};
                el.facet_edge(j, en);
                boundary[en[0]] = true;
                boundary[en[1]] = true;
                if (en[2] > -1) boundary[en[2]] = true;
            }
        }
    }
    
    // count number of non-boundary nodes
    int neq = 0;
    vector<int> glm(ps->Nodes(),-1);
    for (int i=0; i<ps->Nodes(); ++i)
        if (!boundary[i]) glm[i] = neq++;
    if (neq == 0)
    {
        feLogError("Unable to set parabolic fluid normal velocity\n");
        return false;
    }
    
    // create a linear solver
    LinearSolver*		plinsolve;	//!< the linear solver
    FEGlobalMatrix*		pK;			//!< stiffness matrix
    FECoreKernel& fecore = FECoreKernel::GetInstance();
    plinsolve = fecore.CreateLinearSolver(0, "skyline");
    if (plinsolve == 0)
    {
		feLogError("Unknown solver type selected\n");
        return false;
    }
    
    SparseMatrix* pS = plinsolve->CreateSparseMatrix(REAL_SYMMETRIC);
    pK = new FEGlobalMatrix(pS);
    if (pK == 0)
    {
		feLogError("Failed allocating stiffness matrix\n\n");
        return false;
    }
    // build matrix profile for normal velocity at non-boundary nodes
    pK->build_begin(neq);
    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        vector<int> elm(el.Nodes(),-1);
        for (int j=0; j<el.Nodes(); ++j)
            elm[j] = glm[el.m_lnode[j]];
        pK->build_add(elm);
    }
    pK->build_end();
    pS->Zero();
    
    // create global vector
    vector<double> v;           //!< normal velocity solution
    vector<double> rhs;         //!< right-hand-side
    vector<double> Fr;          //!< reaction forces
    v.assign(neq, 0);
    rhs.assign(neq, 0);
    Fr.assign(neq, 0);
    FEModel pfem;
    FEGlobalVector pR(pfem, rhs, Fr);
    
    // calculate the global matrix and vector
    FEElementMatrix ke;
    vector<double> fe;
    vector<int> lm;
    
    for (int m=0; m<ps->Elements(); ++m)
    {
        // get the surface element
        FESurfaceElement& el = ps->Element(m);
        
        int neln = el.Nodes();

        // get the element stiffness matrix
        ke.resize(neln, neln);
        lm.resize(neln);
        fe.resize(neln);
        vector<vec3d> gradN(neln);
        
        // calculate stiffness
        int nint = el.GaussPoints();
        
        // gauss weights
        double* w = el.GaussWeights();
        
        // nodal coordinates
        FEMesh& mesh = *ps->GetMesh();
        vec3d rt[FEElement::MAX_NODES];
        for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
        
        // repeat over integration points
        ke.zero();
        zero(fe);
        vec3d gcnt[2];
        for (int n=0; n<nint; ++n)
        {
            double* N = el.H(n);
            double* Gr = el.Gr(n);
            double* Gs = el.Gs(n);
            ps->ContraBaseVectors(el, n, gcnt);
            
            vec3d dxr(0,0,0), dxs(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                dxr += rt[i]*Gr[i];
                dxs += rt[i]*Gs[i];
                gradN[i] = gcnt[0]*Gr[i] + gcnt[1]*Gs[i];
            }
            
            double da = (dxr ^ dxs).norm();
            
            // calculate stiffness component
            for (int i=0; i<neln; ++i) {
                fe[i] += N[i]*w[n]*da;
                for (int j=0; j<neln; ++j)
                    ke[i][j] += (gradN[i]*gradN[j])*w[n]*da;
            }
        }
        
        // get the element's LM vector
        for (int j=0; j<el.Nodes(); ++j)
            lm[j] = glm[el.m_lnode[j]];
        
        // assemble element matrix in global stiffness matrix
		ke.SetIndices(lm);
        pK->Assemble(ke);
        pR.Assemble(lm, fe);
    }

    // solve linear system
    plinsolve->PreProcess();
    plinsolve->Factor();
    if (plinsolve->BackSolve(v, rhs) == false)
    {
		feLogError("Unable to solve for parabolic fluid normal velocity\n");
        return false;
    }
    plinsolve->Destroy();
    
    // set the nodal normal velocity scale factors
    for (int i=0; i<ps->Nodes(); ++i) {
        if (glm[i] == -1) m_VN[i] = 0;
        else m_VN[i] = v[glm[i]];
    }
    
    // evaluate net area and volumetric flow rate
    double A = 0, Q = 0;
    for (int m=0; m<ps->Elements(); ++m)
    {
        // get the surface element
        FESurfaceElement& el = ps->Element(m);
        
        int neln = el.Nodes();
        int nint = el.GaussPoints();
        double* w = el.GaussWeights();
        
        // nodal coordinates
        FEMesh& mesh = *ps->GetMesh();
        vec3d rt[FEElement::MAX_NODES];
        for (int j=0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            double* N = el.H(n);
            double* Gr = el.Gr(n);
            double* Gs = el.Gs(n);

            double vn = 0;
            vec3d dxr(0,0,0), dxs(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                vn += N[i]*m_VN[el.m_lnode[i]];
                dxr += rt[i]*Gr[i];
                dxs += rt[i]*Gs[i];
            }
            
            double da = (dxr ^ dxs).norm();
            
            for (int i=0; i<neln; ++i) {
                A += N[i]*w[n]*da;
                Q += N[i]*vn*w[n]*da;
            }
        }
    }
    
    // normalize nodal velocity cards
    double vbar = Q/A;
    for (int i=0; i<ps->Nodes(); ++i) m_VN[i] /= vbar;
    
    return true;
}

//! serializatsion
void FEFluidNormalVelocity::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
	ar & m_VN;
	ar & m_nu;
}
