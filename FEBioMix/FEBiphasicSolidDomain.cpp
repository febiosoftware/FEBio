#include "FEBiphasicSolidDomain.h"
#include "FECore/FEMesh.h"
#include "FECore/log.h"
#include <FECore/FEDataExport.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEBiphasicSolidDomain::FEBiphasicSolidDomain(FEModel* pfem) : FESolidDomain(pfem), FEBiphasicDomain(pfem)
{
	EXPORT_DATA(PLT_FLOAT, FMT_NODE, &m_nodePressure, "NPR fluid pressure");
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::SetMaterial(FEMaterial* pmat)
{
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
            if (el.m_bitfc.size()>0 && el.m_bitfc[i] && node.m_ID[m_dofQ] > -1)
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
			pt.m_r0 = r0;
			pt.m_rt = rt;

			pt.m_J = defgrad(el, pt.m_F, j);

            pb.m_Jp = pt.m_J;
            
            pb.m_p = p;
            pb.m_gradp = gradient(el, pn, j);
            pb.m_gradpp = pb.m_gradp;
            pb.m_phi0p = pb.m_phi0;
            
            mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
bool FEBiphasicSolidDomain::Initialize()
{
	// initialize base class
	FESolidDomain::Initialize();
    
    // error flag (set true on error)
    bool bmerr = false;
    
    // initialize body forces
	FEModel& fem = *GetFEModel();
	m_pMat->m_bf.clear();
    for (int j=0; j<fem.BodyLoads(); ++j)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(j));
        if (pbf) m_pMat->m_bf.push_back(pbf);
    }

	// initialize local coordinate systems (can I do this elsewhere?)
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
	}

    // check for initially inverted elements
    for (int i=0; i<Elements(); ++i)
    {
        FESolidElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int n=0; n<nint; ++n)
        {
            double J0 = detJ0(el, n);
            if (J0 <= 0)
            {
                felog.printf("**************************** E R R O R ****************************\n");
                felog.printf("Negative jacobian detected at integration point %d of element %d\n", n+1, el.GetID());
                felog.printf("Jacobian = %lg\n", J0);
                felog.printf("Did you use the right node numbering?\n");
                felog.printf("Nodes:");
                for (int l=0; l<el.Nodes(); ++l)
                {
                    felog.printf("%d", el.m_node[l]+1);
                    if (l+1 != el.Nodes()) felog.printf(","); else felog.printf("\n");
                }
                felog.printf("*******************************************************************\n\n");
                bmerr = true;
            }
        }
    }
    
	// allocate nodal pressures
	m_nodePressure.resize(Nodes(), 0.0);
    
	return (bmerr == false);
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
				node.m_ID[m_dofX] = DOF_ACTIVE;
				node.m_ID[m_dofY] = DOF_ACTIVE;
				node.m_ID[m_dofZ] = DOF_ACTIVE;
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
                node.m_ID[m_dofQ] = DOF_ACTIVE;
            else
                node.m_ID[m_dofP] = DOF_ACTIVE;
        }
    }
}

//-----------------------------------------------------------------------------
//! Unpack the element LM data. 
void FEBiphasicSolidDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*7);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

        // first the displacement dofs
        lm[4*i  ] = id[m_dofX];
        lm[4*i+1] = id[m_dofY];
        lm[4*i+2] = id[m_dofZ];
        
        // now the pressure dofs
        lm[4*i+3] = id[m_dofP];
        
        // rigid rotational dofs
        // TODO: Do I really need this?
        lm[4*N + 3*i  ] = id[m_dofRU];
        lm[4*N + 3*i+1] = id[m_dofRV];
        lm[4*N + 3*i+2] = id[m_dofRW];
    }
    
    // substitute interface dofs for solid-shell interfaces
    for (int i=0; i<el.m_bitfc.size(); ++i)
    {
        if (el.m_bitfc[i]) {
            FENode& node = m_pMesh->Node(el.m_node[i]);
            vector<int>& id = node.m_ID;
            
            // first the back-face displacement dofs
            lm[4*i  ] = id[m_dofU];
            lm[4*i+1] = id[m_dofV];
            lm[4*i+2] = id[m_dofW];
            
            // now the pressure dof (if the shell has it)
            if (id[m_dofQ] > -1) lm[4*i+3] = id[m_dofQ];
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::Reset()
{
	// reset base class data
	FESolidDomain::Reset();

	// get the biphasic material
	FEBiphasic* pmb = m_pMat;

	// initialize all element data
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];
		
		// get the number of integration points
		int nint = el.GaussPoints();
		
		// loop over the integration points
		for (int n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
			// initialize referential solid volume fraction
			pt.m_phi0 = pmb->m_phi0;
		}
	}
}


//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::InternalForces(FEGlobalVector& R)
{
	int NE = m_Elem.size();
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
		ElementInternalForce(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		//#pragma omp critical
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FEBiphasicSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
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
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        s = pt.m_s;
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // next we get the determinant
        double Jp = bpt.m_Jp;
        double J = pt.m_J;
        
        // and then finally
        double divv = ((J-Jp)/dt)/J;
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solvent supply
        double phiwhat = m_pMat->SolventSupply(mp);
        
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
            fe[4*i  ] -= fu.x*detJt;
            
            fe[4*i+1] -= fu.y*detJt;
            
            fe[4*i+2] -= fu.z*detJt;
            
            fe[4*i+3] -= dt*(w*gradN + (phiwhat - divv)*H[i])*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::InternalForcesSS(FEGlobalVector& R)
{
    int NE = m_Elem.size();
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
        //#pragma omp critical
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements (steady-state)

void FEBiphasicSolidDomain::ElementInternalForceSS(FESolidElement& el, vector<double>& fe)
{
    int i, n;
    
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    
    vec3d gradN, GradN;
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
        
        // calculate the jacobian
        detJt = invjact(el, Ji, n);
        
        detJt *= gw[n];
        
        // get the stress vector for this integration point
        s = pt.m_s;
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        H = el.H(n);
        
        // get the flux
        vec3d& w = bpt.m_w;
        
        // get the solvent supply
        double phiwhat = m_pMat->SolventSupply(mp);
        
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
            fe[4*i  ] -= fu.x*detJt;
            
            fe[4*i+1] -= fu.y*detJt;
            
            fe[4*i+2] -= fu.z*detJt;
            
            fe[4*i+3] -= dt*(w*gradN + phiwhat*H[i])*detJt;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrix(FESolver* psolver, bool bsymm)
{
	// repeat over all solid elements
	int NE = m_Elem.size();
    
    #pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> lm;
		
		FESolidElement& el = m_Elem[iel];
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffness(el, ke, bsymm);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
		UnpackLM(el, lm);

        // assemble element matrix in global stiffness matrix
        #pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::StiffnessMatrixSS(FESolver* psolver, bool bsymm)
{
	// repeat over all solid elements
	int NE = m_Elem.size();

	#pragma omp parallel for shared(NE)
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> lm;
		
		FESolidElement& el = m_Elem[iel];
		
		// allocate stiffness matrix
		int neln = el.Nodes();
		int ndof = neln*4;
		ke.resize(ndof, ndof);
		
		// calculate the element stiffness matrix
		ElementBiphasicStiffnessSS(el, ke, bsymm);
		
		// TODO: the problem here is that the LM array that is returned by the UnpackLM
		// function does not give the equation numbers in the right order. For this reason we
		// have to create a new lm array and place the equation numbers in the right order.
		// What we really ought to do is fix the UnpackLM function so that it returns
		// the LM vector in the right order for poroelastic elements.
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		#pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates element stiffness matrix for element iel
//!
bool FEBiphasicSolidDomain::ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm)
{
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Bp-matrix
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    double tau = m_pMat->m_tau;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n);
        
        // contravariant basis vectors in spatial frame
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp + (pt.m_gradp - pt.m_gradpp)*(tau/dt);
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], c, gradN[j]))*tmp;
                
                ke[4*i  ][4*j  ] += Kuu[0][0]; ke[4*i  ][4*j+1] += Kuu[0][1]; ke[4*i  ][4*j+2] += Kuu[0][2];
                ke[4*i+1][4*j  ] += Kuu[1][0]; ke[4*i+1][4*j+1] += Kuu[1][1]; ke[4*i+1][4*j+2] += Kuu[1][2];
                ke[4*i+2][4*j  ] += Kuu[2][0]; ke[4*i+2][4*j+1] += Kuu[2][1]; ke[4*i+2][4*j+2] += Kuu[2][2];
            }
        
        // calculate the kpp matrix
        tmp = detJ*gw[n]*dt;
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                ke[4*i+3][4*j+3] += (H[i]*H[j]*Phip - gradN[i]*(K*gradN[j])*(1+tau/dt))*tmp;
            }
        
        if (!bsymm) {
            // calculate the kup matrix
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    tmp = detJ*gw[n]*H[j];
                    ke[4*i  ][4*j+3] -= tmp*gradN[i].x;
                    ke[4*i+1][4*j+3] -= tmp*gradN[i].y;
                    ke[4*i+2][4*j+3] -= tmp*gradN[i].z;
                }
            }
            
            // calculate the kpu matrix
            //			tmp = detJ*gw[n];
            tmp = detJ*gw[n]*dt;
            mat3ds Q = Phie*ept.m_J + mat3dd(phiwhat - 1./dt);
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    vec3d vt = ((vdotTdotv(-gradN[i], dKdE, gradN[j]).transpose()*gradp)
                                +(Q*gradN[j])*H[i])*tmp;
                    ke[4*i+3][4*j  ] += vt.x;
                    ke[4*i+3][4*j+1] += vt.y;
                    ke[4*i+3][4*j+2] += vt.z;
                }
            }
            
        } else {
            // calculate the kup matrix and let kpu be its symmetric part
            tmp = detJ*gw[n];
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    ke[4*i  ][4*j+3] -= tmp*H[j]*gradN[i].x;
                    ke[4*i+1][4*j+3] -= tmp*H[j]*gradN[i].y;
                    ke[4*i+2][4*j+3] -= tmp*H[j]*gradN[i].z;
                    
                    ke[4*i+3][4*j  ] -= tmp*H[i]*gradN[j].x;
                    ke[4*i+3][4*j+1] -= tmp*H[i]*gradN[j].y;
                    ke[4*i+3][4*j+2] -= tmp*H[i]*gradN[j].z;
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
    int i, j, n;
    
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    double *Gr, *Gs, *Gt, *H;
    
    // jacobian
    double Ji[3][3], detJ;
    
    // Bp-matrix
    vector<vec3d> gradN(neln);
    double tmp;
    
    // gauss-weights
    double* gw = el.GaussWeights();
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    
    // zero stiffness matrix
    ke.zero();
    
    // loop over gauss-points
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint >());
        FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
        
        // calculate jacobian
        detJ = invjact(el, Ji, n);
        
        // contravariant basis vectors in spatial frame
        vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
        vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
        vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
        
        H = el.H(n);
        
        Gr = el.Gr(n);
        Gs = el.Gs(n);
        Gt = el.Gt(n);
        
        for (i=0; i<neln; ++i)
        {
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
        }
        
        // get stress tensor
        mat3ds s = ept.m_s;
        
        // get elasticity tensor
        tens4ds c = m_pMat->Tangent(mp);
        
        // get the fluid flux and pressure gradient
        vec3d gradp = pt.m_gradp;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                mat3d Kuu = (mat3dd(gradN[i]*(s*gradN[j])) + vdotTdotv(gradN[i], c, gradN[j]))*tmp;
                
                ke[4*i  ][4*j  ] += Kuu[0][0]; ke[4*i  ][4*j+1] += Kuu[0][1]; ke[4*i  ][4*j+2] += Kuu[0][2];
                ke[4*i+1][4*j  ] += Kuu[1][0]; ke[4*i+1][4*j+1] += Kuu[1][1]; ke[4*i+1][4*j+2] += Kuu[1][2];
                ke[4*i+2][4*j  ] += Kuu[2][0]; ke[4*i+2][4*j+1] += Kuu[2][1]; ke[4*i+2][4*j+2] += Kuu[2][2];
            }
        
        // calculate the kpp matrix
        tmp = detJ*gw[n]*dt;
        for (i=0; i<neln; ++i)
            for (j=0; j<neln; ++j)
            {
                ke[4*i+3][4*j+3] += (H[i]*H[j]*Phip - gradN[i]*(K*gradN[j]))*tmp;
            }
        
        if (!bsymm) {
            // calculate the kup matrix
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    tmp = detJ*gw[n]*H[j];
                    ke[4*i  ][4*j+3] -= tmp*gradN[i].x;
                    ke[4*i+1][4*j+3] -= tmp*gradN[i].y;
                    ke[4*i+2][4*j+3] -= tmp*gradN[i].z;
                }
            }
            
            // calculate the kpu matrix
            //			tmp = detJ*gw[n];
            tmp = detJ*gw[n]*dt;
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    vec3d vt = ((vdotTdotv(-gradp, dKdE, gradN[j]).transpose()*(gradN[i]))
                                +(mat3dd(phiwhat) + Phie*ept.m_J)*gradN[j]*H[i])*tmp;
                    ke[4*i+3][4*j  ] += vt.x;
                    ke[4*i+3][4*j+1] += vt.y;
                    ke[4*i+3][4*j+2] += vt.z;
                }
            }
            
        } else {
            // calculate the kup matrix and let kpu be its symmetric part
            tmp = detJ*gw[n];
            for (i=0; i<neln; ++i) {
                for (j=0; j<neln; ++j)
                {
                    ke[4*i  ][4*j+3] -= tmp*H[j]*gradN[i].x;
                    ke[4*i+1][4*j+3] -= tmp*H[j]*gradN[i].y;
                    ke[4*i+2][4*j+3] -= tmp*H[j]*gradN[i].z;
                    
                    ke[4*i+3][4*j  ] -= tmp*H[i]*gradN[j].x;
                    ke[4*i+3][4*j+1] -= tmp*H[i]*gradN[j].y;
                    ke[4*i+3][4*j+2] -= tmp*H[i]*gradN[j].z;
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
				if (NegativeJacobian::m_boutput) e.print();
			}
		}
	}
	// if we encountered an error, we request a running restart
	if (berr)
	{
		if (NegativeJacobian::m_boutput == false) felog.printbox("ERROR", "Negative jacobian was detected.");
		throw DoRunningRestart();
	}

	// also update the nodal pressures
	UpdateNodalPressures();
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::UpdateElementStress(int iel)
{
   // extract the elastic component
    FEElasticMaterial* pme = m_pMat->GetElasticMaterial();

	// get the solid element
	FESolidElement& el = m_Elem[iel];
		
	// get the number of integration points
	int nint = el.GaussPoints();
		
	// get the number of nodes
	int neln = el.Nodes();

	// get the nodal data
	FEMesh& mesh = *m_pMesh;
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	double pn[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
        FENode& node = mesh.Node(el.m_node[j]);
		r0[j] = node.m_r0;
		rt[j] = node.m_rt;
        if (el.m_bitfc.size()>0 && el.m_bitfc[j] && node.m_ID[m_dofQ] > -1)
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
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);
			
		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);
			
		// poroelasticity data
		FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			
		// evaluate fluid pressure at gauss-point
		ppt.m_p = el.Evaluate(pn, n);
			
		// calculate the gradient of p at gauss-point
		ppt.m_gradp = gradient(el, pn, n);
			
		// for biphasic materials also update the fluid flux
		ppt.m_w = FluidFlux(mp);
		ppt.m_pa = m_pMat->Pressure(mp);
			
		// calculate the stress at this material point
		pt.m_s = m_pMat->Stress(mp);
	}
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NE = (int)m_Elem.size();
#pragma omp parallel for
    for (int i=0; i<NE; ++i)
    {
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FESolidElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 4*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces
        ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the body forces

void FEBiphasicSolidDomain::ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe)
{
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
    double rhoTw = m_pMat->FluidDensity();
    
    // jacobian
    double detJ;
    double *H;
    double* gw = el.GaussWeights();
    vec3d b;
    
    // number of nodes
    int neln = el.Nodes();
    
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;
        rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;
    }
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
        pt.m_r0 = el.Evaluate(r0, n);
        pt.m_rt = el.Evaluate(rt, n);
        
        detJ = detJt(el, n)*gw[n];
        
        // get the force
        b = BF.force(mp);
        
        // evaluate apparent solid and fluid densities and mixture density
        double phiw = m_pMat->Porosity(mp);
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        H = el.H(n);
        
        for (int i=0; i<neln; ++i)
        {
            fe[4*i  ] -= H[i]*rho*b.x*detJ;
            fe[4*i+1] -= H[i]*rho*b.y*detJ;
            fe[4*i+2] -= H[i]*rho*b.z*detJ;
        }
    }
}

//-----------------------------------------------------------------------------
void FEBiphasicSolidDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
    FEBiphasic* pmb = dynamic_cast<FEBiphasic*>(GetMaterial()); assert(pmb);
    
    // element stiffness matrix
    matrix ke;
    vector<int> lm;
    
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
    for (int iel=0; iel<NE; ++iel)
    {
        FESolidElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
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
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEBiphasicSolidDomain::ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke)
{
    int neln = el.Nodes();
    
    // get true solid and fluid densities
    double rhoTs = m_pMat->SolidDensity();
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
        double rhos = (1-phiw)*rhoTs;
        double rhow = phiw*rhoTw;
        double rho = rhos + rhow;
        
        // evaluate the permeability and its derivatives
        mat3ds K = m_pMat->Permeability(mp);
        tens4ds dKdE = m_pMat->GetPermeability()->Tangent_Permeability_Strain(mp);
        
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
                
                kpu = (vdotTdotv(gradN[i], dKdE, gradN[j]).transpose()*b
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
    int nbf = fem.BodyLoads();
    if (nbf) {
        vec3d b(0,0,0);
        for (int i=0; i<nbf; ++i)
		{
			FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
			if (pbf->IsActive())
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
