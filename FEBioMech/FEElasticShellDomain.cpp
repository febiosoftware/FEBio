#include "stdafx.h"
#include "FEElasticShellDomain.h"
#include "FEElasticMaterial.h"
#include "FEBodyForce.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <math.h>

//-----------------------------------------------------------------------------
FEElasticShellDomain::FEElasticShellDomain(FEModel* pfem) : FEShellDomain(&pfem->GetMesh()), FEElasticDomain(pfem)
{
	m_pMat = 0;
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
bool FEElasticShellDomain::Initialize()
{
	// initialize base class
	FEShellDomain::Initialize();

	// error flag (set true on error)
	bool bmerr = false;

    // get the elements material
    if (m_pMat)
    {
        FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
        if (pme)
        {
            // assign local coordinate system to each integration point
            for (size_t i=0; i<m_Elem.size(); ++i)
            {
                FEShellElement& el = m_Elem[i];
                for (int n=0; n<el.GaussPoints(); ++n) pme->SetLocalCoordinateSystem(el, n, *(el.GetMaterialPoint(n)));
            }
        }
    }
    
	// check for initially inverted shells
	for (int i=0; i<Elements(); ++i)
	{
		FEShellElement& el = Element(i);
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

	return (bmerr == false);
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.m_bexclude == false)
		{
			if (node.m_rid < 0)
			{
				node.m_ID[m_dofX] = DOF_ACTIVE;
				node.m_ID[m_dofY] = DOF_ACTIVE;
				node.m_ID[m_dofZ] = DOF_ACTIVE;
			}

			if (node.m_bshell)
			{
				node.m_ID[m_dofU] = DOF_ACTIVE;
				node.m_ID[m_dofV] = DOF_ACTIVE;
				node.m_ID[m_dofW] = DOF_ACTIVE;
                node.m_ID[m_dofAU] = DOF_ACTIVE;
                node.m_ID[m_dofAV] = DOF_ACTIVE;
                node.m_ID[m_dofAW] = DOF_ACTIVE;
            }
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FEElasticShellDomain::CoBaseVectors(FEShellElement& el, int n, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rt;
        D[i] = ni.m_d0 + ni.get_vec3d(m_dofX, m_dofY, m_dofZ) - ni.get_vec3d(m_dofU, m_dofV, m_dofW);
    }
    
    double eta = el.gt(n);
    
    double* Mr = el.Hr(n);
    double* Ms = el.Hs(n);
    double* M  = el.H(n);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (i=0; i<neln; ++i)
    {
        g[0] += (r[i] - D[i]*(1-eta)/2)*Mr[i];
        g[1] += (r[i] - D[i]*(1-eta)/2)*Ms[i];
        g[2] += D[i]*(M[i]/2);
    }
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FEElasticShellDomain::ContraBaseVectors(FEShellElement& el, int n, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));

}

//-----------------------------------------------------------------------------
// jacobian with respect to current frame
double FEElasticShellDomain::detJ(FEShellElement& el, int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
double FEElasticShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
    CoBaseVectors(el, n, gcov);
    ContraBaseVectors0(el, n, Gcnt);
    
    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
	double J = F.det();
	if (J <= 0) throw NegativeJacobian(el.GetID(), n, J, &el);

	return J;
}

//-----------------------------------------------------------------------------
//! evaluate a vector function over the shell
vec3d FEElasticShellDomain::evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n)
{
    const double *M;
    double eta;
    vec3d gcnt[3];
    double Mu, Md;
    vec3d v(0,0,0);
    
    eta = el.gt(n);
    
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        Mu = (1+eta)/2*M[i];
        Md = (1-eta)/2*M[i];
        v += vn[i]*Mu + dvn[i]*Md;
    }
    
    return v;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the current frame at 
//! integration point n. The inverse jacobian is return in Ji. The return value
//! is the determinant of the jacobian (not the inverse!)
double FEElasticShellDomain::invjact(FEShellElement& el, double Ji[3][3], int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate the inverse of the jacobian
    double deti = 1.0 / det;
    
    Ji[0][0] =  deti*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Ji[1][0] =  deti*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    Ji[2][0] =  deti*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    Ji[0][1] =  deti*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    Ji[1][1] =  deti*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Ji[2][1] =  deti*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    
    Ji[0][2] =  deti*(J[0][1]*J[1][2] - J[1][1]*J[0][2]);
    Ji[1][2] =  deti*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    Ji[2][2] =  deti*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
// Calculates the forces due to the stress
void FEElasticShellDomain::InternalForces(FEGlobalVector& R)
{
    int NS = m_Elem.size();
#pragma omp parallel for shared (NS)
    for (int i=0; i<NS; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // create the element force vector and initialize to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate element's internal force
        ElementInternalForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble the residual
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for shell elements
//! Note that we use a one-point gauss integration rule for the thickness
//! integration. This will integrate linear functions exactly.

void FEElasticShellDomain::ElementInternalForce(FEShellElement& el, vector<double>& fe)
{
	int i, n;

	// jacobian matrix determinant
	double detJt;

	const double* Mr, *Ms, *M;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();
	double eta;
    
    vec3d gcnt[3];

	// repeat for all integration points
	for (n=0; n<nint; ++n)
	{
		FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
        detJt = detJ(el, n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.m_s;

		eta = el.gt(n);

		Mr = el.Hr(n);
		Ms = el.Hs(n);
		M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt);

		for (i=0; i<neln; ++i)
		{
            vec3d gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
            vec3d gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
            vec3d gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
            vec3d fu = s*gradMu;
            vec3d fd = s*gradMd;
            
            // calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[6*i  ] -= fu.x*detJt;
			fe[6*i+1] -= fu.y*detJt;
			fe[6*i+2] -= fu.z*detJt;

			fe[6*i+3] -= fd.x*detJt;
			fe[6*i+4] -= fd.y*detJt;
			fe[6*i+5] -= fd.z*detJt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& BF)
{
    int NS = m_Elem.size();
#pragma omp parallel for
    for (int i=0; i<NS; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // create the element force vector and initialize to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // apply body forces to shells
        ElementBodyForce(BF, el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble the residual
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! Calculates element body forces for shells

void FEElasticShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
    double dens = m_pMat->Density();
    
    // integration weights
    double* gw = el.GaussWeights();
    double eta;
    double *M, detJt;
    
    // loop over integration points
    int nint = el.GaussPoints();
    int neln = el.Nodes();
    
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        
        // calculate the jacobian
        detJt = detJ0(el, n)*gw[n];
        
        M  = el.H(n);
        eta = el.gt(n);
        
        // get the force
        vec3d f = BF.force(mp);
        
        for (int i=0; i<neln; ++i)
        {
            vec3d fu = f*(dens*M[i]*(1+eta)/2*detJt);
            vec3d fd = f*(dens*M[i]*(1-eta)/2*detJt);
            
            fe[6*i  ] -= fu.x;
            fe[6*i+1] -= fu.y;
            fe[6*i+2] -= fu.z;
            
            fe[6*i+3] -= fd.x;
            fe[6*i+4] -= fd.y;
            fe[6*i+5] -= fd.z;
        }
    }
}

//-----------------------------------------------------------------------------
// Calculate inertial forces \todo Why is F no longer needed?
void FEElasticShellDomain::InertialForces(FEGlobalVector& R, vector<double>& F)
{
    FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(GetMaterial()); assert(pme);
    double d = pme->Density();
    
    const int MN = FEElement::MAX_NODES;
    vec3d at[MN], aqt[MN];
    
    // acceleration at integration point
    vec3d a;
    
    // loop over all elements
    int NE = Elements();
    
    for (int iel=0; iel<NE; ++iel)
    {
        vector<double> fe;
        vector<int> lm;
        
        FEShellElement& el = Element(iel);
        
        int nint = el.GaussPoints();
        int neln = el.Nodes();
        
        fe.assign(6*neln,0);
        
        // get the nodal accelerations
        for (int i=0; i<neln; ++i)
        {
            at[i] = m_pMesh->Node(el.m_node[i]).m_at;
            aqt[i] = m_pMesh->Node(el.m_node[i]).get_vec3d(m_dofAU, m_dofAV, m_dofAW);
        }
        
        // evaluate the element inertial force vector
        for (int n=0; n<nint; ++n)
        {
            double J0 = detJ0(el, n)*el.GaussWeights()[n];
            
            // get the acceleration for this integration point
            a = evaluate(el, at, aqt, n);
            
            double* M = el.H(n);
            double eta = el.gt(n);
            
            for (int i=0; i<neln; ++i)
            {
                vec3d fu = a*(d*M[i]*(1+eta)/2*J0);
                vec3d fd = a*(d*M[i]*(1-eta)/2*J0);
                
                fe[6*i  ] -= fu.x;
                fe[6*i+1] -= fu.y;
                fe[6*i+2] -= fu.z;
                
                fe[6*i+3] -= fd.x;
                fe[6*i+4] -= fd.y;
                fe[6*i+5] -= fd.z;
            }
        }
        
        // get the element degrees of freedom
        UnpackLM(el, lm);
        
        // assemble fe into R
        R.Assemble(el.m_node, lm, fe);
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness due to body forces
void FEElasticShellDomain::ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke)
{
    int i, j, i6, j6;
    int neln = el.Nodes();
    
    // don't forget to multiply with the density
    double dens = m_pMat->Density();
    
    // jacobian
    double detJ;
    double *M;
    double* gw = el.GaussWeights();
    mat3ds K;
    
    double Mu[FEElement::MAX_NODES], Md[FEElement::MAX_NODES];
    
    // loop over integration points
    int nint = el.GaussPoints();
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        detJ = detJ0(el, n)*gw[n];
        
        // get the stiffness
        K = BF.stiffness(mp)*dens*detJ;
        
        M = el.H(n);
        
        double eta = el.gt(n);
        
        for (i=0; i<neln; ++i)
        {
            Mu[i] = M[i]*(1+eta)/2;
            Md[i] = M[i]*(1-eta)/2;
        }
        
        for (i=0, i6=0; i<neln; ++i, i6 += 6)
        {
            for (j=0, j6 = 0; j<neln; ++j, j6 += 6)
            {
                mat3d Kuu = K*(Mu[i]*Mu[j]);
                mat3d Kud = K*(Mu[i]*Md[j]);
                mat3d Kdu = K*(Md[i]*Mu[j]);
                mat3d Kdd = K*(Md[i]*Md[j]);
                
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
}

//-----------------------------------------------------------------------------

void FEElasticShellDomain::StiffnessMatrix(FESolver* psolver)
{
    // repeat over all shell elements
    int NS = m_Elem.size();
#pragma omp parallel for shared (NS)
    for (int iel=0; iel<NS; ++iel)
    {
        matrix ke;
        vector<int> lm;
        
        FEShellElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 6*el.Nodes();
        ke.resize(ndof, ndof);
        
        // calculate the element stiffness matrix
        ElementStiffness(iel, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
        
    }
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::MassMatrix(FESolver* psolver, double scale)
{
    // repeat over all solid elements
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEShellElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 6*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementMassMatrix(el, ke, scale);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::BodyForceStiffness(FESolver* psolver, FEBodyForce& bf)
{
    // repeat over all shell elements
    int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
    for (int iel=0; iel<NE; ++iel)
    {
        // element stiffness matrix
        matrix ke;
        vector<int> lm;
        
        FEShellElement& el = m_Elem[iel];
        
        // create the element's stiffness matrix
        int ndof = 6*el.Nodes();
        ke.resize(ndof, ndof);
        ke.zero();
        
        // calculate inertial stiffness
        ElementBodyForceStiffness(bf, el, ke);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element matrix in global stiffness matrix
#pragma omp critical
        psolver->AssembleStiffness(el.m_node, lm, ke);
    }
}

//-----------------------------------------------------------------------------
//! Calculates the shell element stiffness matrix

void FEElasticShellDomain::ElementStiffness(int iel, matrix& ke)
{
    FEShellElement& el = Element(iel);
    
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
    
    // calculate element stiffness matrix
    ke.zero();
    for (n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        
        // calculate the jacobian
        detJt = detJ(el, n);
        
        detJt *= gw[n];
        
        // get the stress and elasticity for this integration point
        mat3ds s = pt.m_s;
        tens4ds C = m_pMat->Tangent(mp);
        
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
//! calculates element inertial stiffness matrix
void FEElasticShellDomain::ElementMassMatrix(FEShellElement& el, matrix& ke, double a)
{
    // Get the current element's data
    const int nint = el.GaussPoints();
    const int neln = el.Nodes();
    
    // weights at gauss points
    const double *gw = el.GaussWeights();
    
    // density
    double D = m_pMat->Density();
    
    // calculate element stiffness matrix
    for (int n=0; n<nint; ++n)
    {
        // shape functions
        double* M = el.H(n);
        
        // Jacobian
        double J0 = detJ0(el, n)*gw[n];
        
        // parametric coordinate through thickness
        double eta = el.gt(n);
        
        for (int i=0; i<neln; ++i)
            for (int j=0; j<neln; ++j)
            {
                double Kuu = (1+eta)/2*M[i]*(1+eta)/2*M[j]*a*D*J0;
                double Kud = (1+eta)/2*M[i]*(1-eta)/2*M[j]*a*D*J0;
                double Kdu = (1-eta)/2*M[i]*(1+eta)/2*M[j]*a*D*J0;
                double Kdd = (1-eta)/2*M[i]*(1-eta)/2*M[j]*a*D*J0;
                
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
    }
    
}

//-----------------------------------------------------------------------------
//! Calculates body forces for shells

void FEElasticShellDomain::ElementBodyForce(FEModel& fem, FEShellElement& el, vector<double>& fe)
{
    int NF = fem.BodyLoads();
    for (int nf = 0; nf < NF; ++nf)
    {
        FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(nf));
        if (pbf)
        {
            double dens0 = m_pMat->Density();
            
            // integration weights
            double* gw = el.GaussWeights();
            double eta;
            double *M, detJt;
            
            // loop over integration points
            int nint = el.GaussPoints();
            int neln = el.Nodes();
            
            for (int n=0; n<nint; ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
                
                // calculate density in current configuration
                double dens = dens0/pt.m_J;
                
                // calculate the jacobian
                detJt = detJ(el, n)*gw[n];
                
                M  = el.H(n);
                eta = el.gt(n);
                
                // get the force
                vec3d f = pbf->force(mp);
                
                for (int i=0; i<neln; ++i)
                {
                    vec3d fu = f*(dens*M[i]*(1+eta)/2);
                    vec3d fd = f*(dens*M[i]*(1-eta)/2);
                    
                    fe[6*i  ] -= fu.x*detJt;
                    fe[6*i+1] -= fu.y*detJt;
                    fe[6*i+2] -= fu.z*detJt;
                    
                    fe[6*i+3] -= fd.x*detJt;
                    fe[6*i+4] -= fd.y*detJt;
                    fe[6*i+5] -= fd.z*detJt;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::Update(const FETimeInfo& tp)
{
	FEMesh& mesh = *GetMesh();
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];

	int n;
	for (int i=0; i<(int) m_Elem.size(); ++i)
	{
		// get the solid element
		FEShellElement& el = m_Elem[i];

		// get the number of integration points
		int nint = el.GaussPoints();

		// number of nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (int j=0; j<neln; ++j)
		{
			r0[j] = mesh.Node(el.m_node[j]).m_r0;
			rt[j] = mesh.Node(el.m_node[j]).m_rt;
		}

		// loop over the integration points and calculate
		// the stress at the integration point
		for (n=0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(n));
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

			// material point coordinates
			// TODO: I'm not entirly happy with this solution
			//		 since the material point coordinates are used by most materials.
			pt.m_r0 = el.Evaluate(r0, n);
			pt.m_rt = el.Evaluate(rt, n);

			// get the deformation gradient and determinant
			pt.m_J = defgrad(el, pt.m_F, n);

			// calculate the stress at this material point
			pt.m_s = m_pMat->Stress(mp);
		}
	}
}


//-----------------------------------------------------------------------------
//! Unpack the element. That is, copy element data in traits structure
//! Note that for the shell elements the lm order is different compared
//! to the solid element ordering. This is because for shell elements the
//! nodes have six degrees of freedom each, where for solids they only
//! have 3 dofs.
void FEElasticShellDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*9);
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[6*i  ] = id[m_dofX];
		lm[6*i+1] = id[m_dofY];
		lm[6*i+2] = id[m_dofZ];

		// next the rotational dofs
		lm[6*i+3] = id[m_dofU];
		lm[6*i+4] = id[m_dofV];
		lm[6*i+5] = id[m_dofW];

		// rigid rotational dofs
		lm[6*N + 3*i  ] = id[m_dofRU];
		lm[6*N + 3*i+1] = id[m_dofRV];
		lm[6*N + 3*i+2] = id[m_dofRW];
	}
}
