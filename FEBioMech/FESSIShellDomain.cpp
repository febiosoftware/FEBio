#include "stdafx.h"
#include "FESSIShellDomain.h"
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
FESSIShellDomain::FESSIShellDomain(FEModel* pfem) : FEShellDomain(&pfem->GetMesh())
{
    m_dofx = pfem->GetDOFIndex("x");
    m_dofy = pfem->GetDOFIndex("y");
    m_dofz = pfem->GetDOFIndex("z");
    m_dofu = pfem->GetDOFIndex("u");
    m_dofv = pfem->GetDOFIndex("v");
    m_dofw = pfem->GetDOFIndex("w");
}

//-----------------------------------------------------------------------------
bool FESSIShellDomain::Initialize()
{
    FEShellDomain::Initialize();
    FindSSI();
    return true;
}

//-----------------------------------------------------------------------------
void FESSIShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	FEShellDomain::PreSolveUpdate(timeInfo);
}

//-----------------------------------------------------------------------------
//! Find interfaces between solid element faces and shell elements
//!
void FESSIShellDomain::FindSSI()
{
	// find out if there are solid domains in this model
	vector<FESolidDomain*> psd;
	FEMesh* mesh = GetMesh();
	int ndom = mesh->Domains();
	for (int i = 0; i<ndom; ++i) {
		FEDomain& pdom = mesh->Domain(i);
		FESolidDomain* psdom = dynamic_cast<FESolidDomain*>(&pdom);
		if (psdom) psd.push_back(psdom);
	}
	size_t nsdom = psd.size();

	// if there are no solid domains we're done
	if (nsdom == 0) return;

	FEMesh& m = *GetMesh();
	int nelem = Elements();
	int nf[9], nn;
	vec3d g[3];

	// check all elements in this shell domain
	for (int i = 0; i<nelem; ++i) {
		FEShellElement& el = *dynamic_cast<FEShellElement*>(&ElementRef(i));

		// check all solid domains
		for (int k = 0; k<nsdom; ++k) {

			// check each solid element in this domain
			int nselem = psd[k]->Elements();
			for (int l = 0; l<nselem; ++l) {
				FEElement& sel = psd[k]->ElementRef(l);

				// check all faces of this solid element
				int nfaces = m.Faces(sel);
				for (int j = 0; j<nfaces; ++j) {
					nn = m.GetFace(sel, j, nf);

					bool found = false;
					if (nn == el.Nodes())
					{
						switch (nn)
						{
						case 3: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) found = true; break;
						case 4: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) found = true; break;
						case 6: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) found = true; break;
						case 7: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2])) found = true; break;
						case 8: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) found = true; break;
						case 9: if (el.HasNode(nf[0]) && el.HasNode(nf[1]) && el.HasNode(nf[2]) && el.HasNode(nf[3])) found = true; break;
						default:
							assert(false);
						}
					}
					if (found) {
						// check interface side
						// get outward normal to solid element face
						vec3d n0 = mesh->Node(nf[0]).m_r0, n1 = mesh->Node(nf[1]).m_r0, n2 = mesh->Node(nf[2]).m_r0;
						vec3d nsld = (n1 - n0) ^ (n2 - n1);
						// get outward normal to shell face
						CoBaseVectors0(el, 0, g);
						vec3d nshl = g[2];
						nshl.unit();
						// compare normals
						if (nsld*nshl > 0) {
							// store result
							sel.m_bitfc.resize(sel.Nodes(), false);
							for (int n = 0; n<nn; ++n) {
								int m = sel.FindNode(nf[n]);
								if (m > -1) sel.m_bitfc[m] = true;
							}
						}
					}
				}
			}
		}
	}
}


//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors0(FEShellElement& el, int n, vec3d g[3])
{
	int i;

	int neln = el.Nodes();

	// current nodal coordinates and directors
	vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
	for (i = 0; i<neln; ++i)
	{
		FENode& ni = m_pMesh->Node(el.m_node[i]);
		r[i] = ni.m_r0;
		D[i] = ni.m_d0;
	}

	double eta = el.gt(n);

	double* Mr = el.Hr(n);
	double* Ms = el.Hs(n);
	double* M = el.H(n);

	// initialize covariant basis vectors
	g[0] = g[1] = g[2] = vec3d(0, 0, 0);

	for (i = 0; i<neln; ++i)
	{
		g[0] += (r[i] - D[i] * (1 - eta) / 2)*Mr[i];
		g[1] += (r[i] - D[i] * (1 - eta) / 2)*Ms[i];
		g[2] += D[i] * (M[i] / 2);
	}
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors0(FEShellElement& el, int n, vec3d gcnt[3])
{
	vec3d gcov[3];
	CoBaseVectors0(el, n, gcov);

	mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
		gcov[0].y, gcov[1].y, gcov[2].y,
		gcov[0].z, gcov[1].z, gcov[2].z);
	mat3d Ji = J.inverse();

	gcnt[0] = vec3d(Ji(0, 0), Ji(0, 1), Ji(0, 2));
	gcnt[1] = vec3d(Ji(1, 0), Ji(1, 1), Ji(1, 2));
	gcnt[2] = vec3d(Ji(2, 0), Ji(2, 1), Ji(2, 2));

}

//-----------------------------------------------------------------------------
double FESSIShellDomain::invjac0(FEShellElement& el, double Ji[3][3], int n)
{
	vec3d gcov[3];
	CoBaseVectors0(el, n, gcov);

	mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
		gcov[0].y, gcov[1].y, gcov[2].y,
		gcov[0].z, gcov[1].z, gcov[2].z);

	double det = J.det();

	// make sure the determinant is positive
	if (det <= 0) throw NegativeJacobian(el.GetID(), n + 1, det);

	// calculate the inverse of the jacobian
	double deti = 1.0 / det;

	Ji[0][0] = deti*(J[1][1] * J[2][2] - J[1][2] * J[2][1]);
	Ji[1][0] = deti*(J[1][2] * J[2][0] - J[1][0] * J[2][2]);
	Ji[2][0] = deti*(J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	Ji[0][1] = deti*(J[0][2] * J[2][1] - J[0][1] * J[2][2]);
	Ji[1][1] = deti*(J[0][0] * J[2][2] - J[0][2] * J[2][0]);
	Ji[2][1] = deti*(J[0][1] * J[2][0] - J[0][0] * J[2][1]);

	Ji[0][2] = deti*(J[0][1] * J[1][2] - J[1][1] * J[0][2]);
	Ji[1][2] = deti*(J[0][2] * J[1][0] - J[0][0] * J[1][2]);
	Ji[2][2] = deti*(J[0][0] * J[1][1] - J[0][1] * J[1][0]);

	return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FESSIShellDomain::detJ0(FEShellElement &el, int n)
{
	vec3d gcov[3];
	CoBaseVectors0(el, n, gcov);

	mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
		gcov[0].y, gcov[1].y, gcov[2].y,
		gcov[0].z, gcov[1].z, gcov[2].z);
	return J.det();
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors(FEShellElement& el, int n, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rt;
        D[i] = ni.m_d0 + ni.get_vec3d(m_dofx, m_dofy, m_dofz) - ni.get_vec3d(m_dofu, m_dofv, m_dofw);
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
void FESSIShellDomain::ContraBaseVectors(FEShellElement& el, int n, vec3d gcnt[3])
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
double FESSIShellDomain::detJ(FEShellElement& el, int n)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
double FESSIShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
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
vec3d FESSIShellDomain::evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n)
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
double FESSIShellDomain::invjact(FEShellElement& el, double Ji[3][3], int n)
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
//! evaluate a scalar function over the shell
double FESSIShellDomain::evaluate(FEShellElement& el, double* pn, double* dpn, int n)
{
    const double *M;
    double eta;
    vec3d gcnt[3];
    double Mu, Md;
    double p = 0;
    
    eta = el.gt(n);
    
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        Mu = (1+eta)/2*M[i];
        Md = (1-eta)/2*M[i];
        p += Mu*pn[i] + Md*dpn[i];
    }
    
    return p;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of a scalar function over the shell
vec3d FESSIShellDomain::gradient(FEShellElement& el, double* pn, double* dpn, int n)
{
    const double* Mr, *Ms, *M;
    double eta;
    vec3d gcnt[3];
    vec3d gradM, gradMu, gradMd;
    vec3d gradp(0,0,0);
    
    eta = el.gt(n);
    
    Mr = el.Hr(n);
    Ms = el.Hs(n);
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
        gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
        gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        gradp += gradMu*pn[i] + gradMd*dpn[i];
    }
    
    return gradp;
}

//-----------------------------------------------------------------------------
//! evaluate a scalar function over the shell
double FESSIShellDomain::evaluate(FEShellElement& el, vector<double> pn, vector<double> dpn, int n)
{
    const double *M;
    double eta;
    vec3d gcnt[3];
    double Mu, Md;
    double p = 0;
    
    eta = el.gt(n);
    
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        Mu = (1+eta)/2*M[i];
        Md = (1-eta)/2*M[i];
        p += Mu*pn[i] + Md*dpn[i];
    }
    
    return p;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of a scalar function over the shell
vec3d FESSIShellDomain::gradient(FEShellElement& el, vector<double> pn, vector<double> dpn, int n)
{
    const double* Mr, *Ms, *M;
    double eta;
    vec3d gcnt[3];
    vec3d gradM, gradMu, gradMd;
    vec3d gradp(0,0,0);
    
    eta = el.gt(n);
    
    Mr = el.Hr(n);
    Ms = el.Hs(n);
    M  = el.H(n);
    
    ContraBaseVectors(el, n, gcnt);
    
    int neln = el.Nodes();
    
    for (int i=0; i<neln; ++i)
    {
        gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
        gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
        gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        gradp += gradMu*pn[i] + gradMd*dpn[i];
    }
    
    return gradp;
}
