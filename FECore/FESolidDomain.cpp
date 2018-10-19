#include "stdafx.h"
#include "FESolidDomain.h"
#include "FEMaterial.h"
#include "tools.h"

//-----------------------------------------------------------------------------
FESolidDomain::FESolidDomain(FEModel* pfem) : FEDomain(FE_DOMAIN_SOLID, pfem)
{
    m_dofx = pfem->GetDOFIndex("x");
    m_dofy = pfem->GetDOFIndex("y");
    m_dofz = pfem->GetDOFIndex("z");
    m_dofsx = pfem->GetDOFIndex("sx");
    m_dofsy = pfem->GetDOFIndex("sy");
    m_dofsz = pfem->GetDOFIndex("sz");
    m_dofsxp = pfem->GetDOFIndex("sxp");
    m_dofsyp = pfem->GetDOFIndex("syp");
    m_dofszp = pfem->GetDOFIndex("szp");
}

//-----------------------------------------------------------------------------
void FESolidDomain::Create(int nsize, int elemType)
{
	// allocate elements
    m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i)
	{
		FESolidElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	// set element type
    if (elemType != -1)
        for (int i=0; i<nsize; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
void FESolidDomain::CopyFrom(FEMeshPartition* pd)
{
    FESolidDomain* psd = dynamic_cast<FESolidDomain*>(pd);
    m_Elem = psd->m_Elem;
	for (int i=0; i<m_Elem.size(); ++i) m_Elem[i].SetMeshPartition(this);
}

//-----------------------------------------------------------------------------
//! initialize element data
bool FESolidDomain::Init()
{
	// base class first
	if (FEDomain::Init() == false) return false;

	// init solid element data
	double Ji[3][3];
	for (int i = 0; i<(int)m_Elem.size(); ++i)
	{
		FESolidElement& el = Element(i);
		int nint = el.GaussPoints();

		for (int n=0; n<nint; ++n)
		{
			invjac0(el, Ji, n);
			el.m_J0i[n] = mat3d(Ji);
		}
	}

	// nodal coordinates
	const int NELN = FEElement::MAX_NODES;
	vec3d r0[NELN], r[NELN], v[NELN], a[NELN];
	for (int i = 0; i < Elements(); ++i)
	{
		FESolidElement& el = Element(i);
		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			r0[j] = node.m_r0;
		}

		// loop over the integration points and calculate
		// the stress at the integration point
		int nint = el.GaussPoints();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// material point coordinates
			mp.m_r0 = el.Evaluate(r0, n);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
// Reset data
void FESolidDomain::Reset()
{
	for (int i = 0; i<(int)m_Elem.size(); ++i)
    {
        FESolidElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FEMaterialPoint* pt = el.GetMaterialPoint(j);
            if (pt) pt->Init();
        }
    }
}

//-----------------------------------------------------------------------------
//! This function finds the element in which point y lies and returns
//! the isoparametric coordinates in r if an element is found
//! (This has only been implemeneted for hexes!)
FESolidElement* FESolidDomain::FindElement(const vec3d& y, double r[3])
{
    int i, j;
    int NE = Elements();
    vec3d x[FEElement::MAX_NODES];
    for (i=0; i<NE; ++i)
    {
        // get the next element
        FESolidElement& e = Element(i);
        assert(e.Type() == FE_HEX8G8);
        
        // get the element nodal coordinates
        int neln = e.Nodes();
        for (j=0; j<neln; ++j) x[j] = m_pMesh->Node(e.m_node[j]).m_rt;
        
        // first, as a quick check, we see if y lies in the bounding box defined by x
        FEBoundingBox box(x[0]);
        for (j=1; j<neln; ++j) box.add(x[j]);
        
        if (box.IsInside(y))
        {
			// If the point y lies inside the box, we apply a Newton method to find
			// the isoparametric coordinates r
			ProjectToElement(e, y, r);
            
            // see if the point r lies inside the element
            const double eps = 1.0001;
            if ((r[0] >= -eps) && (r[0] <= eps) &&
                (r[1] >= -eps) && (r[1] <= eps) &&
                (r[2] >= -eps) && (r[2] <= eps)) return &e;
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------
void FESolidDomain::ProjectToElement(FESolidElement& el, const vec3d& p, double r[3])
{
	const int MN = FEElement::MAX_NODES;
	vec3d rt[MN];

	// get the element nodal coordinates
	int ne = el.Nodes();
	for (int i = 0; i<ne; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;

    r[0] = r[1] = r[2] = 0;
    const double tol = 1e-5;
    double dr[3], norm;
	double H[MN], Gr[MN], Gs[MN], Gt[MN];
    do
    {
		// evaluate shape functions
		el.shape_fnc(H, r[0], r[1], r[2]);
                
		// evaluate shape function derivatives
		el.shape_deriv(Gr, Gs, Gt, r[0], r[1], r[2]);
                
		// solve for coordinate increment
        double R[3] = {0}, A[3][3] = {0};
        for (int i=0; i<ne; ++i)
        {
            R[0] += rt[i].x*H[i];
            R[1] += rt[i].y*H[i];
            R[2] += rt[i].z*H[i];
                    
            A[0][0] -= rt[i].x*Gr[i]; A[0][1] -= rt[i].x*Gs[i]; A[0][2] -= rt[i].x*Gt[i];
            A[1][0] -= rt[i].y*Gr[i]; A[1][1] -= rt[i].y*Gs[i]; A[1][2] -= rt[i].y*Gt[i];
            A[2][0] -= rt[i].z*Gr[i]; A[2][1] -= rt[i].z*Gs[i]; A[2][2] -= rt[i].z*Gt[i];
        }
        R[0] = p.x - R[0];
        R[1] = p.y - R[1];
        R[2] = p.z - R[2];
                
        solve_3x3(A, R, dr);
        r[0] -= dr[0];
        r[1] -= dr[1];
        r[2] -= dr[2];
                
        norm = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    }
    while (norm > tol);	
}

//-----------------------------------------------------------------------------
void FESolidDomain::ProjectToReferenceElement(FESolidElement& el, const vec3d& p, double r[3])
{
	const int MN = FEElement::MAX_NODES;
	vec3d rt[MN];

	// get the element nodal coordinates
	int ne = el.Nodes();
	for (int i = 0; i<ne; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_r0;

    r[0] = r[1] = r[2] = 0;
    const double tol = 1e-5;
    double dr[3], norm;
	double H[MN], Gr[MN], Gs[MN], Gt[MN];
    do
    {
		// evaluate shape functions
		el.shape_fnc(H, r[0], r[1], r[2]);
                
		// evaluate shape function derivatives
		el.shape_deriv(Gr, Gs, Gt, r[0], r[1], r[2]);
                
		// solve for coordinate increment
        double R[3] = {0}, A[3][3] = {0};
        for (int i=0; i<ne; ++i)
        {
            R[0] += rt[i].x*H[i];
            R[1] += rt[i].y*H[i];
            R[2] += rt[i].z*H[i];
                    
            A[0][0] -= rt[i].x*Gr[i]; A[0][1] -= rt[i].x*Gs[i]; A[0][2] -= rt[i].x*Gt[i];
            A[1][0] -= rt[i].y*Gr[i]; A[1][1] -= rt[i].y*Gs[i]; A[1][2] -= rt[i].y*Gt[i];
            A[2][0] -= rt[i].z*Gr[i]; A[2][1] -= rt[i].z*Gs[i]; A[2][2] -= rt[i].z*Gt[i];
        }
        R[0] = p.x - R[0];
        R[1] = p.y - R[1];
        R[2] = p.z - R[2];
                
        solve_3x3(A, R, dr);
        r[0] -= dr[0];
        r[1] -= dr[1];
        r[2] -= dr[2];
                
        norm = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    }
    while (norm > tol);	
}


//-----------------------------------------------------------------------------
//! get the current nodal coordinates
void FESolidDomain::GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt)
{
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i) rt[i] = m_pMesh->Node(el.m_node[i]).m_rt;

	// check for solid-shell interface nodes
	if (el.m_bitfc.empty() == false)
	{
		for (int i = 0; i<neln; ++i)
		{
			if (el.m_bitfc[i])
			{
				FENode& nd = m_pMesh->Node(el.m_node[i]);
				rt[i] -= nd.m_d0 + nd.get_vec3d(m_dofx, m_dofy, m_dofz) - nd.get_vec3d(m_dofsx, m_dofsy, m_dofsz);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! get the current nodal coordinates
void FESolidDomain::GetCurrentNodalCoordinates(const FESolidElement& el, vec3d* rt, double alpha)
{
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i) {
		FENode& nd = m_pMesh->Node(el.m_node[i]);
		rt[i] = nd.m_rt*alpha + nd.m_rp*(1 - alpha);
	}

	// check for solid-shell interface nodes
	if (el.m_bitfc.empty() == false)
	{
		for (int i = 0; i<neln; ++i)
		{
			if (el.m_bitfc[i]) {
				FENode& nd = m_pMesh->Node(el.m_node[i]);
				rt[i] -= nd.m_d0 + rt[i] - nd.m_r0
					- nd.get_vec3d(m_dofsx, m_dofsy, m_dofsz)*alpha
					- nd.get_vec3d(m_dofsxp, m_dofsyp, m_dofszp)*(1 - alpha);
			}
		}
	}

}

//-----------------------------------------------------------------------------
//! get the reference nodal coordinates
void FESolidDomain::GetReferenceNodalCoordinates(const FESolidElement& el, vec3d* r0)
{
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i) r0[i] = m_pMesh->Node(el.m_node[i]).m_r0;

	// check for solid-shell interface nodes
	if (el.m_bitfc.empty() == false)
	{
		for (int i = 0; i<neln; ++i)
		{
			if (el.m_bitfc[i])
			{
				FENode& nd = m_pMesh->Node(el.m_node[i]);
				r0[i] -= nd.m_d0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! get the previous nodal coordinates
void FESolidDomain::GetPreviousNodalCoordinates(const FESolidElement& el, vec3d* rp)
{
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i) rp[i] = m_pMesh->Node(el.m_node[i]).m_rp;

	// check for solid-shell interface nodes
	if (el.m_bitfc.empty() == false)
	{
		for (int i = 0; i<neln; ++i)
		{
			if (el.m_bitfc[i])
			{
				FENode& nd = m_pMesh->Node(el.m_node[i]);
				rp[i] -= nd.m_d0 + nd.m_rp - nd.m_r0
					- nd.get_vec3d(m_dofsxp, m_dofsyp, m_dofszp);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element el at integration point n.
//! The deformation gradient is returned in F and its determinant is the return
//! value of the function
double FESolidDomain::defgrad(FESolidElement &el, mat3d &F, int n)
{
    // nodal points
    vec3d r[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, r);
    
    // calculate inverse jacobian
//    double Ji[3][3];
//    invjac0(el, Ji, n);
	mat3d& Ji = el.m_J0i[n];

	// shape function derivatives
	double *Grn = el.Gr(n);
	double *Gsn = el.Gs(n);
	double *Gtn = el.Gt(n);

    // calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
    F[1][0] = F[1][1] = F[1][2] = 0;
    F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        double Gri = Grn[i];
        double Gsi = Gsn[i];
        double Gti = Gtn[i];
        
        double x = r[i].x;
        double y = r[i].y;
        double z = r[i].z;
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        double GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
        double GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
        double GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
        
        // calculate deformation gradient F
        F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
        F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
        F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
    }
    
    double D = F.det();
    if (D <= 0) throw NegativeJacobian(el.GetID(), n, D, &el);
    
    return D;
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element at point r,s,t
double FESolidDomain::defgrad(FESolidElement &el, mat3d &F, double r, double s, double t)
{
    // shape function derivatives
    const int NME = FEElement::MAX_NODES;
    double Gr[NME], Gs[NME], Gt[NME];
    el.shape_deriv(Gr, Gs, Gt, r, s, t);
    
    // nodal points
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt);

    // calculate inverse jacobian
    double Ji[3][3];
    invjac0(el, Ji, r, s, t);
    
    // calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
    F[1][0] = F[1][1] = F[1][2] = 0;
    F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        double Gri = Gr[i];
        double Gsi = Gs[i];
        double Gti = Gt[i];
        
        double x = rt[i].x;
        double y = rt[i].y;
        double z = rt[i].z;
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        double GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
        double GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
        double GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
        
        // calculate deformation gradient F
        F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
        F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
        F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
    }
    
    double D = F.det();
    if (D <= 0) throw NegativeJacobian(el.GetID(), -1, D, &el);
    
    return D;
}

//-----------------------------------------------------------------------------
//! Calculate the deformation gradient of element el at integration point n.
//! The deformation gradient is returned in F and its determinant is the return
//! value of the function
double FESolidDomain::defgradp(FESolidElement &el, mat3d &F, int n)
{
    // nodal coordinates
    vec3d r[FEElement::MAX_NODES];
	GetPreviousNodalCoordinates(el, r);
    
    // calculate inverse jacobian
    double Ji[3][3];
    invjac0(el, Ji, n);

	// shape function derivatives
	double *Grn = el.Gr(n);
	double *Gsn = el.Gs(n);
	double *Gtn = el.Gt(n);

    // calculate deformation gradient
    F[0][0] = F[0][1] = F[0][2] = 0;
    F[1][0] = F[1][1] = F[1][2] = 0;
    F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        double Gri = Grn[i];
        double Gsi = Gsn[i];
        double Gti = Gtn[i];
        
        double x = r[i].x;
        double y = r[i].y;
        double z = r[i].z;
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        double GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
        double GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
        double GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
        
        // calculate deformation gradient F
        F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
        F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
        F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
    }
    
    double D = F.det();
    if (D <= 0) throw NegativeJacobian(el.GetID(), n, D, &el);
    
    return D;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjac0(const FESolidElement& el, double Ji[3][3], int n)
{
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
	GetReferenceNodalCoordinates(el, r0);
   
    // calculate Jacobian
    double J[3][3] = {0};
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = el.Gr(n)[i];
        const double& Gsi = el.Gs(n)[i];
        const double& Gti = el.Gt(n)[i];
        
        const double& x = r0[i].x;
        const double& y = r0[i].y;
        const double& z = r0[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate the inverse jacobian
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
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjac0(const FESolidElement& el, double Ji[3][3], double r, double s, double t)
{
    // nodal coordinates
    const int NMAX = FEElement::MAX_NODES;
    vec3d r0[NMAX];
	GetReferenceNodalCoordinates(el, r0);
    
    // evaluate shape function derivatives
    double Gr[NMAX], Gs[NMAX], Gt[NMAX];
    el.shape_deriv(Gr, Gs, Gt, r, s, t);
    
    // calculate Jacobian
    double J[3][3] = {0};
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = Gr[i];
        const double& Gsi = Gs[i];
        const double& Gti = Gt[i];
        
        const double& x = r0[i].x;
        const double& y = r0[i].y;
        const double& z = r0[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), -1, det);
    
    // calculate the inverse jacobian
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
double FESolidDomain::invjac0(const FESolidElement& el, double r, double s, double t, mat3d& J)
{
    double Ji[3][3];
    double J0 = invjac0(el, Ji, r, s, t);
    J = mat3d(Ji[0][0], Ji[0][1], Ji[0][2], Ji[1][0], Ji[1][1], Ji[1][2], Ji[2][0], Ji[2][1], Ji[2][2]);
    return J0;
}

//-----------------------------------------------------------------------------
//! Calculate the inverse jacobian with respect to the current frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjact(FESolidElement& el, double Ji[3][3], int n)
{
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt);

    // calculate jacobian
    double J[3][3] = {0};
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = el.Gr(n)[i];
        const double& Gsi = el.Gs(n)[i];
        const double& Gti = el.Gt(n)[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate inverse jacobian
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
//! Calculate the inverse jacobian with respect to the current frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjact(FESolidElement& el, double Ji[3][3], int n, const double alpha)
{
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt, alpha);
    
    // calculate jacobian
	int neln = el.Nodes();
	double J[3][3] = { 0 };
    for (int i=0; i<neln; ++i)
    {
        const double& Gri = el.Gr(n)[i];
        const double& Gsi = el.Gs(n)[i];
        const double& Gti = el.Gt(n)[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate inverse jacobian
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
//! Calculate the inverse jacobian with respect to the previous time frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjactp(FESolidElement& el, double Ji[3][3], int n)
{
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetPreviousNodalCoordinates(el, rt);
    
    // calculate jacobian
	int neln = el.Nodes();
	double J[3][3] = { 0 };
    for (int i=0; i<neln; ++i)
    {
        const double& Gri = el.Gr(n)[i];
        const double& Gsi = el.Gs(n)[i];
        const double& Gti = el.Gt(n)[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), n+1, det);
    
    // calculate inverse jacobian
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
//! Calculate the inverse jacobian with respect to the reference frame at
//! integration point n. The inverse jacobian is retured in Ji
//! The return value is the determinant of the Jacobian (not the inverse!)
double FESolidDomain::invjact(FESolidElement& el, double Ji[3][3], double r, double s, double t)
{
    // nodal coordinates
    const int NMAX = FEElement::MAX_NODES;
    vec3d rt[NMAX];
	GetCurrentNodalCoordinates(el, rt);

    // evaluate shape function derivatives
    double Gr[NMAX], Gs[NMAX], Gt[NMAX];
    el.shape_deriv(Gr, Gs, Gt, r, s, t);
    
    // calculate Jacobian
    double J[3][3] = {0};
	const int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = Gr[i];
        const double& Gsi = Gs[i];
        const double& Gti = Gt[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), -1, det);
    
    // calculate the inverse jacobian
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
//! calculate gradient of function at integration points
vec3d FESolidDomain::gradient(FESolidElement& el, double* fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    double Gx, Gy, Gz;
    
    vec3d gradf;
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i]+Ji[2][0]*Gtn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i]+Ji[2][1]*Gtn[i];
        Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i]+Ji[2][2]*Gtn[i];
        
        // calculate pressure gradient
        gradf.x += Gx*fn[i];
        gradf.y += Gy*fn[i];
        gradf.z += Gz*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate gradient of function at integration points
vec3d FESolidDomain::gradient(FESolidElement& el, vector<double>& fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    double Gx, Gy, Gz;
    
    vec3d gradf;
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i]+Ji[2][0]*Gtn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i]+Ji[2][1]*Gtn[i];
        Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i]+Ji[2][2]*Gtn[i];
        
        // calculate pressure gradient
        gradf.x += Gx*fn[i];
        gradf.y += Gy*fn[i];
        gradf.z += Gz*fn[i];
    }
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
mat3d FESolidDomain::gradient(FESolidElement& el, vec3d* fn, int n)
{
    double Ji[3][3];
    invjact(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
    vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
    vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
    
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    mat3d gradf;
    gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gradf += fn[i] & (g1*Gr[i] + g2*Gs[i] + g3*Gt[i]);
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate spatial gradient of function at integration points
//! at previous time
mat3d FESolidDomain::gradientp(FESolidElement& el, vec3d* fn, int n)
{
    double Ji[3][3];
    invjactp(el, Ji, n);
				
    vec3d g1(Ji[0][0],Ji[0][1],Ji[0][2]);
    vec3d g2(Ji[1][0],Ji[1][1],Ji[1][2]);
    vec3d g3(Ji[2][0],Ji[2][1],Ji[2][2]);
    
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    mat3d gradf;
    gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gradf += fn[i] & (g1*Gr[i] + g2*Gs[i] + g3*Gt[i]);
    
    return gradf;
}

//-----------------------------------------------------------------------------
//! calculate material gradient of function at integration points
vec3d FESolidDomain::Gradient(FESolidElement& el, double* fn, int n)
{
    double Ji[3][3];
    invjac0(el, Ji, n);
    
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    double Gx, Gy, Gz;
    
    vec3d Gradf;
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        Gx = Ji[0][0]*Grn[i]+Ji[1][0]*Gsn[i]+Ji[2][0]*Gtn[i];
        Gy = Ji[0][1]*Grn[i]+Ji[1][1]*Gsn[i]+Ji[2][1]*Gtn[i];
        Gz = Ji[0][2]*Grn[i]+Ji[1][2]*Gsn[i]+Ji[2][2]*Gtn[i];
        
        // calculate pressure gradient
        Gradf.x += Gx*fn[i];
        Gradf.y += Gy*fn[i];
        Gradf.z += Gz*fn[i];
    }
    
    return Gradf;
}

//-----------------------------------------------------------------------------
//! calculate material gradient of function at integration points
mat3d FESolidDomain::Gradient(FESolidElement& el, vec3d* fn, int n)
{
    vec3d Gcnt[3];
    ContraBaseVectors0(el, n, Gcnt);
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    mat3d Gradf;
    Gradf.zero();
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        Gradf += fn[i] & (Gcnt[0]*Gr[i] + Gcnt[1]*Gs[i] + Gcnt[2]*Gt[i]);
    
    return Gradf;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to current frame
double FESolidDomain::detJt(FESolidElement &el, int n)
{
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt);

    // shape function derivatives
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    // jacobian matrix
    double J[3][3] = {0};
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = Grn[i];
        const double& Gsi = Gsn[i];
        const double& Gti = Gtn[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to current frame
double FESolidDomain::detJt(FESolidElement &el, int n, const double alpha)
{
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt, alpha);

    // shape function derivatives
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    // jacobian matrix
	int neln = el.Nodes();
	double J[3][3] = { 0 };
    for (int i=0; i<neln; ++i)
    {
        const double& Gri = Grn[i];
        const double& Gsi = Gsn[i];
        const double& Gti = Gtn[i];
        
        const double& x = rt[i].x;
        const double& y = rt[i].y;
        const double& z = rt[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
    + J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
    + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
//! Calculate jacobian with respect to reference frame
double FESolidDomain::detJ0(FESolidElement &el, int n)
{
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
	GetReferenceNodalCoordinates(el, r0);

    // shape function derivatives
    double* Grn = el.Gr(n);
    double* Gsn = el.Gs(n);
    double* Gtn = el.Gt(n);
    
    // jacobian matrix
    double J[3][3] = {0};
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
    {
        const double& Gri = Grn[i];
        const double& Gsi = Gsn[i];
        const double& Gti = Gtn[i];
        
        const double& x = r0[i].x;
        const double& y = r0[i].y;
        const double& z = r0[i].z;
        
        J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
        J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
        J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
    }
    
    // calculate the determinant
    double det =  J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
				+ J[0][1]*(J[1][2]*J[2][0] - J[2][2]*J[1][0])
				+ J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    
    return det;
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant basis vectors of a solid element
//! in reference configuration at an integration point

void FESolidDomain::CoBaseVectors0(FESolidElement& el, int j, vec3d g[3])
{
    // get the shape function derivatives
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    double* Ht = el.Gt(j);
    
    // nodal coordinates
    vec3d r0[FEElement::MAX_NODES];
	GetReferenceNodalCoordinates(el, r0);
    
    g[0] = g[1] = g[2] = vec3d(0,0,0);
	int n = el.Nodes();
	for (int i = 0; i<n; ++i)
    {
        g[0] += r0[i]*Hr[i];
        g[1] += r0[i]*Hs[i];
        g[2] += r0[i]*Ht[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the covariant basis vectors of a solid element
//! at an integration point

void FESolidDomain::CoBaseVectors(FESolidElement& el, int j, vec3d g[3])
{
    // get the shape function derivatives
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    double* Ht = el.Gt(j);
    
    // nodal coordinates
    vec3d rt[FEElement::MAX_NODES];
	GetCurrentNodalCoordinates(el, rt);

    g[0] = g[1] = g[2] = vec3d(0,0,0);
	int n = el.Nodes();
	for (int i = 0; i<n; ++i)
    {
        g[0] += rt[i]*Hr[i];
        g[1] += rt[i]*Hs[i];
        g[2] += rt[i]*Ht[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the contravariant basis vectors in ref config
//! of a solid element at an integration point

void FESolidDomain::ContraBaseVectors0(FESolidElement& el, int j, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors0(el, j, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
}

//-----------------------------------------------------------------------------
//! This function calculates the contravariant basis vectors of a solid element
//! at an integration point

void FESolidDomain::ContraBaseVectors(FESolidElement& el, int j, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors(el, j, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of covariant basis
//! vectors of a solid element at an integration point

void FESolidDomain::CoBaseVectorDerivatives(FESolidElement& el, int j, vec3d dg[3][3])
{
    FEMesh& m = *m_pMesh;
    
    // get the nr of nodes
    int n = el.Nodes();
    
    // get the shape function derivatives
    double* Hrr = el.Grr(j); double* Hrs = el.Grs(j); double* Hrt = el.Grt(j);
    double* Hsr = el.Gsr(j); double* Hss = el.Gss(j); double* Hst = el.Gst(j);
    double* Htr = el.Gtr(j); double* Hts = el.Gts(j); double* Htt = el.Gtt(j);
    
    dg[0][0] = dg[0][1] = dg[0][2] = vec3d(0,0,0);  // derivatives of g[0]
    dg[1][0] = dg[1][1] = dg[1][2] = vec3d(0,0,0);  // derivatives of g[1]
    dg[2][0] = dg[2][1] = dg[2][2] = vec3d(0,0,0);  // derivatives of g[2]
    
    for (int i=0; i<n; ++i)
    {
        vec3d rt = m.Node(el.m_node[i]).m_rt;
        dg[0][0] += rt*Hrr[i]; dg[0][1] += rt*Hsr[i]; dg[0][2] += rt*Htr[i];
        dg[1][0] += rt*Hrs[i]; dg[1][1] += rt*Hss[i]; dg[1][2] += rt*Hts[i];
        dg[2][0] += rt*Hrt[i]; dg[2][1] += rt*Hst[i]; dg[2][2] += rt*Htt[i];
    }
}

//-----------------------------------------------------------------------------
//! This function calculates the parametric derivatives of contravariant basis
//! vectors of a solid element at an integration point

void FESolidDomain::ContraBaseVectorDerivatives(FESolidElement& el, int j, vec3d dgcnt[3][3])
{
    vec3d gcnt[3];
    vec3d dgcov[3][3];
    ContraBaseVectors(el, j, gcnt);
    CoBaseVectorDerivatives(el, j, dgcov);
    
    // derivatives of gcnt[0]
    dgcnt[0][0] = -gcnt[0]*(gcnt[0]*dgcov[0][0])-gcnt[1]*(gcnt[0]*dgcov[1][0])-gcnt[2]*(gcnt[0]*dgcov[2][0]);
    dgcnt[0][1] = -gcnt[0]*(gcnt[0]*dgcov[0][1])-gcnt[1]*(gcnt[0]*dgcov[1][1])-gcnt[2]*(gcnt[0]*dgcov[2][1]);
    dgcnt[0][2] = -gcnt[0]*(gcnt[0]*dgcov[0][2])-gcnt[1]*(gcnt[0]*dgcov[1][2])-gcnt[2]*(gcnt[0]*dgcov[2][2]);
    
    // derivatives of gcnt[1]
    dgcnt[1][0] = -gcnt[0]*(gcnt[1]*dgcov[0][0])-gcnt[1]*(gcnt[1]*dgcov[1][0])-gcnt[2]*(gcnt[1]*dgcov[2][0]);
    dgcnt[1][1] = -gcnt[0]*(gcnt[1]*dgcov[0][1])-gcnt[1]*(gcnt[1]*dgcov[1][1])-gcnt[2]*(gcnt[1]*dgcov[2][1]);
    dgcnt[1][2] = -gcnt[0]*(gcnt[1]*dgcov[0][2])-gcnt[1]*(gcnt[1]*dgcov[1][2])-gcnt[2]*(gcnt[1]*dgcov[2][2]);
    
    // derivatives of gcnt[2]
    dgcnt[2][0] = -gcnt[0]*(gcnt[2]*dgcov[0][0])-gcnt[1]*(gcnt[2]*dgcov[1][0])-gcnt[2]*(gcnt[2]*dgcov[2][0]);
    dgcnt[2][1] = -gcnt[0]*(gcnt[2]*dgcov[0][1])-gcnt[1]*(gcnt[2]*dgcov[1][1])-gcnt[2]*(gcnt[2]*dgcov[2][1]);
    dgcnt[2][2] = -gcnt[0]*(gcnt[2]*dgcov[0][2])-gcnt[1]*(gcnt[2]*dgcov[1][2])-gcnt[2]*(gcnt[2]*dgcov[2][2]);
}

//-----------------------------------------------------------------------------
//! calculate the laplacian of a vector function at an integration point
vec3d FESolidDomain::lapvec(FESolidElement& el, vec3d* fn, int n)
{
    vec3d gcnt[3];
    vec3d dgcnt[3][3];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(n); double* Hrs = el.Grs(n); double* Hrt = el.Grt(n);
    double* Hsr = el.Gsr(n); double* Hss = el.Gss(n); double* Hst = el.Gst(n);
    double* Htr = el.Gtr(n); double* Hts = el.Gts(n); double* Htt = el.Gtt(n);
    
    vec3d lapv(0,0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        lapv += fn[i]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*gcnt[0]+
                       (gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*gcnt[1]+
                       (gcnt[0]*Hrt[i]+dgcnt[0][2]*Gr[i])*gcnt[2]+
                       (gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*gcnt[0]+
                       (gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*gcnt[1]+
                       (gcnt[1]*Hst[i]+dgcnt[1][2]*Gs[i])*gcnt[2]+
                       (gcnt[2]*Htr[i]+dgcnt[2][0]*Gt[i])*gcnt[0]+
                       (gcnt[2]*Hts[i]+dgcnt[2][1]*Gt[i])*gcnt[1]+
                       (gcnt[2]*Htt[i]+dgcnt[2][2]*Gt[i])*gcnt[2]);
    
    return lapv;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of the divergence of a vector function at an integration point
vec3d FESolidDomain::gradivec(FESolidElement& el, vec3d* fn, int n)
{
    vec3d gcnt[3];
    vec3d dgcnt[3][3];
    ContraBaseVectors(el, n, gcnt);
    ContraBaseVectorDerivatives(el, n, dgcnt);
				
    double* Gr = el.Gr(n);
    double* Gs = el.Gs(n);
    double* Gt = el.Gt(n);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(n); double* Hrs = el.Grs(n); double* Hrt = el.Grt(n);
    double* Hsr = el.Gsr(n); double* Hss = el.Gss(n); double* Hst = el.Gst(n);
    double* Htr = el.Gtr(n); double* Hts = el.Gts(n); double* Htt = el.Gtt(n);
    
    vec3d gdv(0,0,0);
    int N = el.Nodes();
    for (int i=0; i<N; ++i)
        gdv +=
        gcnt[0]*((gcnt[0]*Hrr[i]+dgcnt[0][0]*Gr[i])*fn[i])+
        gcnt[1]*((gcnt[0]*Hrs[i]+dgcnt[0][1]*Gr[i])*fn[i])+
        gcnt[2]*((gcnt[0]*Hrt[i]+dgcnt[0][2]*Gr[i])*fn[i])+
        gcnt[0]*((gcnt[1]*Hsr[i]+dgcnt[1][0]*Gs[i])*fn[i])+
        gcnt[1]*((gcnt[1]*Hss[i]+dgcnt[1][1]*Gs[i])*fn[i])+
        gcnt[2]*((gcnt[1]*Hst[i]+dgcnt[1][2]*Gs[i])*fn[i])+
        gcnt[0]*((gcnt[2]*Htr[i]+dgcnt[2][0]*Gt[i])*fn[i])+
        gcnt[1]*((gcnt[2]*Hts[i]+dgcnt[2][1]*Gt[i])*fn[i])+
        gcnt[2]*((gcnt[2]*Htt[i]+dgcnt[2][2]*Gt[i])*fn[i]);
    
    return gdv;
}

//-----------------------------------------------------------------------------
//! This function calculates the transpose of the spatial gradient of the shape
//! function spatial gradients, gradT(grad H), at an integration point

void FESolidDomain::gradTgradShape(FESolidElement& el, int j, vector<mat3d>& mn)
{
    int N = el.Nodes();
    vec3d g[3], dg[3][3];
    ContraBaseVectors(el, j, g);
    ContraBaseVectorDerivatives(el, j, dg);
    
    // get the shape functions
    double* Hr = el.Gr(j);
    double* Hs = el.Gs(j);
    double* Ht = el.Gt(j);
    
    // get the shape function derivatives
    double* Hrr = el.Grr(j); double* Hrs = el.Grs(j); double* Hrt = el.Grt(j);
    double* Hsr = el.Gsr(j); double* Hss = el.Gss(j); double* Hst = el.Gst(j);
    double* Htr = el.Gtr(j); double* Hts = el.Gts(j); double* Htt = el.Gtt(j);
    
    for (int i=0; i<N; ++i) {
        mn[i] = ((g[0] & dg[0][0]) + (g[1] & dg[0][1]) + (g[2] & dg[0][2]))*Hr[i]
        + ((g[0] & dg[1][0]) + (g[1] & dg[1][1]) + (g[2] & dg[1][2]))*Hs[i]
        + ((g[0] & dg[2][0]) + (g[1] & dg[2][1]) + (g[2] & dg[2][2]))*Ht[i]
        + (g[0] & g[0])*Hrr[i] + (g[1] & g[0])*Hsr[i] + (g[2] & g[0])*Htr[i]
        + (g[0] & g[1])*Hrs[i] + (g[1] & g[1])*Hss[i] + (g[2] & g[1])*Hts[i]
        + (g[0] & g[2])*Hrt[i] + (g[1] & g[2])*Hst[i] + (g[2] & g[2])*Htt[i];
    }
}

//-----------------------------------------------------------------------------
double FESolidDomain::ShapeGradient(FESolidElement& el, int n, vec3d* GradH)
{
    // calculate jacobian
    double Ji[3][3];
    double detJt = invjact(el, Ji, n);
    
    // evaluate shape function derivatives
    int ne = el.Nodes();
    for (int i = 0; i<ne; ++i)
    {
        double Gr = el.Gr(n)[i];
        double Gs = el.Gs(n)[i];
        double Gt = el.Gt(n)[i];
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        GradH[i].x = Ji[0][0] * Gr + Ji[1][0] * Gs + Ji[2][0] * Gt;
        GradH[i].y = Ji[0][1] * Gr + Ji[1][1] * Gs + Ji[2][1] * Gt;
        GradH[i].z = Ji[0][2] * Gr + Ji[1][2] * Gs + Ji[2][2] * Gt;
    }
    
    return detJt;
}

//-----------------------------------------------------------------------------
double FESolidDomain::ShapeGradient(FESolidElement& el, int n, vec3d* GradH, const double alpha)
{
    // calculate jacobian
    double Ji[3][3];
    double detJt = invjact(el, Ji, n, alpha);
    
    // evaluate shape function derivatives
    int ne = el.Nodes();
    for (int i = 0; i<ne; ++i)
    {
        double Gr = el.Gr(n)[i];
        double Gs = el.Gs(n)[i];
        double Gt = el.Gt(n)[i];
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        GradH[i].x = Ji[0][0] * Gr + Ji[1][0] * Gs + Ji[2][0] * Gt;
        GradH[i].y = Ji[0][1] * Gr + Ji[1][1] * Gs + Ji[2][1] * Gt;
        GradH[i].z = Ji[0][2] * Gr + Ji[1][2] * Gs + Ji[2][2] * Gt;
    }
    
    return detJt;
}

//-----------------------------------------------------------------------------
double FESolidDomain::ShapeGradient0(FESolidElement& el, int n, vec3d* GradH)
{
    // calculate jacobian
    double Ji[3][3];
    double detJ0 = invjac0(el, Ji, n);
    
    // evaluate shape function derivatives
    int ne = el.Nodes();
    for (int i = 0; i<ne; ++i)
    {
        double Gr = el.Gr(n)[i];
        double Gs = el.Gs(n)[i];
        double Gt = el.Gt(n)[i];
        
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        GradH[i].x = Ji[0][0] * Gr + Ji[1][0] * Gs + Ji[2][0] * Gt;
        GradH[i].y = Ji[0][1] * Gr + Ji[1][1] * Gs + Ji[2][1] * Gt;
        GradH[i].z = Ji[0][2] * Gr + Ji[1][2] * Gs + Ji[2][2] * Gt;
    }
    
    return detJ0;
}

//-----------------------------------------------------------------------------
double FESolidDomain::ShapeGradient(FESolidElement& el, double r, double s, double t, vec3d* GradH)
{
    // calculate jacobian
    double Ji[3][3];
    double detJt = invjact(el, Ji, r, s, t);
    
    // shape function derivatives
    double Gr[FEElement::MAX_NODES];
    double Gs[FEElement::MAX_NODES];
    double Gt[FEElement::MAX_NODES];
    el.shape_deriv(Gr, Gs, Gt, r, s, t);
    
    int neln = el.Nodes();
    for (int i = 0; i<neln; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        GradH[i].x = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
        GradH[i].y = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
        GradH[i].z = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];
    }
    
    return detJt;
}

//-----------------------------------------------------------------------------
double FESolidDomain::ShapeGradient0(FESolidElement& el, double r, double s, double t, vec3d* GradH)
{
    // calculate jacobian
    double Ji[3][3];
    double detJ0 = invjac0(el, Ji, r, s, t);
    
    // shape function derivatives
    double Gr[FEElement::MAX_NODES];
    double Gs[FEElement::MAX_NODES];
    double Gt[FEElement::MAX_NODES];
    el.shape_deriv(Gr, Gs, Gt, r, s, t);
    
    int neln = el.Nodes();
    for (int i = 0; i<neln; ++i)
    {
        // calculate global gradient of shape functions
        // note that we need the transposed of Ji, not Ji itself !
        GradH[i].x = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
        GradH[i].y = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
        GradH[i].z = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];
    }
    
    return detJ0;
}

//-----------------------------------------------------------------------------
//! calculate the volume of an element
double FESolidDomain::Volume(FESolidElement& el)
{
	vec3d r0[FEElement::MAX_NODES];

	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i) r0[i] = Node(el.m_lnode[i]).m_r0;

	int nint = el.GaussPoints();
	double *w = el.GaussWeights();
	double V = 0;
	for (int n = 0; n<nint; ++n)
	{
		// shape function derivatives
		double* Grn = el.Gr(n);
		double* Gsn = el.Gs(n);
		double* Gtn = el.Gt(n);

		// jacobian matrix
		double J[3][3] = { 0 };
		for (int i = 0; i<neln; ++i)
		{
			const double& Gri = Grn[i];
			const double& Gsi = Gsn[i];
			const double& Gti = Gtn[i];

			const double& x = r0[i].x;
			const double& y = r0[i].y;
			const double& z = r0[i].z;

			J[0][0] += Gri*x; J[0][1] += Gsi*x; J[0][2] += Gti*x;
			J[1][0] += Gri*y; J[1][1] += Gsi*y; J[1][2] += Gti*y;
			J[2][0] += Gri*z; J[2][1] += Gsi*z; J[2][2] += Gti*z;
		}

		// calculate the determinant
		double detJ0 = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
			+ J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
			+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

		V += detJ0*w[n];
	}

	return V;
}
