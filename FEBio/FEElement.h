// FEElement.h: interface for the FEElement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_)
#define AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElementLibrary.h"
#include "FEElementTraits.h"
#include "FEMaterialPoint.h"
#include "FE_enum.h"
#include "FEException.h"

//-----------------------------------------------------------------------------
//! Base class for all element classes

//! From this class the different element classes are derived.

class FEElement
{
public:
	// default constructor
	FEElement() : m_pT(0) 
	{ 
		static int n = 1;
		m_nID = n++;
		m_gid = -1;
		m_nrigid = -1; 
	}

	virtual ~FEElement() {}

	bool IsRigid() { return (m_nrigid >= 0); }

	// Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits)
	{
		m_pT = ptraits;

		int nint = GaussPoints();
		int neln = Nodes();

		m_node.resize( neln );
	}

	void UnpackTraitsData(int nflag) 
	{ 
		m_pT->m_pel = this;
		m_pT->UnpackData(nflag);
	}

	bool IsUnpacked() { return (m_pT->m_pel == this); }

	// interface to traits members
	vec3d* r0() { return m_pT->r0; }		// material coordintes
	vec3d* rt() { return m_pT->rt; }		// spatial coordinates
	vec3d* vt() { return m_pT->vt; }		// velocities

	double* pt() { return m_pT->pt;}		// nodal pressures

	int GaussPoints() const { return m_pT->nint; } 
	int Nodes() const { return m_pT->neln; } 

	double* H(int n) { return m_pT->H[n]; }		// shape function values

	int Type() const { return m_pT->m_ntype; } 

	vector<int>& LM() { return m_pT->LM; }

	//! Get the element's material ID
	int GetMatID() const { return m_mat; } 

	//! Set the element's material ID
	void SetMatID(int id) { m_mat = id; }

	//! Set the type of the element
	void SetType(int ntype)
	{
		FEElementLibrary::SetElementTraits(*this, ntype);
	}
/*
	bool HasNode(int n)
	{
		int m = m_node.size();
		for (int i=0; i<m; ++i) if (m[i] == n) return true;
		return false;
	}
*/
protected:
	// pointer to element traits
	FEElementTraits*	m_pT;

	int		m_mat;		//!< material index

public:

	int		m_nrigid;		//!< rigid body number that this element is attached to

	vector<int>		m_node;	//!< connectivity

	int	m_nID;				//!< element ID

	int	m_gid;	// part ID (i.e. index of domain this element belongs to)
};

//-----------------------------------------------------------------------------
//!  This class defines a solid element

class FESolidElement : public FEElement
{
public:
	//! default constructor
	FESolidElement(){}

	//! destructor
	virtual ~FESolidElement();

	//! copy constructor
	FESolidElement(const FESolidElement& el)
	{
		// set the traits of the element
		SetTraits(el.m_pT);

		// copy the state information
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i] = el.m_State[i]->Copy();

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy solid element data
		m_eJ = el.m_eJ;
		m_ep = el.m_ep;
		m_Lk = el.m_Lk;
	}

	//! assignment operator
	FESolidElement& operator = (const FESolidElement& el)
	{
		// make sure the element type is the same
		if (m_pT == 0) SetTraits(el.m_pT);
		else assert(m_pT == el.m_pT);

		// copy the state data
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			if (m_State[i]) 
			{
				delete m_State[i];
				m_State[i] = 0;
			}

			if (el.m_State[i]) m_State[i] = el.m_State[i]->Copy();
		}

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy solid element data
		m_eJ = el.m_eJ;
		m_ep = el.m_ep;
		m_Lk = el.m_Lk;

		return (*this);
	}



	virtual void SetTraits(FEElementTraits* ptraits)
	{
		FEElement::SetTraits(ptraits);

		int nint = GaussPoints();

		m_State.assign(nint, 0);
	}


	double* GaussWeights() { return ((FESolidElementTraits*)(m_pT))->gw; }			// weights of integration points

	double* Gr(int n) { return ((FESolidElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) { return ((FESolidElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s
	double* Gt(int n) { return ((FESolidElementTraits*)(m_pT))->Gt[n]; }	// shape function derivative to t

	void jact(double J[3][3], int n)
	{
		MATRIX3& Jt = ((FESolidElementTraits*)(m_pT))->m_Jt[n];
		J[0][0] = Jt[0][0]; J[0][1] = Jt[0][1]; J[0][2] = Jt[0][2];
		J[1][0] = Jt[1][0]; J[1][1] = Jt[1][1]; J[1][2] = Jt[1][2];
		J[2][0] = Jt[2][0]; J[2][1] = Jt[2][1]; J[2][2] = Jt[2][2];
	}

	void invjact(double J[3][3], int n)
	{
		MATRIX3& Jt = ((FESolidElementTraits*)(m_pT))->m_Jti[n];
		J[0][0] = Jt[0][0]; J[0][1] = Jt[0][1]; J[0][2] = Jt[0][2];
		J[1][0] = Jt[1][0]; J[1][1] = Jt[1][1]; J[1][2] = Jt[1][2];
		J[2][0] = Jt[2][0]; J[2][1] = Jt[2][1]; J[2][2] = Jt[2][2];
	}

	void jac0(double J[3][3], int n)
	{
		MATRIX3& J0 = ((FESolidElementTraits*)(m_pT))->m_J0[n];
		J[0][0] = J0[0][0]; J[0][1] = J0[0][1]; J[0][2] = J0[0][2];
		J[1][0] = J0[1][0]; J[1][1] = J0[1][1]; J[1][2] = J0[1][2];
		J[2][0] = J0[2][0]; J[2][1] = J0[2][1]; J[2][2] = J0[2][2];
	}

	void invjac0(double J[3][3], int n)
	{
		MATRIX3& J0 = ((FESolidElementTraits*)(m_pT))->m_J0i[n];
		J[0][0] = J0[0][0]; J[0][1] = J0[0][1]; J[0][2] = J0[0][2];
		J[1][0] = J0[1][0]; J[1][1] = J0[1][1]; J[1][2] = J0[1][2];
		J[2][0] = J0[2][0]; J[2][1] = J0[2][1]; J[2][2] = J0[2][2];
	}

	double detJt(int n) { return ((FESolidElementTraits*)(m_pT))->m_detJt[n]; }
	double detJ0(int n) { return ((FESolidElementTraits*)(m_pT))->m_detJ0[n]; }

	double defgrad(mat3d& F, int n)
	{
		// make sure this element is unpacked
		assert(m_pT->m_pel == this);

		double *Grn = Gr(n);
		double *Gsn = Gs(n);
		double *Gtn = Gt(n);

		vec3d* r = rt();

		double GX, GY, GZ;
		double Ji[3][3];
		invjac0(Ji, n);

		F[0][0] = F[0][1] = F[0][2] = 0;
		F[1][0] = F[1][1] = F[1][2] = 0;
		F[2][0] = F[2][1] = F[2][2] = 0;
		int neln = Nodes();
		for (int i=0; i<neln; ++i)
		{
			const double& Gri = Grn[i];
			const double& Gsi = Gsn[i];
			const double& Gti = Gtn[i];

			const double& x = r[i].x;
			const double& y = r[i].y;
			const double& z = r[i].z;

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			GX = Ji[0][0]*Gri+Ji[1][0]*Gsi+Ji[2][0]*Gti;
			GY = Ji[0][1]*Gri+Ji[1][1]*Gsi+Ji[2][1]*Gti;
			GZ = Ji[0][2]*Gri+Ji[1][2]*Gsi+Ji[2][2]*Gti;
	
			// calculate deformation gradient F
			F[0][0] += GX*x; F[0][1] += GY*x; F[0][2] += GZ*x;
			F[1][0] += GX*y; F[1][1] += GY*y; F[1][2] += GZ*y;
			F[2][0] += GX*z; F[2][1] += GY*z; F[2][2] += GZ*z;
		}

		double D = F.det();
		if (D <= 0) throw NegativeJacobian(m_nID, n, D, this);

		return D;
	}

	//! evaluate scalar field at integration point
	double Evaluate(double* fn, int n)
	{
		double* Hn = H(n);
		double f = 0;
		const int N = Nodes();
		for (int i=0; i<N; ++i) f += Hn[i]*fn[i];
		return f;
	}

	//! evaluate spatial gradient of scalar field at integration point
	vec3d gradient(double* fn, int n)
	{
		double Ji[3][3];
		invjact(Ji, n);
					
		double* Grn = Gr(n);
		double* Gsn = Gs(n);
		double* Gtn = Gt(n);

		double Gx, Gy, Gz;

		vec3d gradf;
		int N = Nodes();
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

	//! intialize element data
	void Init(bool bflag)
	{
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]->Init(bflag);
		m_eJ = 1;
		m_ep = 0;
		m_Lk = 0;
	}

	void SetMaterialPointData(FEMaterialPoint* pmp, int n) { m_State[n] = pmp; }

public:
	double	m_eJ;	//!< average dilatation
	double	m_ep;	//!< average pressure
	double	m_Lk;	//!< Lagrangian multiplier for incompressibility

	//! \todo the element state data will be removed from the element
	//! class in future versions of FEBio.
	vector<FEMaterialPoint*>	m_State;	//!< vector of state variables
};

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FESurfaceElement : public FEElement
{
public:
	FESurfaceElement() { m_nelem = -1; }

	FESurfaceElement(FESurfaceElement& el)
	{
		// set the traits of the element
		SetTraits(el.m_pT);

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy surface element data
		m_nelem = el.m_nelem;
		m_lnode = el.m_lnode;
	}

	FESurfaceElement& operator = (FESurfaceElement& el)
	{
		// make sure the element type is the same
		if (m_pT == 0) SetTraits(el.m_pT);
		else assert(m_pT == el.m_pT);

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy surface element data
		m_nelem = el.m_nelem;
		m_lnode = el.m_lnode;

		return (*this); 
	}

	virtual void SetTraits(FEElementTraits* pt)
	{
		FEElement::SetTraits(pt);
		m_lnode.resize(Nodes());
	}
	double* GaussWeights() { return ((FESurfaceElementTraits*)(m_pT))->gw; }			// weights of integration points

	double* Gr(int n) { return ((FESurfaceElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) { return ((FESurfaceElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

	double eval(double* d, int n)
	{
		double* N = H(n);
		int ne = Nodes();
		double a = 0;
		for (int i=0; i<ne; ++i) a += N[i]*d[i];
		return a;
	}

	double eval(double* d, double r, double s)
	{
		int n = Nodes();
		double H[4];
		shape_fnc(H, r, s);
		double a = 0;
		for (int i=0; i<n; ++i) a += H[i]*d[i];
		return a;
	}

	vec3d eval(vec3d* d, double r, double s)
	{
		int n = Nodes();
		double H[4];
		shape_fnc(H, r, s);
		vec3d a(0,0,0);
		for (int i=0; i<n; ++i) a += d[i]*H[i];
		return a;
	}

	vec3d eval(vec3d* d, int n)
	{
		int ne = Nodes();
		double* N = H(n);
		vec3d a(0,0,0);
		for (int i=0; i<ne; ++i) a += d[i]*N[i];
		return a;
	}

	double eval_deriv1(double* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hr[i]*d[i];
		return s;
	}

	double eval_deriv2(double* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hs[i]*d[i];
		return s;
	}

	double eval_deriv1(double* d, double r, double s)
	{
		double Hr[4], Hs[4];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hr[i]*d[i];
		return a;
	}

	double eval_deriv2(double* d, double r, double s)
	{
		double Hr[4], Hs[4];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hs[i]*d[i];
		return a;
	}

	void shape_fnc(double* H, double r, double s)
	{
		int n = Nodes();
		if (n == 4)
		{
			H[0] = 0.25*(1-r)*(1-s);
			H[1] = 0.25*(1+r)*(1-s);
			H[2] = 0.25*(1+r)*(1+s);
			H[3] = 0.25*(1-r)*(1+s);
		}
		else if (n==3)
		{
			H[0] = 1 - r - s;
			H[1] = r;
			H[2] = s;
		}
		else assert(false);
	}

	void shape_deriv(double* Gr, double* Gs, double r, double s)
	{
		int n = Nodes();
		if (n == 4)
		{
			Gr[0] = -0.25*(1-s); Gs[0] = -0.25*(1-r);
			Gr[1] =  0.25*(1-s); Gs[1] = -0.25*(1+r);
			Gr[2] =  0.25*(1+s); Gs[2] =  0.25*(1+r);
			Gr[3] = -0.25*(1+s); Gs[3] =  0.25*(1-r);
		}
		else if (n == 3)
		{
			Gr[0] = -1; Gs[0] = -1;
			Gr[1] =  1; Gs[1] =  0;
			Gr[2] =  0; Gs[2] =  1;
		}
		else assert(false);
	}

	//! this function projects data from the gauss-points to the nodal points
	void project_to_nodes(double* ai, double* ao)
	{
		int ni = GaussPoints();
		int ne = Nodes();
		assert(ni == ne); // TODO: for now we assume that the number of nodes equals the nr of gauss-points
		matrix& Hi = m_pT->Hi;
		for (int i=0; i<ne; ++i)
		{
			ao[i] = 0;
			for (int j=0; j<ni; ++j) ao[i] += Hi[i][j]*ai[j];
		}
	}

	bool HasNode(int n)
	{ 
		int l = Nodes(); 
		for (int i=0; i<l; ++i) 
			if (m_node[i] == n) return true; 
		return false;
	}

public:
	int		m_nelem;		//!< index of solid or shell element this surface element is a face of
	vector<int>	m_lnode;	//!< local node numbering (compared to m_node which is a global numbering)
};

//-----------------------------------------------------------------------------
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FEShellElement : public FEElement
{
public:
	FEShellElement(){}
	virtual ~FEShellElement();

	//! copy constructor
	FEShellElement(const FEShellElement& el)
	{
		// set the traits of the element
		SetTraits(el.m_pT);

		// we just copy the state information
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]= el.m_State[i]->Copy();

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;
	}

	//! assignment operator
	FEShellElement& operator = (const FEShellElement& el)
	{
		// make sure the element type is the same
		if (m_pT == 0) SetTraits(el.m_pT);
		else assert(m_pT == el.m_pT);

		// copy the state data
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) 
		{
			if (m_State[i]) delete m_State[i];
			m_State[i] = el.m_State[i]->Copy();
		}

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy shell data
		m_eJ = el.m_eJ;
		m_ep = el.m_ep;
		m_Lk = el.m_Lk;
		m_h0 = el.m_h0;

		return (*this);
	}

	virtual void SetTraits(FEElementTraits* ptraits)
	{
		FEElement::SetTraits(ptraits);
		m_h0.resize(Nodes());

		int nint = GaussPoints();

		m_State.assign(nint, 0);
	}

	double* GaussWeights() { return ((FEShellElementTraits*)(m_pT))->gw; }	// weights of integration points

	double* Hr(int n) { return ((FEShellElementTraits*)(m_pT))->Hr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FEShellElementTraits*)(m_pT))->Hs[n]; }	// shape function derivative to s

	vec3d* D0() { return ((FEShellElementTraits*)(m_pT))->D0; }
	vec3d* Dt() { return ((FEShellElementTraits*)(m_pT))->Dt; }

	void Init(bool bflag)
	{
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]->Init(bflag);
		m_eJ = 1;
		m_ep = 0;
		m_Lk = 0;
	}

	double gr(int n) { return ((FEShellElementTraits*)(m_pT))->gr[n]; }
	double gs(int n) { return ((FEShellElementTraits*)(m_pT))->gs[n]; }
	double gt(int n) { return ((FEShellElementTraits*)(m_pT))->gt[n]; }

	void invjac0(double J[3][3], int n)
	{
		MATRIX3& Jt = ((FEShellElementTraits*)(m_pT))->m_J0i[n];
		J[0][0] = Jt[0][0]; J[0][1] = Jt[0][1]; J[0][2] = Jt[0][2];
		J[1][0] = Jt[1][0]; J[1][1] = Jt[1][1]; J[1][2] = Jt[1][2];
		J[2][0] = Jt[2][0]; J[2][1] = Jt[2][1]; J[2][2] = Jt[2][2];
	}

	void invjact(double J[3][3], int n)
	{
		MATRIX3& Jt = ((FEShellElementTraits*)(m_pT))->m_Jti[n];
		J[0][0] = Jt[0][0]; J[0][1] = Jt[0][1]; J[0][2] = Jt[0][2];
		J[1][0] = Jt[1][0]; J[1][1] = Jt[1][1]; J[1][2] = Jt[1][2];
		J[2][0] = Jt[2][0]; J[2][1] = Jt[2][1]; J[2][2] = Jt[2][2];
	}

	double detJ0(int n) { return ((FEShellElementTraits*)(m_pT))->m_detJ0[n]; }
	double detJt(int n) { return ((FEShellElementTraits*)(m_pT))->m_detJt[n]; }

	double defgrad(mat3d& F, int n)
	{
		// make sure this element is unpacked
		assert(m_pT->m_pel == this);

		double* Hrn = Hr(n);
		double* Hsn = Hs(n);
		double* Hn  = H(n);
		double NX, NY, NZ, MX, MY, MZ;
		double za;

		vec3d* r = rt();
		vec3d* D = Dt();

		double g = gt(n);

		double Ji[3][3];
		invjac0(Ji, n);

		F[0][0] = F[0][1] = F[0][2] = 0;
		F[1][0] = F[1][1] = F[1][2] = 0;
		F[2][0] = F[2][1] = F[2][2] = 0;
		int neln = Nodes();
		for (int i=0; i<neln; ++i)
		{
			const double& Hri = Hrn[i];
			const double& Hsi = Hsn[i];
			const double& Hi  = Hn[i];

			const double& x = r[i].x;
			const double& y = r[i].y;
			const double& z = r[i].z;

			const double& dx = D[i].x;
			const double& dy = D[i].y;
			const double& dz = D[i].z;

			za = 0.5*g*m_h0[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			NX = Ji[0][0]*Hri+Ji[1][0]*Hsi;
			NY = Ji[0][1]*Hri+Ji[1][1]*Hsi;
			NZ = Ji[0][2]*Hri+Ji[1][2]*Hsi;

			MX = za*Ji[0][0]*Hri + za*Ji[1][0]*Hsi + Ji[2][0]*0.5*m_h0[i]*Hi;
			MY = za*Ji[0][1]*Hri + za*Ji[1][1]*Hsi + Ji[2][1]*0.5*m_h0[i]*Hi;
			MZ = za*Ji[0][2]*Hri + za*Ji[1][2]*Hsi + Ji[2][2]*0.5*m_h0[i]*Hi;

			// calculate deformation gradient F
			F[0][0] += NX*x + MX*dx; F[0][1] += NY*x + MY*dx; F[0][2] += NZ*x + MZ*dx;
			F[1][0] += NX*y + MX*dy; F[1][1] += NY*y + MY*dy; F[1][2] += NZ*y + MZ*dy;
			F[2][0] += NX*z + MX*dz; F[2][1] += NY*z + MY*dz; F[2][2] += NZ*z + MZ*dz;
		}

		double V = F.det();
		if (V <= 0) throw NegativeJacobian(m_nID, n, V, this);

		return V;
	}

	void SetMaterialPointData(FEMaterialPoint* pmp, int n) { m_State[n] = pmp; }

public:
	double	m_eJ;	//!< average dilatation
	double	m_ep;	//!< average pressure
	double	m_Lk;	//!< Lagrangian multiplier for incompressibility

	vector<double>	m_h0;	//!< initial shell thicknesses

	vector<FEMaterialPoint*>	m_State;
};

//-----------------------------------------------------------------------------

class FETrussElement : public FEElement
{
public:
	FETrussElement(){}
	virtual ~FETrussElement();

	FETrussElement(FETrussElement& el)
	{
		// set the traits of the element
		SetTraits(el.m_pT);

		// we just copy the state information
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) m_State[i]= el.m_State[i]->Copy();

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// truss data
		m_a0 = el.m_a0;
	}

	FETrussElement& operator = (FETrussElement& el) 
	{
		// make sure the element type is the same
		if (m_pT == 0) SetTraits(el.m_pT);
		else assert(m_pT == el.m_pT);

		// copy the state data
		int nint = GaussPoints();
		for (int i=0; i<nint; ++i) 
		{
			if (m_State[i]) delete m_State[i];
			m_State[i] = el.m_State[i]->Copy();
		}

		// copy base class data
		m_mat = el.m_mat;
		m_nrigid = el.m_nrigid;
		m_node = el.m_node;
		m_nID = el.m_nID;
		m_gid = el.m_gid;

		// copy truss data
		m_a0 = el.m_a0;

		return (*this); 
	}

	//! intialize element data
	void Init(bool bflag)
	{
		m_State[0]->Init(bflag);
	}

	double Length0()
	{
		assert(m_pT->m_pel == this);
		vec3d a0 = r0()[0];
		vec3d b0 = r0()[1];
		
		return (b0 - a0).norm();
	}

	double Length()
	{
		assert(m_pT->m_pel == this);
		vec3d a = rt()[0];
		vec3d b = rt()[1];
		
		return (b - a).norm();
	}

	vec3d Normal()
	{
		assert(m_pT->m_pel == this);
		assert(m_pT->m_pel == this);
		vec3d a = rt()[0];
		vec3d b = rt()[1];
		vec3d n = b - a;
		n.unit();
		return n;
	}

	double Volume0()
	{
		double L = Length0();
		return m_a0*L;
	}

	virtual void SetTraits(FEElementTraits* ptraits)
	{
		FEElement::SetTraits(ptraits);
		m_State.assign(1, 0);
	}

	void SetMaterialPointData(FEMaterialPoint* pmp, int n) { m_State[n] = pmp; }

public:
	double	m_a0;	// cross-sectional area

	vector<FEMaterialPoint*>	m_State;
};

#endif // !defined(AFX_FEELEMENT_H__2EE38101_58E2_4FEB_B214_BB71B6FB15FB__INCLUDED_)
