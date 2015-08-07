#include "stdafx.h"
#include "FEPressureLoad.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEPressureLoad::LOAD::LOAD()
{ 
	lc = -1;
	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = 1.0;
}

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAMETER_LIST(FEPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_blinear , FE_PARAM_BOOL  , "linear"  );
	ADD_PARAMETER(m_pressure, FE_PARAM_DOUBLE, "pressure");
	ADD_PARAMETER(m_bsymm   , FE_PARAM_BOOL  , "symmetric_stiffness");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_blinear = false;
	m_pressure = 1.0;
	m_bsymm = true;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureLoad::Create(int n)
{
	m_PC.resize(n); 
}

//-----------------------------------------------------------------------------
//! \deprecated This function is only used by the 1.2 file reader and is to be 
//! considered obsolete.
bool FEPressureLoad::SetAttribute(const char* szatt, const char* szval)
{
	if (strcmp(szatt, "type") == 0)
	{
		if      (strcmp(szval, "linear"   ) == 0) SetLinear(true );
		else if (strcmp(szval, "nonlinear") == 0) SetLinear(false);
		else return false;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPressureLoad::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& pc = PressureLoad(nface);
	if      (strcmp(szatt, "id") == 0) {}
	else if (strcmp(szatt, "lc") == 0) pc.lc = atoi(szval) - 1;
	else if (strcmp(szatt, "scale") == 0)
	{
		double s = atof(szval);
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;
		pc.s[8] = s;
	}
	else return false;

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	// choose the symmetric of unsymmetric formulation
	if (m_bsymm) 
		SymmetricPressureStiffness(el, ke, tn);
	else
		UnsymmetricPressureStiffness(el, ke, tn);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// traction at integration point
	double tr;
	
	vec3d dxr, dxs;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

	vec3d kab;

	ke.zero();

	double* N, *Gr, *Gs;

	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
					   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*tr;

				ke[3*i  ][3*j  ] +=      0;
				ke[3*i  ][3*j+1] += -kab.z;
				ke[3*i  ][3*j+2] +=  kab.y;

				ke[3*i+1][3*j  ] +=  kab.z;
				ke[3*i+1][3*j+1] +=      0;
				ke[3*i+1][3*j+2] += -kab.x;

				ke[3*i+2][3*j  ] += -kab.y;
				ke[3*i+2][3*j+1] +=  kab.x;
				ke[3*i+2][3*j+2] +=      0;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureLoad::UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// traction at integration point
	double tr;
	
	vec3d dxr, dxs;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

	vec3d kab;

	ke.zero();

	double* N, *Gr, *Gs;

	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
        mat3d rw; rw.skew(dxr);
        mat3d sw; sw.skew(dxs);
        
		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
                mat3d Kab = (rw*Gs[j] - sw*Gr[j])*(tr*N[i]*w[n]);
                ke[3*i  ][3*j  ] -=  Kab(0,0);
                ke[3*i  ][3*j+1] -=  Kab(0,1);
                ke[3*i  ][3*j+2] -=  Kab(0,2);
                
                ke[3*i+1][3*j  ] -=  Kab(1,0);
                ke[3*i+1][3*j+1] -=  Kab(1,1);
                ke[3*i+1][3*j+2] -=  Kab(1,2);
                
                ke[3*i+2][3*j  ] -=  Kab(2,0);
                ke[3*i+2][3*j+1] -=  Kab(2,1);
                ke[3*i+2][3*j+2] -=  Kab(2,2);
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPressureLoad::PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// traction at integration points
	double tr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		f = (dxr ^ dxs)*tr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPressureLoad::LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// traction at integration points
	double tr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += r0[i]*Gr[i];
			dxs += r0[i]*Gs[i];
		}

		f = (dxr ^ dxs)*tr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void FEPressureLoad::Serialize(DumpFile& ar)
{
	FESurfaceLoad::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << (int) m_PC.size();
		for (int i=0; i< (int) m_PC.size(); ++i)
		{
			LOAD& pc = m_PC[i];
			ar << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
			ar << pc.s[4] << pc.s[5] << pc.s[6] << pc.s[7] << pc.s[8];
		}
	}
	else
	{
		int n;
		ar >> n;
		m_PC.resize(n);
		// pressure forces
		for (int i=0; i<n; ++i)
		{
			LOAD& pc = m_PC[i];
			ar >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
			ar >> pc.s[4] >> pc.s[5] >> pc.s[6] >> pc.s[7] >> pc.s[8];
		}
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::StiffnessMatrix(FESolver* psolver)
{
	FEModel& fem = psolver->GetFEModel();

	matrix ke;
	vector<int> lm;

	int npr = m_PC.size();
	for (int m=0; m<npr; ++m)
	{
		LOAD& pc = m_PC[m];
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// calculate nodal normal tractions
		int neln = el.Nodes();
		vector<double> tn(neln);

		if (m_blinear == false)
		{
			double g = m_pressure;
			if (pc.lc >= 0) g *= fem.GetLoadCurve(pc.lc)->Value();

			// evaluate the prescribed traction.
			// note the negative sign. This is because this boundary condition uses the 
			// convention that a positive pressure is compressive
			for (int j=0; j<neln; ++j) tn[j] = -g*pc.s[j];

			// get the element stiffness matrix
			int ndof = 3*neln;
			ke.resize(ndof, ndof);

			// calculate pressure stiffness
			PressureStiffness(el, ke, tn);

			// get the element's LM vector
			m_psurf->UnpackLM(el, lm);

			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::Residual(FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();

	vector<double> fe;
	vector<int> lm;

	int npr = m_PC.size();
	for (int i=0; i<npr; ++i)
	{
		LOAD& pc = m_PC[i];
		FESurfaceElement& el = m_psurf->Element(i);

		// calculate nodal normal tractions
		int neln = el.Nodes();
		vector<double> tn(neln);

		double g = m_pressure;
		if (pc.lc >= 0) g *= fem.GetLoadCurve(pc.lc)->Value();

		// evaluate the prescribed traction.
		// note the negative sign. This is because this boundary condition uses the 
		// convention that a positive pressure is compressive
		for (int j=0; j<el.Nodes(); ++j) tn[j] = -g*pc.s[j];
		
		int ndof = 3*neln;
		fe.resize(ndof);

		if (m_blinear) LinearPressureForce(el, fe, tn); else PressureForce(el, fe, tn);

		// get the element's LM vector
		m_psurf->UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
