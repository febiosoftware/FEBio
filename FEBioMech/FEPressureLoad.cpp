#include "stdafx.h"
#include "FEPressureLoad.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEPressureLoad, FESurfaceTraction)
	ADD_PARAMETER(m_pressure, "pressure");
	ADD_PARAMETER(m_bsymm   , "symmetric_stiffness");
	ADD_PARAMETER(m_blinear, "linear");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceTraction(pfem)
{ 
	m_pressure = 0.0;
	m_bsymm = true;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_pressure.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::ElementStiffness(FESurfaceElement& el, matrix& ke)
{
	// choose the symmetric of unsymmetric formulation
	if (m_bsymm) 
		SymmetricPressureStiffness(el, ke);
	else
		UnsymmetricPressureStiffness(el, ke);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	GetNodalCoordinates(el, rt);

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i)
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure(pt);

        if (m_bshellb) P = -P;
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
					   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*P;

				ke.add(3*i, 3*j, mat3da(kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureLoad::UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	GetNodalCoordinates(el, rt);

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure(pt);

        if (m_bshellb) P = -P;
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d Kab = (dxr*Gs[j] - dxs*Gr[j])*(P*N[i]*w[n]);
				ke.sub(3*i, 3*j, mat3da(Kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

vec3d FEPressureLoad::Traction(const FESurfaceMaterialPoint& pt)
{
	// evaluate pressure at this material point
	double P = -m_pressure(pt);

	// force vector
	vec3d N = (pt.dr ^ pt.ds); N.unit();
	vec3d f = N*P;

	return f;
}

//-----------------------------------------------------------------------------

void FEPressureLoad::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
}
