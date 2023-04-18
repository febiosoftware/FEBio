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
#include "FEVolumeConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEDataExport.h>
#include <FECore/DumpStream.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
FEVolumeSurface::FEVolumeSurface(FEModel* fem) : FESurface(fem)
{
	m_Lp = 0.0;
	m_p  = 0.0;
	m_V0 = 0.0;
	m_Vt = 0.0;

	// define class exports
	EXPORT_DATA(PLT_FLOAT, FMT_REGION, &m_p, "volume pressure");
}

//-----------------------------------------------------------------------------
bool FEVolumeSurface::Init()
{
	if (FESurface::Init() == false) return false;

	// evaluate the initial volume
	m_V0 = Volume();
	m_Vt = m_V0;

	return (m_V0 != 0.0);
}

//-----------------------------------------------------------------------------
void FEVolumeSurface::CopyFrom(FEVolumeSurface& s)
{
	m_Node = s.m_Node;

	// create elements
	int NE = s.Elements();
	Create(NE);
	for (int i=0; i<NE; ++i) Element(i) = s.Element(i);

	// copy surface data
	m_Lp = s.m_Lp;
	m_p  = s.m_p;
	m_V0 = s.m_V0;
	m_Vt = s.m_Vt;
}

//-----------------------------------------------------------------------------
void FEVolumeSurface::Serialize(DumpStream& ar)
{
	FESurface::Serialize(ar);
	ar & m_Lp & m_p & m_V0 & m_Vt;
}

//-----------------------------------------------------------------------------
//! Calculate the initial volume
double FEVolumeSurface::Volume()
{
	// get the mesh
	FEMesh& mesh = *GetMesh();

	// loop over all elements
	double vol = 0.0;
	int NE = Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& el = Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

		// loop over integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// evaluate the position vector at this point
			vec3d r = el.eval(x, n);

			// calculate the tangent vectors
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int j=0; j<neln; ++j) 
			{
				dxr += x[j]*Gr[j];
				dxs += x[j]*Gs[j];
			}

			// update volume
			vol += w[n]*(r*(dxr^dxs));
		}
	}
	return vol/3.0;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEVolumeConstraint, FESurfaceConstraint);
	ADD_PARAMETER(m_blaugon, "laugon" ); 
	ADD_PARAMETER(m_atol   , "augtol" );
	ADD_PARAMETER(m_eps    , "penalty");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor. Set default parameter values
FEVolumeConstraint::FEVolumeConstraint(FEModel* pfem) : FESurfaceConstraint(pfem)
{
	m_s = new FEVolumeSurface(pfem);

	m_eps = 0.0;
	m_atol = 0.0;
	m_blaugon = false;
	m_binit = false;	// will be set to true during activation

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
FEVolumeConstraint::~FEVolumeConstraint()
{
	delete m_s;
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::CopyFrom(FENLConstraint* plc)
{
	// cast to a periodic boundary
	FEVolumeConstraint& vc = dynamic_cast<FEVolumeConstraint&>(*plc);

	// copy parameters
	GetParameterList() = vc.GetParameterList();

	// copy nodes
	m_s->CopyFrom(*vc.m_s);
}

//-----------------------------------------------------------------------------
//! Returns the surface
FESurface* FEVolumeConstraint::GetSurface()
{
	return m_s;
}

//-----------------------------------------------------------------------------
double FEVolumeConstraint::EnclosedVolume() const
{
	return m_s->m_Vt;
}

//-----------------------------------------------------------------------------
double FEVolumeConstraint::Pressure() const
{
	return m_s->m_p;
}

//-----------------------------------------------------------------------------
//! Initializes data structures. 
void FEVolumeConstraint::Activate()
{
	// don't forget to call base class
	FENLConstraint::Activate();

	// initialize the surface
	if (m_binit == false) m_s->Init();

	// set flag that initial volume is calculated
	m_binit = true;
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// We don't do anything here since the connectivity of a surface
	// is implied by the domain to which it is attached.
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEVolumeSurface& s = *m_s;

	FEMesh& mesh = *s.GetMesh();

	vector<double> fe;
	vector<int> lm;

	// get the pressure
	double p = s.m_p;

	// loop over all elements
	int NE = s.Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& el = s.Element(i);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

		// allocate element residual vector
		int ndof = 3*neln;
		fe.resize(ndof);
		zero(fe);

		// loop over all integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// calculate the tangent vectors
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int j=0; j<neln; ++j) 
			{
				dxr += x[j]*Gr[j];
				dxs += x[j]*Gs[j];
			}

			// evaluate the "normal" vector
			vec3d v = (dxr ^ dxs)*w[n]*p;

			// evaluate the element forces
			double* H = el.H(n);
			for (int j=0; j<neln; ++j)
			{
				fe[3*j  ] += H[j]*v.x;
				fe[3*j+1] += H[j]*v.y;
				fe[3*j+2] += H[j]*v.z;
			}
		}

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEVolumeSurface& s = *m_s;

	FEMesh& mesh = *s.GetMesh();

	// get the pressure
	double p = s.m_p;

	// element stiffness matrix
	vector<int> lm;
	vector<double> fe;

	// loop over all elements
	int NE = s.Elements();
	vec3d x[FEElement::MAX_NODES];
	for (int l=0; l<NE; ++l)
	{
		// get the next element
		FESurfaceElement& el = s.Element(l);

		FEElementMatrix ke(el);

		// get the nodal coordinates
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

		// allocate the stiffness matrix
		int ndof = 3*neln;
		ke.resize(ndof, ndof);
		ke.zero();
		fe.resize(ndof);
		zero(fe);

		// repeat over integration points
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n=0; n<nint; ++n)
		{
			// calculate tangent vectors
			double* N = el.H(n);
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int j=0; j<neln; ++j) 
			{
				dxr += x[j]*Gr[j];
				dxs += x[j]*Gs[j];
			}

			// calculate pressure contribution
			vec3d v = (dxr ^ dxs)*w[n]*m_eps;
			for (int i=0; i<neln; ++i)
			{
				fe[3*i  ] += N[i]*v.x;
				fe[3*i+1] += N[i]*v.y;
				fe[3*i+2] += N[i]*v.z;
			}
			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					ke[3*i  ][3*j  ] += fe[3*i  ]*fe[3*j  ];
					ke[3*i  ][3*j+1] += fe[3*i  ]*fe[3*j+1];
					ke[3*i  ][3*j+2] += fe[3*i  ]*fe[3*j+2];

					ke[3*i+1][3*j  ] += fe[3*i+1]*fe[3*j  ];
					ke[3*i+1][3*j+1] += fe[3*i+1]*fe[3*j+1];
					ke[3*i+1][3*j+2] += fe[3*i+1]*fe[3*j+2];

					ke[3*i+2][3*j  ] += fe[3*i+2]*fe[3*j  ];
					ke[3*i+2][3*j+1] += fe[3*i+2]*fe[3*j+1];
					ke[3*i+2][3*j+2] += fe[3*i+2]*fe[3*j+2];
				}

		
			// calculate displacement contribution
			vec3d kab;
			for (int i=0; i<neln; ++i)
				for (int j=0; j<neln; ++j)
				{
					kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
						   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*p;

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


		// get the element's LM vector
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
bool FEVolumeConstraint::Augment(int naug, const FETimeInfo& tp)
{
	FEVolumeSurface& s = *m_s;

	// make sure we are augmenting
	if ((m_blaugon == false) || (m_atol <= 0.0)) return true;

	feLog("\nvolume constraint:\n");

	double Dp = m_eps*(s.m_Vt - s.m_V0);
	double Lp = s.m_p;
	double err = fabs(Dp/Lp);
	feLog("\tpressure: %lg\n", Lp);
	feLog("\tnorm : %lg (%lg)\n", err, m_atol);
	feLog("\tvolume ratio: %lg\n", s.m_Vt / s.m_V0);

	// check convergence
	if (err < m_atol) return true;

	// update Lagrange multiplier (and pressure variable)
	s.m_Lp = Lp;
	s.m_p = Lp + Dp;

	return false;
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	m_s->Serialize(ar);
	ar & m_binit;
}

//-----------------------------------------------------------------------------
void FEVolumeConstraint::Reset()
{
}

//-----------------------------------------------------------------------------
// This function is called when the FE model's state needs to be updated.
void FEVolumeConstraint::Update()
{
	FEVolumeSurface& s = *m_s;

	// calculate the current volume
	s.m_Vt = s.Volume();

	// update pressure variable
	s.m_p = s.m_Lp + m_eps*(s.m_Vt - s.m_V0);
}
