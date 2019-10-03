/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FESolutesDomain.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioFluidSolutes.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FESolutesDomain::FESolutesDomain(FEModel* pfem) : FESolidDomain(pfem), m_dof(pfem)
{
	m_pMat = 0;
	m_btrans = true;

	m_dofC = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
	m_dofAC = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION_TDERIV), 0);

	// list the degrees of freedom
	// (This allows the FEDomain base class to handle several tasks such as UnpackLM)
//	vector<int> dof;
//	SetDOFList(dof);
}

//-----------------------------------------------------------------------------
//! get the total dofs
const FEDofList& FESolutesDomain::GetDOFList() const
{
	return m_dof;
}

//-----------------------------------------------------------------------------
//! Assign material
void FESolutesDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	if (pmat)
	{
		m_pMat = dynamic_cast<FESolutesMaterial*>(pmat);
		assert(m_pMat);
	}
	else m_pMat = 0;
}

//-----------------------------------------------------------------------------
bool FESolutesDomain::Init()
{
	// initialize base class
	if (FESolidDomain::Init() == false) return false;

	const int nsol = m_pMat->Solutes();

	// set the active degrees of freedom list
	m_dof.Clear();
	for (int i = 0; i<nsol; ++i)
	{
		int m = m_pMat->GetSolute(i)->GetSoluteDOF();
		m_dof.AddDof(m_dofC + m);
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESolutesDomain::Reset()
{
	// reset base class
	FESolidDomain::Reset();

	const int nsol = m_pMat->Solutes();

	for (int i = 0; i<(int)m_Elem.size(); ++i)
	{
		// get the solid element
		FESolidElement& el = m_Elem[i];

		// get the number of integration points
		int nint = el.GaussPoints();

		// loop over the integration points
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FESolutesMaterial::Point& ps = *(mp.ExtractData<FESolutesMaterial::Point>());

			// initialize solutes
			ps.m_nsol = nsol;
			ps.m_c.assign(nsol, 0);
			ps.m_cdot.assign(nsol, 0);
			ps.m_gradc.assign(nsol, vec3d(0, 0, 0));
			ps.m_j.assign(nsol, vec3d(0, 0, 0));
		}
	}
}

//-----------------------------------------------------------------------------
void FESolutesDomain::Activate()
{
	const int nsol = m_pMat->Solutes();
	for (int i = 0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			for (int isol = 0; isol<nsol; ++isol)
				node.set_active(m_dofC + m_pMat->GetSolute(isol)->GetSoluteDOF());
		}
	}
}

//-----------------------------------------------------------------------------
void FESolutesDomain::InitMaterialPoints()
{
	const int nsol = m_pMat->Solutes();
	FEMesh& m = *GetMesh();

	const int NE = FEElement::MAX_NODES;
	vector< vector<double> > c0(nsol, vector<double>(NE));
	vector<int> sid(nsol);
	for (int j = 0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();

	for (int j = 0; j<(int)m_Elem.size(); ++j)
	{
		// get the solid element
		FESolidElement& el = m_Elem[j];

		// get the number of nodes
		int neln = el.Nodes();
		// get initial values of fluid pressure and solute concentrations
		for (int i = 0; i<neln; ++i)
		{
			FENode& ni = m.Node(el.m_node[i]);
			for (int isol = 0; isol<nsol; ++isol)
				c0[isol][i] = ni.get(m_dofC + sid[isol]);
		}

		// get the number of integration points
		int nint = el.GaussPoints();

		// loop over the integration points
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FESolutesMaterial::Point& ps = *(mp.ExtractData<FESolutesMaterial::Point>());

			// initialize solutes
			ps.m_nsol = nsol;

			// initialize effective solute concentrations
			for (int isol = 0; isol<nsol; ++isol) {
				ps.m_c[isol] = el.Evaluate(c0[isol], n);
				ps.m_gradc[isol] = gradient(el, c0[isol], n);
			}

			for (int isol = 0; isol<nsol; ++isol)
				ps.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FESolutesDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], r0, v;
	FEMesh& m = *GetMesh();
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int i = 0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
		}

		int n = el.GaussPoints();
		for (int j = 0; j<n; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mp.m_r0 = el.Evaluate(x0, j);
			mp.Update(timeInfo);
		}
	}
}

//-----------------------------------------------------------------------------
void FESolutesDomain::InternalForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	int NE = (int)m_Elem.size();
#pragma omp parallel for shared (NE)
	for (int i = 0; i<NE; ++i)
	{
		// element force vector
		vector<double> fe;
		vector<int> lm;

		// get the element
		FESolidElement& el = m_Elem[i];

		// get the element force vector and initialize it to zero
		int nsol = m_pMat->Solutes();
		int ndof = nsol*el.Nodes();
		fe.assign(ndof, 0);

		// calculate internal force vector
		ElementInternalForce(el, fe, tp);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element 'fe'-vector into global R vector
		//#pragma omp critical
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
//! calculates the internal equivalent nodal forces for solid elements

void FESolutesDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe, const FETimeInfo& tp)
{
	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3];

	// get the time step size
	double dt = tp.timeIncrement;

	// number of solutes
	const int nsol = m_pMat->Solutes();

	// gradient of shape functions
	int neln = el.Nodes();
	vector<vec3d> gradN(neln);

	// repeat for all integration points
	int nint = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n = 0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FESolutesMaterial::Point& spt = *(mp.ExtractData<FESolutesMaterial::Point>());

		// calculate the jacobian
		double detJ = invjac0(el, Ji, n)*gw[n];

		vec3d g1(Ji[0][0], Ji[0][1], Ji[0][2]);
		vec3d g2(Ji[1][0], Ji[1][1], Ji[1][2]);
		vec3d g3(Ji[2][0], Ji[2][1], Ji[2][2]);

		double* H = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);
		double* Gt = el.Gt(n);

		// evaluate spatial gradient of shape functions
		for (int i = 0; i<neln; ++i)
		{
			gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];
		}

		for (int i = 0; i<neln; ++i)
		{
			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			for (int isol = 0; isol<nsol; ++isol)
				fe[nsol*i + isol] -= (spt.m_j[isol] * gradN[i] - H[i] * (spt.m_cdot[isol] + spt.m_c[isol] * spt.m_divf))*detJ*dt;
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FESolutesDomain::ElementStiffness(FESolidElement &el, matrix &ke, const FETimeInfo& tp)
{
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int nsol = m_pMat->Solutes();

	// gradient of shape functions
	vector<vec3d> gradN(neln);

	double dt = tp.timeIncrement;
	double ksi = tp.alpham / (tp.gamma*tp.alphaf);

	// jacobian
	double Ji[3][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n = 0; n<nint; ++n)
	{
		// calculate jacobian
		double detJ = invjac0(el, Ji, n)*gw[n];

		vec3d g1(Ji[0][0], Ji[0][1], Ji[0][2]);
		vec3d g2(Ji[1][0], Ji[1][1], Ji[1][2]);
		vec3d g3(Ji[2][0], Ji[2][1], Ji[2][2]);

		double* H = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);
		double* Gt = el.Gt(n);

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FESolutesMaterial::Point& spt = *(mp.ExtractData<FESolutesMaterial::Point>());

		// evaluate spatial gradient of shape functions
		for (int i = 0; i<neln; ++i)
			gradN[i] = g1*Gr[i] + g2*Gs[i] + g3*Gt[i];

		// evaluate stiffness matrix
		for (int i = 0; i<neln; ++i)
		{
			for (int j = 0; j<neln; ++j)
			{
				for (int isol = 0; isol<nsol; ++isol) {
					double d0 = m_pMat->GetSolute(isol)->m_pDiff->Free_Diffusivity(mp);
					double d0p = m_pMat->GetSolute(isol)->m_pDiff->Tangent_Free_Diffusivity_Concentration(mp, isol);
					double kcc = -(H[i] * ((ksi / dt + spt.m_divf)*H[j] + gradN[j] * spt.m_vft) + (gradN[j] * d0 + spt.m_gradc[isol] * H[j] * d0p)*gradN[i]);
					ke[i*nsol + isol][j*nsol + isol] += kcc*detJ*dt;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESolutesDomain::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// repeat over all solid elements
	int NE = (int)m_Elem.size();

#pragma omp parallel for shared (NE)
	for (int iel = 0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);

		// create the element's stiffness matrix
		int nsol = m_pMat->Solutes();
		int ndof = nsol*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate material stiffness
		ElementStiffness(el, ke, tp);

		// get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
#pragma omp critical
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FESolutesDomain::Update(const FETimeInfo& tp)
{
	bool berr = false;
	int NE = (int)m_Elem.size();
#pragma omp parallel for shared(NE, berr)
	for (int i = 0; i<NE; ++i)
	{
		try
		{
			UpdateElementStress(i, tp);
		}
		catch (NegativeJacobian e)
		{
#pragma omp critical
			{
				// reset the logfile mode
				berr = true;
				if (NegativeJacobian::DoOutput()) feLogError(e.what());
			}
		}
	}

	// if we encountered an error, we request a running restart
	if (berr)
	{
		if (NegativeJacobian::DoOutput() == false) feLogError("Negative jacobian was detected.");
		throw DoRunningRestart();
	}
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
void FESolutesDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
	// time constants
	double alphaf = tp.alphaf;
	double alpham = tp.alpham;

	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// number of nodes
	int neln = el.Nodes();

	// number of solutes
	const int nsol = m_pMat->Solutes();

	// nodal coordinates
	const int NELN = FEElement::MAX_NODES;
	vector< vector<double> > ct(nsol, vector<double>(NELN));
	vector< vector<double> > cp(nsol, vector<double>(NELN));
	vector< vector<double> > act(nsol, vector<double>(NELN));
	vector< vector<double> > acp(nsol, vector<double>(NELN));
	vector<int> sid(nsol);
	for (int j = 0; j<nsol; ++j) sid[j] = m_pMat->GetSolute(j)->GetSoluteDOF();
	for (int j = 0; j<neln; ++j) {
		FENode& node = m_pMesh->Node(el.m_node[j]);
		for (int k = 0; k<nsol; ++k) {
			ct[k][j] = node.get(m_dofC + sid[k]);
			cp[k][j] = node.get_prev(m_dofC + sid[k]);
			act[k][j] = node.get(m_dofAC + sid[k]);
			acp[k][j] = node.get_prev(m_dofAC + sid[k]);
		}
	}

	// loop over the integration points and update
	// velocity, velocity gradient, acceleration
	// stress and pressure at the integration point
	int nint = el.GaussPoints();
	for (int n = 0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FESolutesMaterial::Point& spt = *(mp.ExtractData<FESolutesMaterial::Point>());

		// material point data
		for (int isol = 0; isol < nsol; ++isol) {
			spt.m_c[isol] = el.Evaluate(ct[isol], n)*alphaf + el.Evaluate(cp[isol], n)*(1 - alphaf);
			spt.m_gradc[isol] = gradient(el, ct[isol], n)*alphaf + gradient(el, cp[isol], n)*(1 - alphaf);
			spt.m_cdot[isol] = spt.m_gradc[isol] * spt.m_vft;
			if (m_btrans) spt.m_cdot[isol] += el.Evaluate(act[isol], n)*alpham + el.Evaluate(acp[isol], n)*(1 - alpham);
		}

		// calculate the solute flux
		for (int isol = 0; isol < nsol; ++isol)
			spt.m_j[isol] = m_pMat->SoluteFlux(mp, isol);
	}
}
