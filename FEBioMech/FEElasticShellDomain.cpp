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
#include "FEElasticShellDomain.h"
#include "FEElasticMaterial.h"
#include "FEBodyForce.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <math.h>
#include <FECore/FESolidDomain.h>
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
FEElasticShellDomain::FEElasticShellDomain(FEModel* pfem) : FESSIShellDomain(pfem), FEElasticDomain(pfem), m_dofV(pfem), m_dofSV(pfem), m_dofSA(pfem), m_dofR(pfem)
{
	m_pMat = 0;
    m_alphaf = m_beta = 1;
    m_alpham = 2;

	m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCTIY));
	m_dofSV.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY));
	m_dofSA.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION));
	m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
}

//-----------------------------------------------------------------------------
FEElasticShellDomain& FEElasticShellDomain::operator = (FEElasticShellDomain& d)
{ 
	m_Elem = d.m_Elem; 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::SetMaterial(FEMaterial* pmat)
{
	FEDomain::SetMaterial(pmat);
	m_pMat = dynamic_cast<FESolidMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::Activate()
{
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(m_dofU[0]);
				node.set_active(m_dofU[1]);
				node.set_active(m_dofU[2]);

				if (node.HasFlags(FENode::SHELL))
				{
					node.set_active(m_dofSU[0]);
					node.set_active(m_dofSU[1]);
					node.set_active(m_dofSU[2]);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEElasticShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    m_alphaf = timeInfo.alphaf;
    m_alpham = timeInfo.alpham;
    m_beta = timeInfo.beta;
    
    vec3d r0, rt;
    for (size_t i=0; i<m_Elem.size(); ++i)
    {
        FEShellElement& el = m_Elem[i];
        int n = el.GaussPoints();
        for (int j=0; j<n; ++j)
        {
            FEMaterialPoint& mp = *el.GetMaterialPoint(j);
            FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
            pt.m_Wp = pt.m_Wt;
            
            mp.Update(timeInfo);
        }
    }
}

//-----------------------------------------------------------------------------
// Calculates the forces due to the stress
void FEElasticShellDomain::InternalForces(FEGlobalVector& R)
{
    int NS = (int)m_Elem.size();
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
        R.Assemble(el.m_node, lm, fe, true);
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
        detJt = detJ(el, n, m_alphaf)*gw[n];

		// get the stress vector for this integration point
		mat3ds& s = pt.m_s;

		eta = el.gt(n);

		Mr = el.Hr(n);
		Ms = el.Hs(n);
		M  = el.H(n);
        
        ContraBaseVectors(el, n, gcnt, m_alphaf);

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
    int NS = (int)m_Elem.size();
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
        R.Assemble(el.m_node, lm, fe, true);
    }
}

//-----------------------------------------------------------------------------
//! Calculates element body forces for shells

void FEElasticShellDomain::ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe)
{
    FEParamDouble& density = m_pMat->Density();
    
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
		double dens = density(mp);
        
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
    int NE = (int)m_Elem.size();
    for (int i=0; i<NE; ++i)
    {
        // element force vector
        vector<double> fe;
        vector<int> lm;
        
        // get the element
        FEShellElement& el = m_Elem[i];
        
        // get the element force vector and initialize it to zero
        int ndof = 6*el.Nodes();
        fe.assign(ndof, 0);
        
        // calculate internal force vector
        ElementInertialForce(el, fe);
        
        // get the element's LM vector
        UnpackLM(el, lm);
        
        // assemble element 'fe'-vector into global R vector
        R.Assemble(el.m_node, lm, fe, true);
    }
}

//-----------------------------------------------------------------------------
void FEElasticShellDomain::ElementInertialForce(FEShellElement& el, vector<double>& fe)
{
    int nint = el.GaussPoints();
    int neln = el.Nodes();

	FEParamDouble& density = m_pMat->Density();
    
    // evaluate the element inertial force vector
    for (int n=0; n<nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        double dens = density(mp);
        double J0 = detJ0(el, n)*el.GaussWeights()[n];
        
        double* M = el.H(n);
        double eta = el.gt(n);
        
        for (int i=0; i<neln; ++i)
        {
            vec3d fu = pt.m_a*(dens*M[i]*(1+eta)/2*J0);
            vec3d fd = pt.m_a*(dens*M[i]*(1-eta)/2*J0);
            
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
//! This function calculates the stiffness due to body forces
void FEElasticShellDomain::ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke)
{
    int i, j, i6, j6;
    int neln = el.Nodes();
    
    // don't forget to multiply with the density
    FEParamDouble& density = m_pMat->Density();
    
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
        detJ = detJ0(el, n)*gw[n]*m_alphaf;
        
        // get the stiffness
        K = BF.stiffness(mp)*density(mp)*detJ;
        
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
    int NS = (int)m_Elem.size();
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
        detJt = detJ(el, n, m_alphaf)*gw[n]*m_alphaf;
        
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
    FEParamDouble& density = m_pMat->Density();
    
    // calculate element stiffness matrix
    for (int n=0; n<nint; ++n)
    {
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		double D = density(mp);

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
            FEParamDouble& density = m_pMat->Density();
            
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

				double dens0 = density(mp);
                
                // calculate the jacobian
                detJt = detJ0(el, n)*gw[n];
                
                M  = el.H(n);
                eta = el.gt(n);
                
                // get the force
                vec3d f = pbf->force(mp);
                
                for (int i=0; i<neln; ++i)
                {
                    vec3d fu = f*(dens0*M[i]*(1+eta)/2);
                    vec3d fd = f*(dens0*M[i]*(1-eta)/2);
                    
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
	FESSIShellDomain::Update(tp);

    double dt = tp.timeIncrement;
    
    const int NELN = FEElement::MAX_NODES;
	vec3d r0[NELN], s0[NELN], r[NELN], s[NELN];
    vec3d v[NELN], w[NELN];
    vec3d a[NELN], b[NELN];

	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(m_pMat);

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
            FENode& node = m_pMesh->Node(el.m_node[j]);
            r0[j] = node.m_r0;
            s0[j] = node.m_r0 - node.m_d0;
            r[j] = node.m_rt*m_alphaf + node.m_rp*(1-m_alphaf);
            s[j] = s0[j] + node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2])*m_alphaf + node.get_vec3d_prev(m_dofSU[0], m_dofSU[1], m_dofSU[2])*(1-m_alphaf);
            v[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*m_alphaf + node.m_vp*(1-m_alphaf);
            w[j] = node.get_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2])*m_alphaf + node.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2])*(1-m_alphaf);
            a[j] = node.m_at*m_alpham + node.m_ap*(1-m_alpham);
            b[j] = node.get_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2])*m_alpham + node.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2])*(1-m_alpham);
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
			pt.m_r0 = evaluate(el, r0, s0, n);
			pt.m_rt = evaluate(el, r, s, n);

            // get the deformation gradient and determinant at intermediate time
            double Jt, Jp;
            mat3d Ft, Fp;
            Jt = defgrad(el, Ft, n);
            Jp = defgradp(el, Fp, n);
            pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
            pt.m_J = pt.m_F.det();
            mat3d Fi = pt.m_F.inverse();
            pt.m_L = (Ft - Fp)*Fi/dt;
            pt.m_v = evaluate(el, v, w, n);
            pt.m_a = evaluate(el, a, b, n);
            
            // evaluate strain energy at current time
            FEElasticMaterialPoint et = pt;
            et.m_F = Ft;
            et.m_J = Jt;
            pt.m_Wt = pme->StrainEnergyDensity(et);
            
            // calculate the stress at this material point
            pt.m_s = m_pMat->Stress(mp);
            
            // adjust stress for strain energy conservation
            if (m_alphaf == 0.5) {
                mat3ds D = pt.m_L.sym();
                double D2 = D.dotdot(D);
                if (D2 > 0)
                    pt.m_s += D*(((pt.m_Wt-pt.m_Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
            }
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
		lm[6*i  ] = id[m_dofU[0]];
		lm[6*i+1] = id[m_dofU[1]];
		lm[6*i+2] = id[m_dofU[2]];

		// next the shell displacement dofs
		lm[6*i+3] = id[m_dofSU[0]];
		lm[6*i+4] = id[m_dofSU[1]];
		lm[6*i+5] = id[m_dofSU[2]];

		// rigid rotational dofs
		lm[6*N + 3*i  ] = id[m_dofR[0]];
		lm[6*N + 3*i+1] = id[m_dofR[1]];
		lm[6*N + 3*i+2] = id[m_dofR[2]];
	}
}
