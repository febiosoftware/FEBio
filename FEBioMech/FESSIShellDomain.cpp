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
#include "FESSIShellDomain.h"
#include "FECore/FEElemElemList.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModelParam.h>
#include <FECore/FEDataStream.h>
#include <FECore/log.h>
#include <FECore/FEMesh.h>
#include "FEBioMech.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESSIShellDomain, FEShellDomainNew)
	ADD_PARAMETER(m_bnodalnormals, "shell_normal_nodal");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESSIShellDomain::FESSIShellDomain(FEModel* pfem) : FEShellDomainNew(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofR(pfem)
{
    m_bnodalnormals = true;

    // TODO: Can this be done in Init, since there is no error checking
    if (pfem)
    {
        m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
        m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
        m_dofR.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
    }
}

//-----------------------------------------------------------------------------
bool FESSIShellDomain::Init()
{
	FEShellDomain::Init();
    FindSSI();

	// check for initially inverted shells
	try {
		for (int i = 0; i < Elements(); ++i)
		{
			FEShellElementNew& el = ShellElement(i);
			int neln = el.Nodes();
			int nint = el.GaussPoints();
			el.m_E.resize(nint, mat3ds(0, 0, 0, 0, 0, 0));

			for (int n = 0; n < nint; ++n)
			{
				double J0 = detJ0(el, n);
			}
		}
	}
	catch (NegativeJacobian e)
	{
		feLogError("Zero or negative jacobians detected at integration points of domain: %s\n", GetName().c_str());
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESSIShellDomain::Serialize(DumpStream& ar)
{
	FEShellDomainNew::Serialize(ar);
	ar & m_bnodalnormals;
}

//-----------------------------------------------------------------------------
//! Calculate all shell normals (i.e. the shell directors).
//! And find shell nodes
void FESSIShellDomain::InitShells()
{
    FEShellDomain::InitShells();
    
    if (!m_bnodalnormals) {
        FEMesh& mesh = *GetMesh();
        for (int i = 0; i<Elements(); ++i)
        {
            FEShellElementNew& el = ShellElement(i);
            int ne = el.Nodes();
            vec3d gr(0,0,0), gs(0,0,0);
            double Mr[FEElement::MAX_NODES], Ms[FEElement::MAX_NODES];
            el.shape_deriv(Mr, Ms, 0, 0);
            for (int j = 0; j<ne; ++j)
            {
                gr += mesh.Node(el.m_node[j]).m_r0*Mr[j];
                gs += mesh.Node(el.m_node[j]).m_r0*Ms[j];
            }
            vec3d d0 = gr ^ gs;
            d0.unit();
            for (int j = 0; j<ne; ++j)
            {
                el.m_d0[j] = d0 * el.m_h0[j];
            }
        }
    }

	FEMesh& mesh = *GetMesh();
	for (int i = 0; i < Elements(); ++i)
	{
		FEShellElementNew& el = ShellElement(i);
		int ni = el.GaussPoints();
		vec3d g0[3];
		for (int j = 0; j < ni; ++j)
		{
			//NOTE: calculate covariant first since contravariant depends on it!
			CoBaseVectors0(el, j, g0);
			el.m_gt[0][j] = el.m_gp[0][j] = el.m_g0[0][j] = g0[0];
			el.m_gt[1][j] = el.m_gp[1][j] = el.m_g0[1][j] = g0[1];
			el.m_gt[2][j] = el.m_gp[2][j] = el.m_g0[2][j] = g0[2];

			ContraBaseVectors0(el, j, g0);
			el.m_G0[0][j] = el.m_Gt[0][j] = g0[0];
			el.m_G0[1][j] = el.m_Gt[1][j] = g0[1];
			el.m_G0[2][j] = el.m_Gt[2][j] = g0[2];
		}
	}
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
	FEMesh& mesh = *GetMesh();
	int ndom = mesh.Domains();
	for (int i = 0; i<ndom; ++i)
	{
		FEDomain& pdom = mesh.Domain(i);
		FESolidDomain* psdom = dynamic_cast<FESolidDomain*>(&pdom);
		if (psdom) psd.push_back(psdom);
	}
	size_t nsdom = psd.size();

	// if there are no solid domains we're done
	if (nsdom == 0) return;

	// tag all nodes that belong to this shell domain
	vector<int> tag(mesh.Nodes(), 0);
	for (int i=0; i<Nodes(); ++i) tag[m_Node[i]] = 1;

	int nelem = Elements();
	int nf[9];
	vec3d g[3];

	// loop over all solid domains
	for (int i=0; i<nsdom; ++i)
	{
		// get the next domain
		FESolidDomain& di = *psd[i];

		// see if there are any potential candidates
		// NOTE: We can't use the domain's node list since it might not be initialized yet
		//       since shell domains are initialized before solid domains
		bool hasNodes = false;
		for (int j = 0; j < di.Elements(); ++j)
		{
			FEElement& el = di.Element(j);
			int ne = el.Nodes();
			for (int k = 0; k < ne; ++k)
			{
				FENode& node = mesh.Node(el.m_node[k]);
				if (tag[el.m_node[k]] == 1)
				{
					hasNodes = true;
					break;
				}
			}
		}

		// if the solid domain shares nodes with this shell domain,
		// there might be matches in this domain
		if (hasNodes)
		{
			// build a node-element list for faster lookup
			FENodeElemList NEL; NEL.Create(di);

			// check all shell elements
			for (int l = 0; l < nelem; ++l)
			{
				FEShellElement& sel = Element(l);

				int n0 = sel.m_node[0];
				int nval = NEL.Valence(n0);

				// see if we can match any shells
				FEElement** ppe = NEL.ElementList(n0);
				for (int j = 0; j < nval; ++j)
				{
					FESolidElement& elj = dynamic_cast<FESolidElement&>(*ppe[j]);

					// loop over all its faces
					int nfaces = elj.Faces();
					for (int k = 0; k < nfaces; ++k)
					{
						int nn = elj.GetFace(k, nf);

						// see if it matches this face
						bool found = false;
						if (nn == sel.Nodes())
						{
							switch (nn)
							{
							case 3: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2])) found = true; break;
							case 4: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2]) && sel.HasNode(nf[3])) found = true; break;
							case 6: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2])) found = true; break;
							case 7: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2])) found = true; break;
							case 8: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2]) && sel.HasNode(nf[3])) found = true; break;
							case 9: if (sel.HasNode(nf[0]) && sel.HasNode(nf[1]) && sel.HasNode(nf[2]) && sel.HasNode(nf[3])) found = true; break;
							default:
								assert(false);
							}
						}

						// see if we found a match
						if (found)
						{
							// check interface side
							// get outward normal to solid element face
							vec3d n0 = mesh.Node(nf[0]).m_r0, n1 = mesh.Node(nf[1]).m_r0, n2 = mesh.Node(nf[2]).m_r0;
							vec3d nsld = (n1 - n0) ^ (n2 - n1);
							// get outward normal to shell face
							CoBaseVectors0(sel, 0, g);
							vec3d nshl = g[2];

							// compare normals
							if (nsld*nshl > 0)
							{
								// set solid element attached to shell back face
								sel.m_elem[0] = elj.GetID();

								// store result
								if (m_bnodalnormals) {
									elj.m_bitfc.resize(elj.Nodes(), false);
									for (int n = 0; n < nn; ++n) {
										int m = elj.FindNode(nf[n]);
										if (m > -1) elj.m_bitfc[m] = true;
									}
								}
							}
							else
							{
								// set solid element attached to shell front face
								sel.m_elem[1] = elj.GetID();
							}

							// the same element cannot be both front and back of a shell
							// so we can leave the loop over the faces
							break;
						}
					}
				}
			}
		}
	}
    
    // check for elements that only have one or two shell nodes
    // but don't share a whole face

    // create the node element list
    FENodeElemList NEL;
    NEL.Create(mesh);
    
    for (int i = 0; i<ndom; ++i) {
        FEDomain& pdom = mesh.Domain(i);
        FEShellDomain* psdom = dynamic_cast<FEShellDomain*>(&pdom);
        if (psdom) {
            // find the solid domain attached to the back of these shells
            FESolidDomain* sldmn = nullptr;
            for (int j=0; j<psdom->Elements(); ++j) {
                FEShellElement& el1 = psdom->Element(j);
                // identify solid domain at back of shell domain
                if (el1.m_elem[0] != -1) {
                    FEElement* sel = mesh.FindElementFromID(el1.m_elem[0]);
                    if (sel) sldmn = dynamic_cast<FESolidDomain*>(sel->GetMeshPartition());
                    break;
                }
            }
            if (sldmn) {
                // for each node in this shell domain, check the solid elements it belongs to
                for (int j=0; j<psdom->Nodes(); ++j) {
                    FENode& node = psdom->Node(j);
                    int nid = node.GetID() - 1;
                    int nval = NEL.Valence(nid);
                    FEElement** pe = NEL.ElementList(nid);
                    for (int k=0; k<nval; ++k)
                    {
                        // get the element
                        FEElement& el = *pe[k];
                        // check that it belongs to the solid domain at the back of the shell domain
                        if ((el.GetMeshPartition() == sldmn) && m_bnodalnormals)
						{
							FESolidElement& sel = dynamic_cast<FESolidElement&>(el);
                            if (sel.m_bitfc.size() == 0)
                                sel.m_bitfc.resize(el.Nodes(), false);
                            int lid = sel.FindNode(nid);
                            sel.m_bitfc[lid] = true;
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
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
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
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors0(FEShellElement& el, double r, double s, double t, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r0[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i = 0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r0[i] = ni.m_r0;
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
    }
    
    double eta = t;
    
    double M[FEElement::MAX_NODES];
    double Mr[FEElement::MAX_NODES];
    double Ms[FEElement::MAX_NODES];
    el.shape_fnc(M, r, s);
    el.shape_deriv(Mr, Ms, r, s);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0, 0, 0);
    
    for (i = 0; i<neln; ++i)
    {
        g[0] += (r0[i] - D[i] * (1 - eta) / 2)*Mr[i];
        g[1] += (r0[i] - D[i] * (1 - eta) / 2)*Ms[i];
        g[2] += D[i] * (M[i] / 2);
    }
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors0(FEShellElement& el, int n, vec3d gcnt[3])
{
	vec3d gcov[3];
//	CoBaseVectors0(el, n, gcov);
	gcov[0] = el.m_g0[0][n];
	gcov[1] = el.m_g0[1][n];
	gcov[2] = el.m_g0[2][n];

	mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
		gcov[0].y, gcov[1].y, gcov[2].y,
		gcov[0].z, gcov[1].z, gcov[2].z);
	mat3d Ji = J.inverse();

	gcnt[0] = vec3d(Ji(0, 0), Ji(0, 1), Ji(0, 2));
	gcnt[1] = vec3d(Ji(1, 0), Ji(1, 1), Ji(1, 2));
	gcnt[2] = vec3d(Ji(2, 0), Ji(2, 1), Ji(2, 2));

}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors0(FEShellElement& el, double r, double s, double t, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors0(el, r, s, t, gcov);
    
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
//	CoBaseVectors0(el, n, gcov);
	gcov[0] = el.m_g0[0][n];
	gcov[1] = el.m_g0[1][n];
	gcov[2] = el.m_g0[2][n];

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
double FESSIShellDomain::invjac0(FEShellElement& el, double Ji[3][3], double r, double s, double t)
{
    vec3d gcov[3];
    CoBaseVectors0(el, r, s, t, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    
    double det = J.det();
    
    // make sure the determinant is positive
    if (det <= 0) throw NegativeJacobian(el.GetID(), -1, det);
    
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
//	CoBaseVectors0(el, n, gcov);
	gcov[0] = el.m_g0[0][n];
	gcov[1] = el.m_g0[1][n];
	gcov[2] = el.m_g0[2][n];

	mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
		gcov[0].y, gcov[1].y, gcov[2].y,
		gcov[0].z, gcov[1].z, gcov[2].z);
	return J.det();
}

//-----------------------------------------------------------------------------
// jacobian with respect to current frame
double FESSIShellDomain::detJ0(FEShellElement& el, double r, double s, double t)
{
    vec3d Gcov[3];
    CoBaseVectors0(el, r, s, t, Gcov);
    
    mat3d J = mat3d(Gcov[0].x, Gcov[1].x, Gcov[2].x,
                    Gcov[0].y, Gcov[1].y, Gcov[2].y,
                    Gcov[0].z, Gcov[1].z, Gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors(FEShellElement& el, int n, vec3d g[3])
{
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (int i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rt;
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
        D[i] +=  ni.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]) - ni.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
    }
    
    double eta = el.gt(n);
    
    double* Mr = el.Hr(n);
    double* Ms = el.Hs(n);
    double* M  = el.H(n);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (int i=0; i<neln; ++i)
    {
        g[0] += (r[i] - D[i]*(1-eta)/2)*Mr[i];
        g[1] += (r[i] - D[i]*(1-eta)/2)*Ms[i];
        g[2] += D[i]*(M[i]/2);
    }
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point at previous time
void FESSIShellDomain::CoBaseVectorsP(FEShellElement& el, int n, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // previous time nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rp;
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
        D[i] += ni.m_rp - ni.m_r0 - ni.get_vec3d_prev(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
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
//! calculates covariant basis vectors at an integration point at intermediate time
void FESSIShellDomain::CoBaseVectors(FEShellElement& el, int n, vec3d g[3], const double alpha)
{
    vec3d gt[3], gp[3];
//    CoBaseVectors(el, n, gt);
	gt[0] = el.m_gt[0][n];
	gt[1] = el.m_gt[1][n];
	gt[2] = el.m_gt[2][n];
//    CoBaseVectorsP(el, n, gp);
	gp[0] = el.m_gp[0][n];
	gp[1] = el.m_gp[1][n];
	gp[2] = el.m_gp[2][n];
    for (int i=0; i<3; ++i)
        g[i] = gt[i]*alpha + gp[i]*(1-alpha);
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d rt[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        rt[i] = ni.m_rt;
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
        D[i] += ni.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]) - ni.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
    }
    
    double eta = t;
    
    double M[FEElement::MAX_NODES];
    double Mr[FEElement::MAX_NODES];
    double Ms[FEElement::MAX_NODES];
    el.shape_fnc(M, r, s);
    el.shape_deriv(Mr, Ms, r, s);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (i=0; i<neln; ++i)
    {
        g[0] += (rt[i] - D[i]*(1-eta)/2)*Mr[i];
        g[1] += (rt[i] - D[i]*(1-eta)/2)*Ms[i];
        g[2] += D[i]*(M[i]/2);
    }
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FESSIShellDomain::CoBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3], const double alpha)
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d rt[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        rt[i] = ni.m_rt*alpha + ni.m_rp*(1-alpha);
        D[i] = m_bnodalnormals ? ni.m_d0 : el.m_d0[i];
        D[i] += (ni.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]) - ni.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]))*alpha
        + (ni.get_vec3d_prev(m_dofU[0], m_dofU[1], m_dofU[2]) - ni.get_vec3d_prev(m_dofSU[0], m_dofSU[1], m_dofSU[2]))*(1-alpha);
    }
    
    double eta = t;
    
    double M[FEElement::MAX_NODES];
    double Mr[FEElement::MAX_NODES];
    double Ms[FEElement::MAX_NODES];
    el.shape_fnc(M, r, s);
    el.shape_deriv(Mr, Ms, r, s);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (i=0; i<neln; ++i)
    {
        g[0] += (rt[i] - D[i]*(1-eta)/2)*Mr[i];
        g[1] += (rt[i] - D[i]*(1-eta)/2)*Ms[i];
        g[2] += D[i]*(M[i]/2);
    }
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors(FEShellElement& el, int n, vec3d gcnt[3])
{
    vec3d gcov[3];
//    CoBaseVectors(el, n, gcov);
	gcov[0] = el.m_gt[0][n];
	gcov[1] = el.m_gt[1][n];
	gcov[2] = el.m_gt[2][n];
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
    
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors(FEShellElement& el, int n, vec3d gcnt[3], const double alpha)
{
    vec3d gcov[3];
    CoBaseVectors(el, n, gcov, alpha);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
    
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors(FEShellElement& el, double r, double s, double t, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors(el, r, s, t, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
    
}

//-----------------------------------------------------------------------------
//! calculates contravariant basis vectors at an integration point
void FESSIShellDomain::ContraBaseVectors(FEShellElement& el, double r, double s, double t, vec3d gcnt[3], const double alpha)
{
    vec3d gcov[3];
    CoBaseVectors(el, r, s, t, gcov, alpha);
    
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
//    CoBaseVectors(el, n, gcov);
	gcov[0] = el.m_gt[0][n];
	gcov[1] = el.m_gt[1][n];
	gcov[2] = el.m_gt[2][n];
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
// jacobian with respect to intermediate time frame
double FESSIShellDomain::detJ(FEShellElement& el, int n, const double alpha)
{
    vec3d gcovt[3], gcovp[3], gcov[3];
//    CoBaseVectors(el, n, gcovt);
	gcovt[0] = el.m_gt[0][n];
	gcovt[1] = el.m_gt[1][n];
	gcovt[2] = el.m_gt[2][n];
//    CoBaseVectorsP(el, n, gcovp);
	gcovp[0] = el.m_gp[0][n];
	gcovp[1] = el.m_gp[1][n];
	gcovp[2] = el.m_gp[2][n];
	for (int i=0; i<3; ++i)
        gcov[i] = gcovt[i]*alpha + gcovp[i]*(1-alpha);

    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
// jacobian with respect to current frame
double FESSIShellDomain::detJ(FEShellElement& el, double r, double s, double t)
{
    vec3d gcov[3];
    CoBaseVectors(el, r, s, t, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
double FESSIShellDomain::defgrad(FEShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
//    CoBaseVectors(el, n, gcov);
	gcov[0] = el.m_gt[0][n];
	gcov[1] = el.m_gt[1][n];
	gcov[2] = el.m_gt[2][n];
//	ContraBaseVectors0(el, n, Gcnt);
	Gcnt[0] = el.m_G0[0][n];
	Gcnt[1] = el.m_G0[1][n];
	Gcnt[2] = el.m_G0[2][n];

    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    double J = F.det();
    if (J <= 0) throw NegativeJacobian(el.GetID(), n, J, &el);
    
    return J;
}

//-----------------------------------------------------------------------------
double FESSIShellDomain::defgradp(FEShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
//    CoBaseVectorsP(el, n, gcov);
	gcov[0] = el.m_gp[0][n];
	gcov[1] = el.m_gp[1][n];
	gcov[2] = el.m_gp[2][n];
//    ContraBaseVectors0(el, n, Gcnt);
	Gcnt[0] = el.m_G0[0][n];
	Gcnt[1] = el.m_G0[1][n];
	Gcnt[2] = el.m_G0[2][n];

    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    double J = F.det();
    if (J <= 0) throw NegativeJacobian(el.GetID(), n, J, &el);
    
    return J;
}

//-----------------------------------------------------------------------------
double FESSIShellDomain::defgrad(FEShellElement& el, mat3d& F, double r, double s, double t)
{
    vec3d gcov[3], Gcnt[3];
    CoBaseVectors(el, r, s, t, gcov);
    ContraBaseVectors0(el, r, s, t, Gcnt);
    
    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    double J = F.det();
    if (J <= 0) throw NegativeJacobian(el.GetID(), -1, J, &el);
    
    return J;
}

//-----------------------------------------------------------------------------
//! evaluate a vector function over the shell
vec3d FESSIShellDomain::evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n)
{
	vec3d gcnt[3];
//	ContraBaseVectors(el, n, gcnt);
	gcnt[0] = el.m_Gt[0][n];
	gcnt[1] = el.m_Gt[1][n];
	gcnt[2] = el.m_Gt[2][n];

	double eta = el.gt(n);
	const double* M = el.H(n);
	vec3d v(0, 0, 0);
	int neln = el.Nodes();
	for (int i=0; i<neln; ++i)
    {
        double Mu = (1+eta)/2*M[i];
        double Md = (1-eta)/2*M[i];
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
//    CoBaseVectors(el, n, gcov);
	gcov[0] = el.m_gt[0][n];
	gcov[1] = el.m_gt[1][n];
	gcov[2] = el.m_gt[2][n];
    
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
	vec3d gcnt[3];
//    ContraBaseVectors(el, n, gcnt);
	gcnt[0] = el.m_Gt[0][n];
	gcnt[1] = el.m_Gt[1][n];
	gcnt[2] = el.m_Gt[2][n];

	double eta = el.gt(n);
	const double* M = el.H(n);
	double p = 0;
	int neln = el.Nodes();
    for (int i=0; i<neln; ++i)
    {
        double Mu = (1+eta)/2*M[i];
        double Md = (1-eta)/2*M[i];
        p += Mu*pn[i] + Md*dpn[i];
    }
    
    return p;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of a scalar function over the shell
vec3d FESSIShellDomain::gradient(FEShellElement& el, double* pn, double* dpn, int n)
{
	vec3d gcnt[3];
//    ContraBaseVectors(el, n, gcnt);
	gcnt[0] = el.m_Gt[0][n];
	gcnt[1] = el.m_Gt[1][n];
	gcnt[2] = el.m_Gt[2][n];

	double eta = el.gt(n);

	const double* Mr = el.Hr(n);
	const double* Ms = el.Hs(n);
	const double* M = el.H(n);

	vec3d gradp(0, 0, 0);
	int neln = el.Nodes();
    for (int i=0; i<neln; ++i)
    {
        vec3d gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
		vec3d gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
		vec3d gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        gradp += gradMu*pn[i] + gradMd*dpn[i];
    }
    
    return gradp;
}

//-----------------------------------------------------------------------------
//! evaluate a scalar function over the shell
double FESSIShellDomain::evaluate(FEShellElement& el, vector<double> pn, vector<double> dpn, int n)
{
	vec3d gcnt[3];
//    ContraBaseVectors(el, n, gcnt);
	gcnt[0] = el.m_Gt[0][n];
	gcnt[1] = el.m_Gt[1][n];
	gcnt[2] = el.m_Gt[2][n];

	double eta = el.gt(n);
	const double* M = el.H(n);

	double p = 0;
	int neln = el.Nodes();
    for (int i=0; i<neln; ++i)
    {
        double Mu = (1+eta)/2*M[i];
        double Md = (1-eta)/2*M[i];
        p += Mu*pn[i] + Md*dpn[i];
    }
    
    return p;
}

//-----------------------------------------------------------------------------
//! calculate the gradient of a scalar function over the shell
vec3d FESSIShellDomain::gradient(FEShellElement& el, vector<double> pn, vector<double> dpn, int n)
{
	vec3d gcnt[3];
//    ContraBaseVectors(el, n, gcnt);
	gcnt[0] = el.m_Gt[0][n];
	gcnt[1] = el.m_Gt[1][n];
	gcnt[2] = el.m_Gt[2][n];

	double eta = el.gt(n);
	const double* Mr = el.Hr(n);
	const double* Ms = el.Hs(n);
	const double* M = el.H(n);

    int neln = el.Nodes();
	vec3d gradp(0, 0, 0);
	for (int i=0; i<neln; ++i)
    {
        vec3d gradM = gcnt[0]*Mr[i] + gcnt[1]*Ms[i];
        vec3d gradMu = (gradM*(1+eta) + gcnt[2]*M[i])/2;
        vec3d gradMd = (gradM*(1-eta) - gcnt[2]*M[i])/2;
        gradp += gradMu*pn[i] + gradMd*dpn[i];
    }
    
    return gradp;
}

void FESSIShellDomain::Update(const FETimeInfo& tp)
{
	int NS = Elements();
	FEMesh& mesh = *GetMesh();
	int NE = Elements();
#pragma omp parallel for
	for (int i = 0; i < NE; ++i)
	{
		FEShellElement& e = Element(i);

		int n = e.Nodes();
		for (int j = 0; j<n; ++j)
		{
			FENode& nj = mesh.Node(e.m_node[j]);
			vec3d D = nj.m_dt;
			double h = D.norm();

			e.m_ht[j] = h;
		}

		vec3d g[3];
		int ni = e.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			CoBaseVectors(e, j, g);
			e.m_gt[0][j] = g[0];
			e.m_gt[1][j] = g[1];
			e.m_gt[2][j] = g[2];

			CoBaseVectorsP(e, j, g);
			e.m_gp[0][j] = g[0];
			e.m_gp[1][j] = g[1];
			e.m_gp[2][j] = g[2];

			// NOTE: calculate covariant vectors first since contravariant depends on them
			ContraBaseVectors(e, j, g);
			e.m_Gt[0][j] = g[0];
			e.m_Gt[1][j] = g[1];
			e.m_Gt[2][j] = g[2];
		}
	};
}

//=================================================================================================
template <class T> void _writeIntegratedElementValueT(FESSIShellDomain& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint& mp)> fnc)
{
	for (int i = 0; i<dom.Elements(); ++i)
	{
		FEShellElement& el = dom.Element(i);
		double* gw = el.GaussWeights();

		// integrate
		T ew(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			T vj = fnc(mp);
			double detJ = dom.detJ0(el, j)*gw[j];
			ew += vj*detJ;
		}
		ar << ew;
	}
}

void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint& mp)> fnc) { _writeIntegratedElementValueT<double>(dom, ar, fnc); }
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint& mp)> fnc) { _writeIntegratedElementValueT<vec3d >(dom, ar, fnc); }
