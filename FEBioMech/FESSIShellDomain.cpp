#include "stdafx.h"
#include "FESSIShellDomain.h"
#include "FECore/FEElemElemList.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModelParam.h>
#include <FECore/FEDataStream.h>

//-----------------------------------------------------------------------------
FESSIShellDomain::FESSIShellDomain(FEModel* pfem) : FEShellDomainNew(&pfem->GetMesh())
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
bool FESSIShellDomain::Init()
{
	FEShellDomain::Init();
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
		FESolidDomain& di = *psd[i];

		// identify all candidate elements by checking the tags on the nodes
		vector<int> elem; elem.reserve(di.Elements());
		for (int j=0; j<di.Elements(); ++j)
		{
			FESolidElement& el = di.Element(j);
			int ne = el.Nodes();
			for (int k=0; k<ne; ++k)
			{
				if (tag[el.m_node[k]] == 1)
				{
					elem.push_back(j);
					break;
				}
			}
		}

		// see if we can match any shells
		for (int j = 0; j<elem.size(); ++j)
		{
			FESolidElement& elj = di.Element(elem[j]);

			// loop over all its faces
			int nfaces = mesh.Faces(elj);
			for (int k=0; k<nfaces; ++k)
			{
				int nn = mesh.GetFace(elj, k, nf);

				// check all shell elements
				for (int l = 0; l<nelem; ++l)
				{
					FEShellElement& sel = Element(l);

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
							elj.m_bitfc.resize(elj.Nodes(), false);
							for (int n = 0; n<nn; ++n) {
								int m = elj.FindNode(nf[n]);
                                if (m > -1) elj.m_bitfc[m] = true;
							}
						}
                        else 
						{
                            // set solid element attached to shell front face
                            sel.m_elem[1] = elj.GetID();
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
                    if (sel) sldmn = dynamic_cast<FESolidDomain*>(sel->GetDomain());
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
                        if (el.GetDomain() == sldmn)
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

/*
//-----------------------------------------------------------------------------
//! Find interfaces between solid element faces and shell elements
//!
void FESSIShellDomain::FindSSI()
{
	// find out if there are solid domains in this model
	vector<FESolidDomain*> psd;
	FEMesh& mesh = *GetMesh();
	int ndom = mesh.Domains();
	for (int i = 0; i<ndom; ++i) {
		FEDomain& pdom = mesh.Domain(i);
		FESolidDomain* psdom = dynamic_cast<FESolidDomain*>(&pdom);
		if (psdom) psd.push_back(psdom);
	}
	size_t nsdom = psd.size();

	// if there are no solid domains we're done
	if (nsdom == 0) return;

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
				int nfaces = mesh.Faces(sel);
				for (int j = 0; j<nfaces; ++j) {
					nn = mesh.GetFace(sel, j, nf);

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
						vec3d n0 = mesh.Node(nf[0]).m_r0, n1 = mesh.Node(nf[1]).m_r0, n2 = mesh.Node(nf[2]).m_r0;
						vec3d nsld = (n1 - n0) ^ (n2 - n1);
						// get outward normal to shell face
						CoBaseVectors0(el, 0, g);
						vec3d nshl = g[2];
						nshl.unit();
						// compare normals
						if (nsld*nshl > 0) {
                            // set solid element attached to shell back face
                            el.m_elem[0] = sel.GetID();
							// store result
							sel.m_bitfc.resize(sel.Nodes(), false);
							for (int n = 0; n<nn; ++n) {
								int m = sel.FindNode(nf[n]);
                                if (m > -1) sel.m_bitfc[m] = true;
							}
						}
                        else {
                            // set solid element attached to shell front face
                            el.m_elem[1] = sel.GetID();
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
                    if (sel) sldmn = dynamic_cast<FESolidDomain*>(sel->GetDomain());
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
                        if (el.GetDomain() == sldmn) {
                            if (el.m_bitfc.size() == 0)
                                el.m_bitfc.resize(el.Nodes(), false);
                            int lid = el.FindNode(nid);
                            el.m_bitfc[lid] = true;
                        }
                    }
                }
            }
        }
    }
}
*/

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
        D[i] = ni.m_d0;
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
	CoBaseVectors0(el, n, gcov);

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
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_rt;
        D[i] = ni.m_d0 + ni.get_vec3d(m_dofx, m_dofy, m_dofz) - ni.get_vec3d(m_dofsx, m_dofsy, m_dofsz);
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
        D[i] = ni.m_d0 + ni.m_rp - ni.m_r0 - ni.get_vec3d(m_dofsxp, m_dofsyp, m_dofszp);
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
    CoBaseVectors(el, n, gt);
    CoBaseVectorsP(el, n, gp);
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
        D[i] = ni.m_d0 + ni.get_vec3d(m_dofx, m_dofy, m_dofz) - ni.get_vec3d(m_dofsx, m_dofsy, m_dofsz);
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
// jacobian with respect to intermediate time frame
double FESSIShellDomain::detJ(FEShellElement& el, int n, const double alpha)
{
    vec3d gcovt[3], gcovp[3], gcov[3];
    CoBaseVectors(el, n, gcovt);
    CoBaseVectorsP(el, n, gcovp);
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
    CoBaseVectors(el, n, gcov);
    ContraBaseVectors0(el, n, Gcnt);
    
    F = (gcov[0] & Gcnt[0]) + (gcov[1] & Gcnt[1]) + (gcov[2] & Gcnt[2]);
    double J = F.det();
    if (J <= 0) throw NegativeJacobian(el.GetID(), n, J, &el);
    
    return J;
}

//-----------------------------------------------------------------------------
double FESSIShellDomain::defgradp(FEShellElement& el, mat3d& F, int n)
{
    vec3d gcov[3], Gcnt[3];
    CoBaseVectorsP(el, n, gcov);
    ContraBaseVectors0(el, n, Gcnt);
    
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

void FESSIShellDomain::Update(const FETimeInfo& tp)
{
	int NS = Elements();
	FEMesh& mesh = *GetMesh();
	for (int i = 0; i<NS; ++i)
	{
		FEShellElementNew& e = ShellElement(i);
		int n = e.Nodes();
		for (int j = 0; j<n; ++j)
		{
			FENode& nj = mesh.Node(e.m_node[j]);
			vec3d D = nj.m_d0 + nj.get_vec3d(m_dofx, m_dofy, m_dofz) - nj.get_vec3d(m_dofsx, m_dofsy, m_dofsz);
			double h = D.norm();

			e.m_ht[j] = h;
		}
	}
}


//=================================================================================================
template <class T> void _writeIntegratedElementValueT(FESSIShellDomain& dom, FEValuator<T>& var, FEDataStream& ar)
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
			T vj = var(mp);
			double detJ = dom.detJ0(el, j)*gw[j];
			ew += vj*detJ;
		}
		ar << ew;
	}
}

void writeIntegratedElementValue(FESSIShellDomain& dom, FEValuator<double>& var, FEDataStream& ar) { _writeIntegratedElementValueT<double>(dom, var, ar); }
void writeIntegratedElementValue(FESSIShellDomain& dom, FEValuator<vec3d>&  var, FEDataStream& ar) { _writeIntegratedElementValueT<vec3d >(dom, var, ar); }
