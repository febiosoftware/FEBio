#include "stdafx.h"
#include "FEShellDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"
#include <FESolidDomain.h>

//-----------------------------------------------------------------------------
void FEShellDomain::Create(int nelems, int elemType)
{
	m_Elem.resize(nelems);
	if (elemType != -1)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(elemType);
}

//-----------------------------------------------------------------------------
void FEShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FEShellElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Update(timeInfo);
	}
    
    // check for solid-shell interfaces
    if (!m_binit) {
        FindSSI();
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FEShellElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Init();
	}
}

//-----------------------------------------------------------------------------
//! calculates covariant basis vectors at an integration point
void FEShellDomain::CoBaseVectors0(FEShellElement& el, int n, vec3d g[3])
{
    int i;
    
    int neln = el.Nodes();
    
    // current nodal coordinates and directors
    vec3d r[FEElement::MAX_NODES], D[FEElement::MAX_NODES];
    for (i=0; i<neln; ++i)
    {
        FENode& ni = m_pMesh->Node(el.m_node[i]);
        r[i] = ni.m_r0;
        D[i] = ni.m_d0;
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
void FEShellDomain::ContraBaseVectors0(FEShellElement& el, int n, vec3d gcnt[3])
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    mat3d Ji = J.inverse();
    
    gcnt[0] = vec3d(Ji(0,0),Ji(0,1),Ji(0,2));
    gcnt[1] = vec3d(Ji(1,0),Ji(1,1),Ji(1,2));
    gcnt[2] = vec3d(Ji(2,0),Ji(2,1),Ji(2,2));
    
}

//-----------------------------------------------------------------------------
double FEShellDomain::invjac0(FEShellElement& el, double Ji[3][3], int n)
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
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
//! Calculate jacobian with respect to reference frame
double FEShellDomain::detJ0(FEShellElement &el, int n)
{
    vec3d gcov[3];
    CoBaseVectors0(el, n, gcov);
    
    mat3d J = mat3d(gcov[0].x, gcov[1].x, gcov[2].x,
                    gcov[0].y, gcov[1].y, gcov[2].y,
                    gcov[0].z, gcov[1].z, gcov[2].z);
    return J.det();
}

//-----------------------------------------------------------------------------
void FEShellDomain::Serialize(DumpStream &ar)
{
	if (ar.IsShallow())
	{
		int NEL = (int) m_Elem.size();
		for (int i=0; i<NEL; ++i)
		{
			FEShellElement& el = m_Elem[i];
			int nint = el.GaussPoints();
			for (int j=0; j<nint; ++j) el.GetMaterialPoint(j)->Serialize(ar);
		}
	}
	else
	{
		if (ar.IsSaving())
		{
			ar << m_Node;
		
			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FEShellElement& el = m_Elem[i];
				ar << el.Type();

				ar << el.GetMatID();
				ar << el.GetID();
				ar << el.m_node;

				ar << el.m_h0;

				for (int j=0; j<el.GaussPoints(); ++j) el.GetMaterialPoint(j)->Serialize(ar);
			}
		}
		else
		{
			int n, mat, nid;

			ar >> m_Node;

			FEModel& fem = ar.GetFEModel();
			FEMaterial* pmat = GetMaterial();
			assert(pmat);

			for (size_t i=0; i<m_Elem.size(); ++i)
			{
				FEShellElement& el = m_Elem[i];
				ar >> n;

				el.SetType(n);

				ar >> mat; el.SetMatID(mat);
				ar >> nid; el.SetID(nid);
				ar >> el.m_node;

				ar >> el.m_h0;

				for (int j=0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Find interfaces between solid element faces and shell elements
//!
void FEShellDomain::FindSSI()
{
    // find out if there are solid domains in this model
    vector<FESolidDomain*> psd;
    FEMesh* mesh = GetMesh();
    int ndom = mesh->Domains();
    for (int i=0; i<ndom; ++i) {
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
    for (int i=0; i<nelem; ++i) {
        FEShellElement& el = *dynamic_cast<FEShellElement*>(&ElementRef(i));
        
        // check all solid domains
        for (int k=0; k<nsdom; ++k) {
            
            // check each solid element in this domain
            int nselem = psd[k]->Elements();
            for (int l=0; l<nselem; ++l) {
                FEElement& sel = psd[k]->ElementRef(l);
                
                // check all faces of this solid element
                int nfaces = m.Faces(sel);
                for (int j=0; j<nfaces; ++j) {
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
                            for (int n=0; n<nn; ++n) {
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
