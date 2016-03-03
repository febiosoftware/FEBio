#include "stdafx.h"
#include "FEShellDomain.h"
#include "FEMesh.h"
#include "FEModel.h"

//-----------------------------------------------------------------------------
void FEShellDomain::InitElements()
{
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FEShellElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) el.GetMaterialPoint(j)->Init(false);
	}
}

//-----------------------------------------------------------------------------
void FEShellDomain::Reset()
{
	for (int i=0; i<(int) m_Elem.size(); ++i) m_Elem[i].Init(true);
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
        D[i] = el.m_D0[i];
    }
    
    double eta = el.gt(n);
    
    double* Mr = el.Hr(n);
    double* Ms = el.Hs(n);
    double* M  = el.H(n);
    
    // initialize covariant basis vectors
    g[0] = g[1] = g[2] = vec3d(0,0,0);
    
    for (i=0; i<neln; ++i)
    {
        g[0] += (r[i] + D[i]*eta/2)*Mr[i];
        g[1] += (r[i] + D[i]*eta/2)*Ms[i];
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
				ar << el.m_D0;

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
				ar >> el.m_D0;

				for (int j=0; j<el.GaussPoints(); ++j)
				{
					el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
					el.GetMaterialPoint(j)->Serialize(ar);
				}
			}
		}
	}
}
