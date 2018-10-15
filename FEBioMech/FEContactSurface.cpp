#include "stdafx.h"
#include "FEContactSurface.h"
#include "FECore/FEModel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include <assert.h>

//-----------------------------------------------------------------------------
FEContactSurface::FEContactSurface(FEModel* pfem) : FESurface(pfem), m_pfem(pfem)
{
	m_pSibling = 0; 
	m_dofX = -1;
	m_dofY = -1;
	m_dofZ = -1;
}

//-----------------------------------------------------------------------------
FEContactSurface::~FEContactSurface() { m_pSibling = 0; m_pContactInterface = 0; }

//-----------------------------------------------------------------------------
bool FEContactSurface::Init()
{
	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofX = dofs.GetDOF("x");
	m_dofY = dofs.GetDOF("y");
	m_dofZ = dofs.GetDOF("z");

	return FESurface::Init();
}

//-----------------------------------------------------------------------------
void FEContactSurface::SetSibling(FEContactSurface* ps) { m_pSibling = ps; }

//-----------------------------------------------------------------------------
void FEContactSurface::SetContactInterface(FEContactInterface* ps) { m_pContactInterface = ps; }

//-----------------------------------------------------------------------------
void FEContactSurface::GetVectorGap(int nface, vec3d& pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetContactTraction(int nface, vec3d& pt) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalVectorGap(int nface, vec3d* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactPressure(int nface, double* pg) {}

//-----------------------------------------------------------------------------
void FEContactSurface::GetNodalContactTraction(int nface, vec3d* pt) {}

//-----------------------------------------------------------------------------
vec3d FEContactSurface::GetContactForce() { return vec3d(0,0,0); }

//-----------------------------------------------------------------------------
double FEContactSurface::GetContactArea() { return 0; }

//-----------------------------------------------------------------------------
// Evaluate surface traction from the stress tensor of the attached solid element
void FEContactSurface::GetSurfaceTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    FEElement* e = el.m_elem[0];
    FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    if (se) {
        mat3ds si[FEElement::MAX_INTPOINTS];
        mat3ds so[FEElement::MAX_NODES];
        for (int i=0; i<se->GaussPoints(); ++i) {
            FEMaterialPoint* pt = se->GetMaterialPoint(i);
            FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
            if (ep)
                si[i] = ep->m_s;
            else
                si[i].zero();
        }
        // project stresses from integration points to nodes
        se->project_to_nodes(si, so);
        // only keep the stresses at the nodes of the contact face
        mat3ds sn[FEElement::MAX_NODES];
        for (int i=0; i<el.Nodes(); ++i)
            sn[i] = so[se->FindNode(el.m_node[i])];
        // evaluate tractions at integration points of that face
        vec3d t[FEElement::MAX_INTPOINTS];
        for (int i=0; i<el.GaussPoints(); ++i) {
            double *H = el.H(i);
            t[i] = vec3d(0,0,0);
            for (int j=0; j<el.Nodes(); ++j) {
                vec3d n = SurfaceNormal(el, j);
                t[i] += sn[j]*n*H[j];
            }
        }
        // now save the average traction on that face
        int ni = el.GaussPoints();
        pt = vec3d(0,0,0);
        for (int k=0; k<ni; ++k) pt += t[k];
        pt /= ni;
    }
    else
        pt = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
// Evaluate surface traction from the stress tensor of the attached solid element
void FEContactSurface::GetNodalSurfaceTraction(int nface, vec3d* pt)
{
    FESurfaceElement& el = Element(nface);
	FEElement* e = el.m_elem[0];
	FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    if (se) {
        mat3ds si[FEElement::MAX_INTPOINTS];
        mat3ds so[FEElement::MAX_NODES];
        for (int i=0; i<se->GaussPoints(); ++i) {
            FEMaterialPoint* pt = se->GetMaterialPoint(i);
            FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
            if (ep)
                si[i] = ep->m_s;
            else
                si[i].zero();
        }
        // project stresses from integration points to nodes
        se->project_to_nodes(si, so);
        // only keep the stresses at the nodes of the contact face
        mat3ds sn[FEElement::MAX_NODES];
        for (int i=0; i<el.Nodes(); ++i)
            sn[i] = so[se->FindNode(el.m_node[i])];
        // evaluate tractions at integration points of that face
        vec3d t[FEElement::MAX_INTPOINTS];
        for (int i=0; i<el.GaussPoints(); ++i) {
            double *H = el.H(i);
            t[i] = vec3d(0,0,0);
            for (int j=0; j<el.Nodes(); ++j) {
                vec3d n = SurfaceNormal(el, j);
                t[i] += sn[j]*n*H[j];
            }
        }
        // now project the traction from integration points back to nodes
        el.project_to_nodes(t, pt);
    }
    else
        for (int i=0; i<el.Nodes(); ++i) pt[i] = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
// Evaluate surface traction from the stress tensor of the attached solid element
void FEContactSurface::GetGPSurfaceTraction(int nface, vec3d* pt)
{
    FESurfaceElement& el = Element(nface);
	FEElement* e = el.m_elem[0];
	FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    if (se) {
        mat3ds si[FEElement::MAX_INTPOINTS];
        mat3ds so[FEElement::MAX_NODES];
        for (int i=0; i<se->GaussPoints(); ++i) {
            FEMaterialPoint* pt = se->GetMaterialPoint(i);
            FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
            if (ep)
                si[i] = ep->m_s;
            else
                si[i].zero();
        }
        // project stresses from integration points to nodes
        se->project_to_nodes(si, so);
        // only keep the stresses at the nodes of the contact face
        mat3ds sn[FEElement::MAX_NODES];
        for (int i=0; i<el.Nodes(); ++i)
            sn[i] = so[se->FindNode(el.m_node[i])];
        // evaluate tractions at integration points of that face
        for (int i=0; i<el.GaussPoints(); ++i) {
            double *H = el.H(i);
            pt[i] = vec3d(0,0,0);
            for (int j=0; j<el.Nodes(); ++j) {
                vec3d n = SurfaceNormal(el, j);
                pt[i] += sn[j]*n*H[j];
            }
        }
    }
    else
        for (int i=0; i<el.Nodes(); ++i) pt[i] = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
void FEContactSurface::GetStickStatus(int nface, double& pt) {}

//-----------------------------------------------------------------------------
void FEContactSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];
	}
}
