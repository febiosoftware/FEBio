//
//  FETiedFluidInterface.cpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/6/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FETiedFluidInterface.h"
#include "FEFluid.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedFluidInterface, FEContactInterface)
	ADD_PARAMETER(m_blaugon  , "laugon"             );
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_ptol     , "ptol"               );
	ADD_PARAMETER(m_epst     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_epsn     , "pressure_penalty"   );
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedFluidSurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_vg  = vec3d(0,0,0);
    m_nu  = vec3d(0,0,0);
    m_rs  = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_tv  = vec3d(0,0,0);
    m_Lmp = 0.0;
    m_epst= 1.0;
    m_epsn= 1.0;
    m_pg  = 0.0;
    m_vn  = 0.0;
    m_pme = (FESurfaceElement*)0;
}

//-----------------------------------------------------------------------------
// FETiedFluidSurface
//-----------------------------------------------------------------------------

FETiedFluidSurface::FETiedFluidSurface(FEModel* pfem) : FEContactSurface(pfem)
{
    m_pfem = pfem;
}

//-----------------------------------------------------------------------------
bool FETiedFluidSurface::Init()
{
    // initialize surface data first
    if (FEContactSurface::Init() == false) return false;
    
    DOFS& dofs = GetFEModel()->GetDOFS();
    m_dofWX = dofs.GetDOF("wx");
    m_dofWY = dofs.GetDOF("wy");
    m_dofWZ = dofs.GetDOF("wz");
    m_dofEF = dofs.GetDOF("ef");

    return true;
}


//-----------------------------------------------------------------------------
void FETiedFluidSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    lm.resize(N*4);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        vector<int>& id = node.m_ID;
        
        lm[4*i  ] = id[m_dofWX];
        lm[4*i+1] = id[m_dofWY];
        lm[4*i+2] = id[m_dofWZ];
        lm[4*i+3] = id[m_dofEF];
    }
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FETiedFluidSurface::CreateMaterialPoint()
{
	return new FETiedFluidSurface::Data;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::Serialize(DumpStream& ar)
{
    if (ar.IsShallow())
    {
        if (ar.IsSaving())
        {
            for (int i=0; i<Elements(); ++i)
            {
				FESurfaceElement& el = Element(i);
                int nint = el.GaussPoints();
                for (int j=0; j<nint; ++j)
                {
                    Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
                    ar << d.m_Lmd;
                    ar << d.m_Gap;
                    ar << d.m_vg;
                    ar << d.m_pg;
                    ar << d.m_Lmp;
                    ar << d.m_tv;
                    ar << d.m_vn;
                }
            }
        }
        else
        {
			for (int i = 0; i<Elements(); ++i)
			{
				FESurfaceElement& el = Element(i);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar >> d.m_Lmd;
                    ar >> d.m_Gap;
                    ar >> d.m_vg;
                    ar >> d.m_pg;
                    ar >> d.m_Lmp;
                    ar >> d.m_tv;
                    ar >> d.m_vn;
                }
            }
        }
    }
    else
    {
        // We can serialize the base-class data
        FEContactSurface::Serialize(ar);
        
        // And finally, we serialize the surface data
        if (ar.IsSaving())
        {
			for (int i = 0; i<Elements(); ++i)
			{
				FESurfaceElement& el = Element(i);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar << d.m_Gap;
                    ar << d.m_vg;
                    ar << d.m_nu;
                    ar << d.m_rs;
                    ar << d.m_Lmd;
                    ar << d.m_Lmp;
                    ar << d.m_epst;
                    ar << d.m_epsn;
                    ar << d.m_pg;
                    ar << d.m_tv;
                    ar << d.m_vn;
                }
            }
        }
        else
        {
			for (int i = 0; i<Elements(); ++i)
			{
				FESurfaceElement& el = Element(i);
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					Data& d = static_cast<Data&>(*el.GetMaterialPoint(j));
					ar >> d.m_Gap;
                    ar >> d.m_vg;
                    ar >> d.m_nu;
                    ar >> d.m_rs;
                    ar >> d.m_Lmd;
                    ar >> d.m_Lmp;
                    ar >> d.m_epst;
                    ar >> d.m_epsn;
                    ar >> d.m_pg;
                    ar >> d.m_tv;
                    ar >> d.m_vn;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetVelocityGap(int nface, vec3d& vg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vg = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		vg += d.m_vg;
	}
    vg /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetPressureGap(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		pg += d.m_pg;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetViscousTraction(int nface, vec3d& tv)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    tv = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		tv += d.m_tv;
	}
    tv /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetNormalVelocity(int nface, double& vn)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vn = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		vn += d.m_vn;
	}
    vn /= ni;
}

//-----------------------------------------------------------------------------
// FETiedFluidInterface
//-----------------------------------------------------------------------------

FETiedFluidInterface::FETiedFluidInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_atol = 0.1;
    m_epst = 1;
    m_epsn = 1;
    m_btwo_pass = false;
    m_stol = 0.01;
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_ptol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedFluidInterface::~FETiedFluidInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedFluidInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // get the DOFS
    m_dofWX = fem.GetDOFIndex("wx");
    m_dofWY = fem.GetDOFIndex("wy");
    m_dofWZ = fem.GetDOFIndex("wz");
    m_dofEF = fem.GetDOFIndex("ef");
    
    vector<int> lm(4*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FETiedFluidSurface& ss = (np == 0? m_ss : m_ms);
        
        int ni = 0, k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k, ++ni)
            {
                FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(k));
                FESurfaceElement* pe = pt.m_pme;
                if (pe != 0)
                {
                    FESurfaceElement& me = *pe;
                    int* mn = &me.m_node[0];
                    
                    assign(lm, -1);
                    
                    int nseln = se.Nodes();
                    int nmeln = me.Nodes();
                    
                    for (l=0; l<nseln; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[4*l  ] = id[m_dofWX];
                        lm[4*l+1] = id[m_dofWY];
                        lm[4*l+2] = id[m_dofWZ];
                        lm[4*l+3] = id[m_dofEF];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[4*(l+nseln)  ] = id[m_dofWX];
                        lm[4*(l+nseln)+1] = id[m_dofWY];
                        lm[4*(l+nseln)+2] = id[m_dofWZ];
                        lm[4*(l+nseln)+3] = id[m_dofEF];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPressurePenalty(m_ss);
        if (m_btwo_pass) CalcAutoPressurePenalty(m_ms);
    }
    
    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    InitialProjection(m_ss, m_ms);
    if (m_btwo_pass) InitialProjection(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::CalcAutoPressurePenalty(FETiedFluidSurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPressurePenalty(el, s);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epst = eps;
            if (eps != 0) pt.m_epsn = 1./eps;
        }
    }
}

//-----------------------------------------------------------------------------
double FETiedFluidInterface::AutoPressurePenalty(FESurfaceElement& el, FETiedFluidSurface& s)
{
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // evaluate element surface normal at parametric center
    vec3d t[2];
    s.CoBaseVectors0(el, 0, 0, t);
    vec3d n = t[0] ^ t[1];
    n.unit();
    
    // get the element this surface element belongs to
    FEElement* pe = el.m_elem[0];
    if (pe == 0) return 0.0;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    
    // check that it is a fluid material
    FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
    if (fluid == 0) return 0.0;
    m_pfluid = fluid;

    // get a material point
    FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
    // get the shear viscosity
    double mu = fluid->GetViscous()->ShearViscosity(mp);

    // get the area of the surface element
    double A = s.FaceArea(el);
    
    // get the volume of the volume element
    double V = m.ElementVolume(*pe);
    
    return mu*A/V;
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidInterface::InitialProjection(FETiedFluidSurface& ss, FETiedFluidSurface& ms)
{
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    
    // loop over all integration points
    int n = 0;
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        int nint = el.GaussPoints();
        
        for (int j=0; j<nint; ++j, ++n)
        {
            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // calculate the normal at this integration point
            nu = ss.SurfaceNormal(el, j);
            
            // find the intersection point with the master surface
            pme = np.Project2(r, nu, rs);
            
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_pme = pme;
            pt.m_nu = nu;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                pt.m_Gap = q - r;
            }
            else
            {
                // the node is not in contact
                pt.m_Gap = vec3d(0,0,0);
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for fluid tangential velocity and fluid pressure
void FETiedFluidInterface::ProjectSurface(FETiedFluidSurface& ss, FETiedFluidSurface& ms)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    FESurfaceElement* pme;
    vec3d r;
    
    vec3d  vt[FEElement::MAX_NODES], v1;
    double ps[FEElement::MAX_NODES], p1;
    
    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        
        // get the nodal velocities and pressures
        for (int j=0; j<ne; ++j) {
            vt[j] = mesh.Node(el.m_node[j]).get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
            ps[j] = m_pfluid->Pressure(mesh.Node(el.m_node[j]).get(m_dofEF));
        }

        for (int j=0; j<nint; ++j)
        {
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // get the velocity and pressure at the integration point
            v1 = el.eval(vt, j);
            p1 = el.eval(ps, j);
            
            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // calculate the tangential velocity gap function
                vec3d vm[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) vm[k] = mesh.Node(pme->m_node[k]).get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
                vec3d v2 = pme->eval(vm, pt.m_rs[0], pt.m_rs[1]);
                pt.m_vg = v2 - v1;

                // calculate the pressure gap function
                double pm[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) pm[k] = m_pfluid->Pressure(mesh.Node(pme->m_node[k]).get(m_dofEF));
                double p2 = pme->eval(pm, pt.m_rs[0], pt.m_rs[1]);
                pt.m_pg = p1 - p2;
            }
            else
            {
                // the node is not tied
                pt.m_vg = vec3d(0,0,0);
                pt.m_pg = 0;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedFluidInterface::Update(int niter, const FETimeInfo& tp)
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_ss, m_ms);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss);
    
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    int i, j, k;
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    double N[8*MN];
    
    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get slave and master surface
        FETiedFluidSurface& ss = (np == 0? m_ss : m_ms);
        FETiedFluidSurface& ms = (np == 0? m_ms : m_ss);
        
        // loop over all slave elements
        for (i=0; i<ss.Elements(); ++i)
        {
            // get the surface element
            FESurfaceElement& se = ss.Element(i);
            
            // get the nr of nodes and integration points
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            
            // copy the LM vector; we'll need it later
            ss.UnpackLM(se, sLM);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (j=0; j<nint; ++j)
            {
                // get the base vectors
                vec3d g[2];
                ss.CoBaseVectors(se, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
                // integration weights
                w[j] = se.GaussWeights()[j];
            }
            
            // loop over all integration points
            // note that we are integrating over the current surface
            for (j=0; j<nint; ++j)
            {
				FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));

                // get the master element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the master element
                    FESurfaceElement& me = *pme;
                    
                    // get the nr of master element nodes
                    int nmeln = me.Nodes();
                    
                    // copy LM vector
                    ms.UnpackLM(me, mLM);
                    
                    // calculate degrees of freedom
                    int ndof = 3*(nseln + nmeln);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    for (k=0; k<nseln; ++k)
                    {
                        LM[3*k  ] = sLM[4*k  ];
                        LM[3*k+1] = sLM[4*k+1];
                        LM[3*k+2] = sLM[4*k+2];
                    }
                    
                    for (k=0; k<nmeln; ++k)
                    {
                        LM[3*(k+nseln)  ] = mLM[4*k  ];
                        LM[3*(k+nseln)+1] = mLM[4*k+1];
                        LM[3*(k+nseln)+2] = mLM[4*k+2];
                    }
                    
                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // get slave element shape functions
                    Hs = se.H(j);
                    
                    // get master element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // gap function
                    vec3d dg = pt.m_vg;
                    
                    // lagrange multiplier
                    vec3d Lm = pt.m_Lmd;
                    
                    // penalty
                    double eps = m_epst*pt.m_epst;
                    
                    // viscous traction
                    vec3d tv = Lm + dg*eps;
                    pt.m_tv = tv;
                    
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                    
                    for (k=0; k<nseln; ++k)
                    {
                        N[3*k  ] = Hs[k]*tv.x;
                        N[3*k+1] = Hs[k]*tv.y;
                        N[3*k+2] = Hs[k]*tv.z;
                    }
                    
                    for (k=0; k<nmeln; ++k)
                    {
                        N[3*(k+nseln)  ] = -Hm[k]*tv.x;
                        N[3*(k+nseln)+1] = -Hm[k]*tv.y;
                        N[3*(k+nseln)+2] = -Hm[k]*tv.z;
                    }
                    
                    for (k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
                    
                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                    
                    // do the pressure stuff
                    // calculate nr of pressure dofs
                    ndof = nseln + nmeln;
                    
                    // calculate the flow rate
                    double epsn = m_epsn*pt.m_epsn;
                    
                    double vn = pt.m_Lmp + epsn*pt.m_pg;
                    pt.m_vn = vn;
                    
                    // fill the LM
                    LM.resize(ndof);
                    for (k=0; k<nseln; ++k) LM[k        ] = sLM[4*k+3];
                    for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[4*k+3];
                    
                    // fill the force array
                    fe.resize(ndof);
                    zero(fe);
                    for (k=0; k<nseln; ++k) N[k      ] = -Hs[k];
                    for (k=0; k<nmeln; ++k) N[k+nseln] =  Hm[k];
                    
                    for (k=0; k<ndof; ++k) fe[k] += vn*N[k]*detJ[j]*w[j];

                    // assemble residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    int i, j, k, l;
    vector<int> sLM, mLM, LM, en;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN], pt[MN], dpr[MN], dps[MN];
    matrix ke;
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
        // get the slave and master surface
        FETiedFluidSurface& ss = (np == 0? m_ss : m_ms);
        FETiedFluidSurface& ms = (np == 0? m_ms : m_ss);
        
        // loop over all slave elements
        for (i=0; i<ss.Elements(); ++i)
        {
            // get ths slave element
            FESurfaceElement& se = ss.Element(i);
            
            // get nr of nodes and integration points
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            
            // nodal velocities and pressures
            vec3d vt[FEElement::MAX_NODES];
            double pn[FEElement::MAX_NODES];
            for (j=0; j<nseln; ++j) {
                vt[j] = ss.GetMesh()->Node(se.m_node[j]).get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
                pn[j] = m_pfluid->Pressure(ss.GetMesh()->Node(se.m_node[j]).get(m_dofEF));
            }
            
            // copy the LM vector
            ss.UnpackLM(se, sLM);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (j=0; j<nint; ++j)
            {
                // get the base vectors
                vec3d g[2];
                ss.CoBaseVectors(se, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
                // integration weights
                w[j] = se.GaussWeights()[j];
                
                // pressure
                pt[j] = se.eval(pn, j);
                dpr[j] = se.eval_deriv1(pn, j);
                dps[j] = se.eval_deriv2(pn, j);
            }
            
            // loop over all integration points
            for (j=0; j<nint; ++j)
            {
				FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));

                // get the master element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    FESurfaceElement& me = *pme;
                    
                    // get the nr of master nodes
                    int nmeln = me.Nodes();
                    
                    // nodal pressure
                    vec3d vm[FEElement::MAX_NODES];
                    double pm[FEElement::MAX_NODES];
                    for (k=0; k<nmeln; ++k) {
                        vm[k] = ms.GetMesh()->Node(me.m_node[k]).get_vec3d(m_dofWX, m_dofWY, m_dofWZ);
                        pm[k] = m_pfluid->Pressure(ms.GetMesh()->Node(me.m_node[k]).get(m_dofEF));
                    }
                    
                    // copy the LM vector
                    ms.UnpackLM(me, mLM);
                    
                    int ndpn;    // number of dofs per node
                    int ndof;    // number of dofs in stiffness matrix
                    
                    // calculate degrees of freedom for fluid-on-fluid contact
                    ndpn = 4;
                    ndof = ndpn*(nseln+nmeln);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    
                    for (k=0; k<nseln; ++k)
                    {
                        LM[4*k  ] = sLM[4*k  ];             // wx-dof
                        LM[4*k+1] = sLM[4*k+1];             // wy-dof
                        LM[4*k+2] = sLM[4*k+2];             // wz-dof
                        LM[4*k+3] = sLM[4*k+3];             // ef-dof
                    }
                    for (k=0; k<nmeln; ++k)
                    {
                        LM[4*(k+nseln)  ] = mLM[4*k  ];     // wx-dof
                        LM[4*(k+nseln)+1] = mLM[4*k+1];     // wy-dof
                        LM[4*(k+nseln)+2] = mLM[4*k+2];     // wz-dof
                        LM[4*(k+nseln)+3] = mLM[4*k+3];     // ef-dof
                    }

                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // slave shape functions
                    Hs = se.H(j);
                    
                    // master shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // penalty
                    double eps = m_epst*pt.m_epst;
                    
                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();
                    
                    // a. K-term
                    //------------------------------------
                    
                    for (k=0; k<nseln; ++k) {
                        for (l=0; l<nseln; ++l)
                        {
                            double K = eps*Hs[k]*Hs[l]*detJ[j]*w[j];
                            ke[ndpn*k    ][ndpn*l    ] += K;
                            ke[ndpn*k + 1][ndpn*l + 1] += K;
                            ke[ndpn*k + 2][ndpn*l + 2] += K;
                            
                        }
                        for (l=0; l<nmeln; ++l)
                        {
                            double K = -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
                            ke[ndpn*k    ][ndpn*(nseln+l)    ] += K;
                            ke[ndpn*k + 1][ndpn*(nseln+l) + 1] += K;
                            ke[ndpn*k + 2][ndpn*(nseln+l) + 2] += K;
                        }
                    }
                    
                    for (k=0; k<nmeln; ++k) {
                        for (l=0; l<nseln; ++l)
                        {
                            double K = -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
                            ke[ndpn*(nseln+k)    ][ndpn*l    ] += K;
                            ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += K;
                            ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += K;
                        }
                        for (l=0; l<nmeln; ++l)
                        {
                            double K = eps*Hm[k]*Hm[l]*detJ[j]*w[j];
                            ke[ndpn*(nseln+k)    ][ndpn*(nseln+l)    ] += K;
                            ke[ndpn*(nseln+k) + 1][ndpn*(nseln+l) + 1] += K;
                            ke[ndpn*(nseln+k) + 2][ndpn*(nseln+l) + 2] += K;
                        }
                    }
                    
                    // --- D I L A T A T I O N   S T I F F N E S S ---
                    {
                        double epsn = m_epsn*pt.m_epsn;
                        double K = m_pfluid->m_k;
                        
                        for (k=0; k<nseln; ++k) {
                            for (l=0; l<nseln; ++l)
                                ke[4*k + 3][4*l+3] += -K*epsn*w[j]*detJ[j]*Hs[k]*Hs[l];
                            for (l=0; l<nmeln; ++l)
                                ke[4*k + 3][4*(nseln+l)+3] += K*epsn*w[j]*detJ[j]*Hs[k]*Hm[l];
                        }
                        
                        for (k=0; k<nmeln; ++k) {
                            for (l=0; l<nseln; ++l)
                                ke[4*(nseln+k)+3][4*l + 3] += K*epsn*w[j]*detJ[j]*Hm[k]*Hs[l];
                            for (l=0; l<nmeln; ++l)
                                ke[4*(nseln+k)+3][4*(nseln+l) + 3] += -K*epsn*w[j]*detJ[j]*Hm[k]*Hm[l];
                        }
                        
                    }
                    
                    // assemble the global stiffness
                    psolver->AssembleStiffness(en, LM, ke);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (!m_blaugon) return true;
    
    int i;
    vec3d Ln;
    bool bconv = true;
    
    int NS = m_ss.Elements();
    int NM = m_ms.Elements();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normP0 = 0;
    for (int i=0; i<NS; ++i)
    {
		FESurfaceElement& se = m_ss.Element(i);
        for (int j=0; j<se.GaussPoints(); ++j)
        {
			FETiedFluidSurface::Data& ds = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
            normP0 += ds.m_Lmp*ds.m_Lmp;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& me = m_ms.Element(i);
        for (int j=0; j<me.GaussPoints(); ++j)
        {
			FETiedFluidSurface::Data& dm = static_cast<FETiedFluidSurface::Data&>(*me.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
            normP0 += dm.m_Lmp*dm.m_Lmp;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    double maxpg = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, normP1 = 0, eps, epsn;
    for (i=0; i<NS; ++i)
    {
		FESurfaceElement& se = m_ss.Element(i);
		for (int j = 0; j<se.GaussPoints(); ++j)
		{
			FETiedFluidSurface::Data& ds = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));

            if (ds.m_pme) {
                // update Lagrange multipliers on slave surface
                eps = m_epst*ds.m_epst;
                ds.m_Lmd = ds.m_Lmd + ds.m_vg*eps;
                maxgap = max(maxgap,sqrt(ds.m_vg*ds.m_vg));
                normL1 += ds.m_Lmd*ds.m_Lmd;
                
                epsn = m_epsn*ds.m_epsn;
                ds.m_Lmp = ds.m_Lmp + epsn*ds.m_pg;
                maxpg = max(maxpg,fabs(ds.m_pg));
                normP1 += ds.m_Lmp*ds.m_Lmp;
            }
        }
    }
    
    for (i=0; i<NM; ++i)
    {
		FESurfaceElement& me = m_ms.Element(i);
		for (int j = 0; j<me.GaussPoints(); ++j)
		{
			FETiedFluidSurface::Data& dm = static_cast<FETiedFluidSurface::Data&>(*me.GetMaterialPoint(j));

            if (dm.m_pme) {
                // update Lagrange multipliers on master surface
                eps = m_epst*dm.m_epst;
                dm.m_Lmd = dm.m_Lmd + dm.m_vg*eps;
                maxgap = max(maxgap,sqrt(dm.m_vg*dm.m_vg));
                normL1 += dm.m_Lmd*dm.m_Lmd;
                
                epsn = m_epsn*dm.m_epsn;
                dm.m_Lmp = dm.m_Lmp + epsn*dm.m_pg;
                maxpg = max(maxpg,fabs(dm.m_pg));
                normP1 += dm.m_Lmp*dm.m_Lmp;
            }
        }
    }
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    double pnorm = (normP1 != 0 ? fabs((normP1 - normP0) / normP1) : fabs(normP1 - normP0));
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    if ((m_ptol > 0) && (maxpg > m_ptol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    felog.printf(" tied fluid interface # %d\n", GetID());
    felog.printf("                        CURRENT        REQUIRED\n");
    felog.printf("    V multiplier : %15le", lnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
    felog.printf("    P multiplier        : %15le", pnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
    
    felog.printf("    maximum velocity gap  : %15le", maxgap);
    if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");
    felog.printf("    maximum pressure gap : %15le", maxpg);
    if (m_ptol > 0) felog.printf("%15le\n", m_ptol); else felog.printf("       ***\n");

    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
    
    // serialize pointers
    if (ar.IsShallow() == false)
    {
        if (ar.IsSaving())
        {
            int NE = m_ss.Elements();
            for (int i=0; i<NE; ++i)
            {
				FESurfaceElement& se = m_ss.Element(i);
                for (int j=0; j<se.GaussPoints(); ++j)
                {
					FETiedFluidSurface::Data& ds = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));
					FESurfaceElement* pe = ds.m_pme;
                    if (pe) ar << pe->m_lid; else ar << -1;
                }
            }
        }
        else
        {
			int lid = -1;
			int NE = m_ss.Elements();
			for (int i = 0; i<NE; ++i)
			{
				FESurfaceElement& se = m_ss.Element(i);
				for (int j = 0; j<se.GaussPoints(); ++j)
				{
					FETiedFluidSurface::Data& ds = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(j));
					ar >> lid;
                    ds.m_pme = (lid < 0 ? 0 : &m_ms.Element(lid));
                }
            }
        }
    }
}
