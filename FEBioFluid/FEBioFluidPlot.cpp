#include "stdafx.h"
#include "FEBioFluidPlot.h"
#include "FEFluidDomain3D.h"
#include "FEFluidDomain2D.h"
#include "FEFluid.h"
#include "FEFluidDomain.h"
#include "FEFluidFSIDomain.h"
#include "FEFluidFSI.h"
#include "FEBioPlot/FEBioPlotFile.h"

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotDisplacement::Save(FEMesh& m, FEDataStream& a)
{
    // loop over all nodes
    for (int i=0; i<m.Nodes(); ++i)
    {
        FENode& node = m.Node(i);
        a << node.m_rt - node.m_r0;
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal fluid velocity
bool FEPlotNodalFluidVelocity::Save(FEMesh& m, FEDataStream& a)
{
    int dofWX = m_pfem->GetDOFIndex("wx");
    int dofWY = m_pfem->GetDOFIndex("wy");
    int dofWZ = m_pfem->GetDOFIndex("wz");
    int dofVX = m_pfem->GetDOFIndex("vx");
    int dofVY = m_pfem->GetDOFIndex("vy");
    int dofVZ = m_pfem->GetDOFIndex("vz");

	bool bvel = true;
	if ((dofVX == -1) || (dofVY == -1) || (dofVZ == -1))
	{
		bvel = false;
	}
    
    // loop over all nodes
    for (int i=0; i<m.Nodes(); ++i)
    {
        FENode& node = m.Node(i);
        vec3d vs = (bvel ? node.get_vec3d(dofVX, dofVY, dofVZ) : vec3d(0,0,0));
        vec3d w = node.get_vec3d(dofWX, dofWY, dofWZ);
        a << vs+w;
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal relative fluid velocity
bool FEPlotNodalRelativeFluidVelocity::Save(FEMesh& m, FEDataStream& a)
{
    int dofWX = m_pfem->GetDOFIndex("wx");
    int dofWY = m_pfem->GetDOFIndex("wy");
    int dofWZ = m_pfem->GetDOFIndex("wz");
    
    // loop over all nodes
    for (int i=0; i<m.Nodes(); ++i)
    {
        FENode& node = m.Node(i);
        vec3d vf = node.get_vec3d(dofWX, dofWY, dofWZ);
        a << vf;
    }
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal dilatations
bool FEPlotFluidDilatation::Save(FEMesh& m, FEDataStream& a)
{
    // get the dilatation dof index
    int dof_e = GetFEModel()->GetDOFIndex("ef");

    // loop over all nodes
    for (int i=0; i<m.Nodes(); ++i)
    {
        FENode& node = m.Node(i);
        a << node.get(dof_e);
    }
    return true;
}

//=============================================================================
//                       S U R F A C E    D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotFluidSurfaceForce::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    // Evaluate this field only for a specific domain, by checking domain name
    if (pcs->GetName() != m_szdom) return false;
    
    int NF = pcs->Elements();
    vec3d fn(0,0,0);    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        m_elem.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
            m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_elem[j];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
            FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm);
            if (fluid || fsi) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                mat3ds s(mat3dd(0));
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += pt.m_sf;
                }
                s /= nint;
                
                // Evaluate contribution to net force on surface.
                // Negate the fluid traction since we want the traction on the surface,
                // which is the opposite of the traction on the fluid.
                fn -= s*m_area[j];
            }
        }
    }
    
    // save results
	a << fn;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSurfaceTractionPower::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    // Evaluate this field only for a specific domain, by checking domain name
	if (pcs->GetName() != m_szdom) return false;

    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        m_elem.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
            m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_elem[j];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
            if (fluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                double s = 0;
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += pt.m_vft*(pt.m_sf*m_area[j]);
                }
                s /= nint;
                
                // Evaluate contribution to net traction power on surface.
                fn += s;
            }
        }
    }
    
    // save results
    a << fn;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSurfaceEnergyFlux::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    // Evaluate this field only for a specific domain, by checking domain name
	if (pcs->GetName() != m_szdom) return false;

    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        m_elem.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
            m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_elem[j];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
            if (fluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                double s = 0;
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += fluid->EnergyDensity(mp)*(pt.m_vft*m_area[j]);
                }
                s /= nint;
                
                // Evaluate contribution to net energy flux on surface.
                fn += s;
            }
        }
    }
    
    // save results
    a << fn;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidMassFlowRate::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    // Evaluate this field only for a specific domain, by checking domain name
	if (pcs->GetName() != m_szdom) return false;

    int dofWX = m_pfem->GetDOFIndex("wx");
    int dofWY = m_pfem->GetDOFIndex("wy");
    int dofWZ = m_pfem->GetDOFIndex("wz");
    int dofEF  = m_pfem->GetDOFIndex("ef");
    
    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // initialize on the first pass to identify solid element associated with this surface element
    if (m_binit) {
        m_elem.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
        }
        m_binit = false;
    }
    
    // calculate net fluid mass flow rate
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_elem[j];
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
            FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm);
            if (fluid || fsi) {
                double rhor;
                if (fluid) rhor = fluid->m_rhor;
                else rhor = fsi->Fluid()->m_rhor;
                // get the surface element
                FESurfaceElement& el = pcs->Element(j);
                // extract nodal velocities and dilatation
                int neln = el.Nodes();
                vec3d vt[FEElement::MAX_NODES];
                double et[FEElement::MAX_NODES];
                for (int j=0; j<neln; ++j) {
                    vt[j] = m_pMesh->Node(el.m_node[j]).get_vec3d(dofWX, dofWY, dofWZ);
                    et[j] = m_pMesh->Node(el.m_node[j]).get(dofEF);
                }
                
                // evaluate mass flux across this surface element
                int nint = el.GaussPoints();
                double*	gw = el.GaussWeights();
                vec3d gcov[2];
                for (int n=0; n<nint; ++n) {
                    vec3d v = el.eval(vt, n);
                    double J = 1 + el.eval(et, n);
                    pcs->CoBaseVectors(el, n, gcov);
                    fn += (v*(gcov[0] ^ gcov[1]))*rhor/J*gw[n];
                }
            }
        }
    }
    
    // save results
    a << fn;
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidFlowRate::Save(FESurface &surf, FEDataStream &a)
{
	FESurface* pcs = &surf;
	if (pcs == 0) return false;

	// Evaluate this field only for a specific domain, by checking domain name
	if (pcs->GetName() != m_szdom) return false;

	int NF = pcs->Elements();
	double fn = 0;    // initialize

	FEMesh* m_pMesh = pcs->GetMesh();

	// initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
	if (m_binit) {
		m_area.resize(NF);
		m_elem.resize(NF);
		for (int j = 0; j<NF; ++j)
		{
			FESurfaceElement& el = pcs->Element(j);
			m_area[j] = pcs->SurfaceNormal(el, 0, 0)*pcs->FaceArea(el);
			m_elem[j] = m_pMesh->FindElementFromID(pcs->FindElement(el));
		}
		m_binit = false;
	}

	// calculate net flow rate normal to this surface
	for (int j = 0; j<NF; ++j)
	{
		// get the element this surface element belongs to
		FEElement* pe = m_elem[j];
		if (pe)
		{
			// evaluate the average fluid flux in this element
			int nint = pe->GaussPoints();
			vec3d w(0, 0, 0);
			for (int n = 0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
				FEFluidMaterialPoint* ptf = mp.ExtractData<FEFluidMaterialPoint>();
				if (ptf) w += ptf->m_vft / ptf->m_Jf;
			}
			w /= nint;

			// Evaluate contribution to net flow rate across surface.
			fn += w*m_area[j];
		}
	}

	// save results
	a << fn;

	return true;
}

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------

class FEFluidPressure : public FEValuator<double>
{
public:
	double eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
		return (pt ? pt->m_pf : 0.0);
	}
};

bool FEPlotFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);

	if (dynamic_cast<FEFluidDomain* >(&bd) ||
		dynamic_cast<FEFluidFSIDomain* >(&bd))
	{
		writeAverageElementValue(dom, FEFluidPressure(), a);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// TODO: This plots the exact same value as FEFluidPressure. Delete?
bool FEPlotElasticFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
	writeAverageElementValue(dom, FEFluidPressure(), a);

    return true;
}

//-----------------------------------------------------------------------------
class FEFluidTemperature : public FEValuator<double>
{
public:
	FEFluidTemperature(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		return m_mat->Temperature(const_cast<FEMaterialPoint&>(mp));
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidTemperature::Save(FEDomain &dom, FEDataStream& a)
{
    FESolidDomain& bd = static_cast<FESolidDomain&>(dom);

    if (dynamic_cast<FEFluidDomain* >(&bd) ||
        dynamic_cast<FEFluidFSIDomain* >(&bd))
    {
		FEMaterial* pm = dom.GetMaterial();
		FEFluid* pfluid = dynamic_cast<FEFluid*>(pm);
		if (pfluid == nullptr)
		{
			FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
			if (pfsi) pfluid = pfsi->Fluid();
		}

		if (pfluid)
		{
			writeAverageElementValue(dom, FEFluidTemperature(pfluid), a);
			return true;
		}
    }
    return false;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidVolumeRatio : public FEValuator<double>
{
public:
	FEFluidVolumeRatio(FEModel* fem, FESolidDomain& dom) : m_dom(dom), m_el(0)
	{
		dofEF = fem->GetDOFIndex("ef");
	}

	double eval(const FEMaterialPoint& mp) override
	{
		if (m_el != mp.m_elem)
		{
			m_el = dynamic_cast<FESolidElement*>(mp.m_elem);
			FESolidElement& el = *m_el;
			FEMesh& mesh = *m_dom.GetMesh();
			int neln = el.Nodes();
			for (int j = 0; j<neln; ++j)
				et[j] = mesh.Node(el.m_node[j]).get(dofEF);
		}

		double  Jf = 1 + m_el->Evaluate(et, mp.m_index);
		return Jf;
	}

private:
	FESolidDomain&	m_dom;
	FESolidElement*	m_el;
	int dofEF;

	double et[FEElement::MAX_NODES];
};

bool FEPlotFluidVolumeRatio::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue(dom, FEFluidVolumeRatio(m_pfem, sd), a);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidDensity : public FEValuator<double>
{
public:
	FEFluidDensity(FEModel* fem, FESolidDomain& dom, FEFluid* pm) : m_dom(dom), m_mat(pm), m_el(0)
	{
		dofEF = fem->GetDOFIndex("ef");
	}

	double eval(const FEMaterialPoint& mp) override
	{
		if (m_el != mp.m_elem)
		{
			m_el = dynamic_cast<FESolidElement*>(mp.m_elem);
			FESolidElement& el = *m_el;
			FEMesh& mesh = *m_dom.GetMesh();
			int neln = el.Nodes();
			for (int j = 0; j<neln; ++j)
				et[j] = mesh.Node(el.m_node[j]).get(dofEF);
		}

		double rhor = m_mat->m_rhor;
		double Jf = 1 + m_el->Evaluate(et, mp.m_index);
		return rhor / Jf;
	}

private:
	FESolidDomain&	m_dom;
	FESolidElement*	m_el;
	FEFluid*	m_mat;
	int dofEF;

	double et[FEElement::MAX_NODES];
};

bool FEPlotFluidDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue(dom, FEFluidDensity(m_pfem, sd, pme), a);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidDensityRate : public FEValuator<double>
{
public:
	FEFluidDensityRate(FEModel* fem, FESolidDomain& dom, FEFluid* pm) : m_dom(dom), m_el(0), m_mat(pm)
	{
		dofVX = fem->GetDOFIndex("vx");
		dofVY = fem->GetDOFIndex("vy");
		dofVZ = fem->GetDOFIndex("vz");
		dofEF = fem->GetDOFIndex("ef");
		dofAEF = fem->GetDOFIndex("aef");
	}

	double eval(const FEMaterialPoint& mp)
	{
		if (m_el != mp.m_elem)
		{
			m_el = dynamic_cast<FESolidElement*>(mp.m_elem);

			FESolidElement& el = *m_el;
			FEMesh& mesh = *m_dom.GetMesh();

			int neln = m_el->Nodes();
			for (int j = 0; j<neln; ++j) {
				vt[j] = mesh.Node(el.m_node[j]).get_vec3d(dofVX, dofVY, dofVZ);
				et[j] = mesh.Node(el.m_node[j]).get(dofEF);
				aet[j] = mesh.Node(el.m_node[j]).get(dofAEF);
			}
		}

		FESolidElement& el = *m_el;
		double rhor = m_mat->m_rhor;
		double Jf = 1.0 + el.Evaluate(et, mp.m_index);
		double Jfdot = el.Evaluate(aet, mp.m_index);
		double divvs = m_dom.gradient(el, vt, mp.m_index).trace();

		return rhor / Jf*(divvs - Jfdot / Jf);
	}

private:
	FESolidDomain&	m_dom;
	FESolidElement*	m_el;
	FEFluid* m_mat;

	int dofVX, dofVY, dofVZ;
	int dofEF, dofAEF;
	vec3d vt[FEElement::MAX_NODES];
	double et[FEElement::MAX_NODES];
	double aet[FEElement::MAX_NODES];
};

bool FEPlotFluidDensityRate::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue(sd, FEFluidDensityRate(m_pfem, sd, pme), a);
		return true;
    }
    
    return false;
}

//-----------------------------------------------------------------------------
class FEFluidVelocity : public FEValuator<vec3d>
{
public:
	vec3d eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_vft : vec3d(0.));
	}
};

bool FEPlotFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEFluidVelocity(), a);

	return true;
}

//-----------------------------------------------------------------------------
class FERelativeFluidVelocity : public FEValuator<vec3d>
{
public:
	vec3d eval(const FEMaterialPoint& mp) override
	{
		const FEFSIMaterialPoint* ppt = mp.ExtractData<FEFSIMaterialPoint>();
		return (ppt ? ppt->m_w : vec3d(0.));
	}
};

bool FEPlotRelativeFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
	writeAverageElementValue(dom, FERelativeFluidVelocity(), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidAcceleration : public FEValuator<vec3d>
{
public:
	vec3d eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_aft : vec3d(0.));
	}
};

bool FEPlotFluidAcceleration::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEFluidAcceleration(), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidVorticity : public FEValuator<vec3d>
{
public:
	vec3d eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->Vorticity() : vec3d(0.));
	}
};

bool FEPlotFluidVorticity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEFluidVorticity(), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEElementFluidStress : public FEValuator<mat3ds>
{
public:
	mat3ds eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_sf : mat3ds(0.));
	}
};

//! Store the average stresses for each element.
bool FEPlotElementFluidStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEElementFluidStress(), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEElementFluidRateOfDef : public FEValuator<mat3ds>
{
public:
	mat3ds eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->RateOfDeformation() : mat3ds(0.));
	}
};

//! Store the average stresses for each element.
bool FEPlotElementFluidRateOfDef::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEElementFluidRateOfDef(), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidStressPowerDensity : public FEValuator<double>
{
public:
	double eval(const FEMaterialPoint& mp) override
	{
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? (ppt->m_sf*ppt->m_Lf).trace() : 0.0);
	}
};

bool FEPlotFluidStressPowerDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEFluidStressPowerDensity(), a);

    return true;
}

//-----------------------------------------------------------------------------
class FEFluidHeatSupplyDensity : public FEValuator<double>
{
public:
	FEFluidHeatSupplyDensity(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		FEFluidMaterialPoint* ppt = (mp_noconst.ExtractData<FEFluidMaterialPoint>());
		return (ppt ? -(m_mat->GetViscous()->Stress(mp_noconst)*ppt->m_Lf).trace() : 0.0);
	}
private:
	FEFluid*	m_mat;
};

bool FEPlotFluidHeatSupplyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;
    
    // write solid element data
	writeAverageElementValue(dom, FEFluidHeatSupplyDensity(pme), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidShearViscosity : public FEValuator<double>
{
public:
	FEFluidShearViscosity(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->GetViscous()->ShearViscosity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidShearViscosity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    // write solid element data
	writeAverageElementValue(dom, FEFluidShearViscosity(pme), a);

	return true;
}

//-----------------------------------------------------------------------------
class FEFluidStrainEnergyDensity : public FEValuator<double>
{
public:
	FEFluidStrainEnergyDensity(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->StrainEnergyDensity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    // write solid element data
	writeAverageElementValue(dom, FEFluidStrainEnergyDensity(pme), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidKineticEnergyDensity : public FEValuator<double>
{
public:
	FEFluidKineticEnergyDensity(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->KineticEnergyDensity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidKineticEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    // write solid element data
	writeAverageElementValue(dom, FEFluidKineticEnergyDensity(pme), a);
    
    return true;
}

//-----------------------------------------------------------------------------
class FEFluidEnergyDensity : public FEValuator<double>
{
public:
	FEFluidEnergyDensity(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->EnergyDensity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    // write solid element data
	writeAverageElementValue(dom, FEFluidEnergyDensity(pme), a);

    return true;
}

//-----------------------------------------------------------------------------
class FEFluidElementStrainEnergy : public FEValuator<double>
{
public:
	FEFluidElementStrainEnergy(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->StrainEnergyDensity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidElementStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue(sd, FEFluidElementStrainEnergy(pme), a);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
class FEFluidElementKineticEnergy : public FEValuator<double>
{
public:
	FEFluidElementKineticEnergy(FEFluid* mat) : m_mat(mat) {}
	double eval(const FEMaterialPoint& mp) override
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return m_mat->KineticEnergyDensity(mp_noconst);
	}

private:
	FEFluid*	m_mat;
};

bool FEPlotFluidElementKineticEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue(sd, FEFluidElementKineticEnergy(pme), a);
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidElementCenterOfMass::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	double dens = pme->m_rhor;
        
	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			double* gw = el.GaussWeights();
                
			// integrate zeroth and first mass moments
			vec3d ew = vec3d(0,0,0);
			double m = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEFluidMaterialPoint& pt = *(el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>());
				double detJ = bd.detJ0(el, j)*gw[j];
				ew += pt.m_r0*(dens*detJ);
				m += dens*detJ;
			}
                
			a << ew/m;
		}
		return true;
    }
	return false;
}

//-----------------------------------------------------------------------------
class FEFluidElementLinearMomentum : public FEValuator<vec3d>
{
public:
	FEFluidElementLinearMomentum(FEFluid* pm) : m_mat(pm) {}
	vec3d eval(const FEMaterialPoint& mp) override
	{
		double dens = m_mat->m_rhor;
		const FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
		return pt.m_vft*dens;
	}

private:
	FEFluid* m_mat;
};

bool FEPlotFluidElementLinearMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue(bd, FEFluidElementLinearMomentum(pme), a);
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
class FEFluidElementAngularMomentum : public FEValuator<vec3d>
{
public:
	FEFluidElementAngularMomentum(FEFluid* pm) : m_mat(pm) {}
	vec3d eval(const FEMaterialPoint& mp) override
	{
		double dens = m_mat->m_rhor;
		const FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
		return pt.m_vft*dens;
	}

private:
	FEFluid* m_mat;
};

bool FEPlotFluidElementAngularMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    FEFluidFSI* sme = dynamic_cast<FEFluidFSI*>(dom.GetMaterial());
    if ((pme == 0) && (sme == 0)) return false;
	if (pme == 0) pme = sme->Fluid();
	if (pme == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue(bd, FEFluidElementAngularMomentum(pme), a);
		return true;
    }
    return false;
}
