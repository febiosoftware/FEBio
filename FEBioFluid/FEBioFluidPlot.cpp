#include "stdafx.h"
#include "FEBioFluidPlot.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"
#include "FEBioPlot/FEBioPlotFile.h"

//-----------------------------------------------------------------------------
bool FEPlotFluidDilatation::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidDomain* pd = dynamic_cast<FEFluidDomain*>(&dom);
    if (pd)
    {
		const int dof_E = GetFEModel()->GetDOFIndex("e");
        int N = pd->Nodes();
        for (int i=0; i<N; ++i)
        {
            FENode& node = pd->Node(i);
            a << node.get(dof_E);
        }
        return true;
    }
    return false;
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
    if (strcmp(pcs->GetName(), "") == 0) return false;
    if (strcmp(pcs->GetName(), m_szdom) != 0) return false;
    
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
            m_elem[j] = pcs->FindElement(el);
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
        // get the element this surface element belongs to
        FEElement* pe = m_pMesh->FindElementFromID(m_elem[j]);
        if (pe)
        {
            // get the material
            FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluid* fluid = dynamic_cast<FEFluid*> (pm);
            if (fluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                mat3ds s(mat3dd(0));
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += pt.m_s;
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

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
bool FEPlotElasticFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
	if (dom.Class() != FE_DOMAIN_SOLID) return false;
	FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
	if (dynamic_cast<FEFluidDomain* >(&bd))
	{
		for (int i=0; i<bd.Elements(); ++i)
		{
			FESolidElement& el = bd.Element(i);
			
			// calculate average pressure
			double ew = 0;
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
				
				if (pt) ew += pt->m_p;
			}
			ew /= el.GaussPoints();
			
			a << ew;
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVolumeRatio::Save(FEDomain &dom, FEDataStream& a)
{
    if (dom.Class() != FE_DOMAIN_SOLID) return false;
    FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
    if (dynamic_cast<FEFluidDomain* >(&bd))
    {
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            
            // calculate average pressure
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
                
                if (pt) ew += pt->m_J;
            }
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidDensity::Save(FEDomain &dom, FEDataStream& a)
{
    if (dom.Class() != FE_DOMAIN_SOLID) return false;
    FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
    if (dynamic_cast<FEFluidDomain* >(&bd))
    {
        FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            
            // calculate average pressure
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
                
                if (pt) ew += pme->Density(mp);
            }
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        int BE = bd.Elements();
        for (int i=0; i<BE; ++i)
        {
            FESolidElement& el = bd.Element(i);
            int n = el.GaussPoints();
            vec3d r = vec3d(0,0,0);
            for (int j=0; j<n; ++j)
            {
                FEFluidMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>();
                r += pt.m_vt;
            }
            r /= n;
            
			a << r;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidAcceleration::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        int BE = bd.Elements();
        for (int i=0; i<BE; ++i)
        {
            FESolidElement& el = bd.Element(i);
            int n = el.GaussPoints();
            vec3d r = vec3d(0,0,0);
            for (int j=0; j<n; ++j)
            {
                FEFluidMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>();
                r += pt.m_at;
            }
            
            r /= n;
            
			a << r;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVorticity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
        int BE = bd.Elements();
        for (int i=0; i<BE; ++i)
        {
            FESolidElement& el = bd.Element(i);
            int n = el.GaussPoints();
            vec3d r = vec3d(0,0,0);
            for (int j=0; j<n; ++j)
            {
                FEFluidMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>();
                r += pt.Vorticity();
            }
            
            r /= n;

			a << r;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotElementFluidStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    // write solid element data
    int N = dom.Elements();
    for (int i=0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        int nint = el.GaussPoints();
        double f = 1.0 / (double) nint;
        
        // since the PLOT file requires floats we need to convert
        // the doubles to single precision
        // we output the average stress values of the gauss points
		mat3ds s(0.0);
        for (int j=0; j<nint; ++j)
        {
            FEFluidMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>());
            if (ppt)
            {
                FEFluidMaterialPoint& pt = *ppt;
				s += pt.m_s;
            }
        }
		s *= f;
        
		a << s;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotElementFluidRateOfDef::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    // write solid element data
    int N = dom.Elements();
    for (int i=0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        int nint = el.GaussPoints();
        double f = 1.0 / (double) nint;
        
        // since the PLOT file requires floats we need to convert
        // the doubles to single precision
        // we output the average stress values of the gauss points
		mat3ds s(0.0);
        for (int j=0; j<nint; ++j)
        {
            FEFluidMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>());
            if (ppt)
            {
                FEFluidMaterialPoint& pt = *ppt;
                mat3ds D = pt.RateOfDeformation();
				s += D;
            }
        }
		s *= f;
        
		a << s;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidStressPower::Save(FEDomain &dom, FEDataStream& a)
{
    if (dom.Class() != FE_DOMAIN_SOLID) return false;
    FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
    if (dynamic_cast<FEFluidDomain* >(&bd))
    {
        for (int i=0; i<bd.Elements(); ++i)
        {
            FESolidElement& el = bd.Element(i);
            
            // calculate average pressure
            double ew = 0;
            for (int j=0; j<el.GaussPoints(); ++j)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
                
                if (pt) ew += (pt->m_s*pt->m_L).trace();
            }
            ew /= el.GaussPoints();
            
            a << ew;
        }
        return true;
    }
    return false;
}
