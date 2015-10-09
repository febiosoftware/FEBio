#include "stdafx.h"
#include "FEBioFluidPlot.h"
#include "FEFluidDomain.h"
#include "FEFluid.h"
#include "FEBioPlot/FEBioPlotFile.h"

//-----------------------------------------------------------------------------
bool FEPlotFluidDilatation::Save(FEDomain &dom, vector<float>& a)
{
    FEFluidDomain* pd = dynamic_cast<FEFluidDomain*>(&dom);
    if (pd)
    {
        int N = pd->Nodes();
        for (int i=0; i<N; ++i)
        {
            FENode& node = pd->Node(i);
            a.push_back((float) node.m_et);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotElasticFluidPressure::Save(FEDomain &dom, vector<float>& a)
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
			
			a.push_back((float) ew);
		}
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVolumeRatio::Save(FEDomain &dom, vector<float>& a)
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
            
            a.push_back((float) ew);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidDensity::Save(FEDomain &dom, vector<float>& a)
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
            
            a.push_back((float) ew);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVelocity::Save(FEDomain &dom, vector<float>& a)
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
            
            float f[3];
            f[0] = (float) r.x;
            f[1] = (float) r.y;
            f[2] = (float) r.z;
            
            a.push_back(f[0]);
            a.push_back(f[1]);
            a.push_back(f[2]);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidAcceleration::Save(FEDomain &dom, vector<float>& a)
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
            
            float f[3];
            f[0] = (float) r.x;
            f[1] = (float) r.y;
            f[2] = (float) r.z;
            
            a.push_back(f[0]);
            a.push_back(f[1]);
            a.push_back(f[2]);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVorticity::Save(FEDomain &dom, vector<float>& a)
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
            
            float f[3];
            f[0] = (float) r.x;
            f[1] = (float) r.y;
            f[2] = (float) r.z;
            
            a.push_back(f[0]);
            a.push_back(f[1]);
            a.push_back(f[2]);
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotElementFluidStress::Save(FEDomain& dom, vector<float>& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    // write solid element data
    int N = dom.Elements();
    for (int i=0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        float s[6] = {0};
        int nint = el.GaussPoints();
        double f = 1.0 / (double) nint;
        
        // since the PLOT file requires floats we need to convert
        // the doubles to single precision
        // we output the average stress values of the gauss points
        for (int j=0; j<nint; ++j)
        {
            FEFluidMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>());
            if (ppt)
            {
                FEFluidMaterialPoint& pt = *ppt;
                s[0] += (float) (f*pt.m_s.xx());
                s[1] += (float) (f*pt.m_s.yy());
                s[2] += (float) (f*pt.m_s.zz());
                s[3] += (float) (f*pt.m_s.xy());
                s[4] += (float) (f*pt.m_s.yz());
                s[5] += (float) (f*pt.m_s.xz());
            }
        }
        
        a.push_back(s[0]);
        a.push_back(s[1]);
        a.push_back(s[2]);
        a.push_back(s[3]);
        a.push_back(s[4]);
        a.push_back(s[5]);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotElementFluidRateOfDef::Save(FEDomain& dom, vector<float>& a)
{
    FEFluid* pme = dynamic_cast<FEFluid*>(dom.GetMaterial());
    if (pme == 0) return false;
    
    // write solid element data
    int N = dom.Elements();
    for (int i=0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        float s[6] = {0};
        int nint = el.GaussPoints();
        double f = 1.0 / (double) nint;
        
        // since the PLOT file requires floats we need to convert
        // the doubles to single precision
        // we output the average stress values of the gauss points
        for (int j=0; j<nint; ++j)
        {
            FEFluidMaterialPoint* ppt = (el.GetMaterialPoint(j)->ExtractData<FEFluidMaterialPoint>());
            if (ppt)
            {
                FEFluidMaterialPoint& pt = *ppt;
                mat3ds D = pt.RateOfDeformation();
                s[0] += (float) (f*D.xx());
                s[1] += (float) (f*D.yy());
                s[2] += (float) (f*D.zz());
                s[3] += (float) (f*D.xy());
                s[4] += (float) (f*D.yz());
                s[5] += (float) (f*D.xz());
            }
        }
        
        a.push_back(s[0]);
        a.push_back(s[1]);
        a.push_back(s[2]);
        a.push_back(s[3]);
        a.push_back(s[4]);
        a.push_back(s[5]);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidStressPower::Save(FEDomain &dom, vector<float>& a)
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
            
            a.push_back((float) ew);
        }
        return true;
    }
    return false;
}

