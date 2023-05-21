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
#include "FEBioFluidPlot.h"
#include "FEFluidDomain3D.h"
#include "FEFluidMaterial.h"
#include "FEFluid.h"
#include "FEPolarFluid.h"
#include "FEFluidDomain.h"
#include "FEFluidFSIDomain.h"
#include "FEFluidFSI.h"
#include "FEBiphasicFSIDomain.h"
#include "FEBiphasicFSI.h"
#include "FEMultiphasicFSIDomain.h"
#include "FEMultiphasicFSI.h"
#include "FEThermoFluid.h"
#include <FECore/FEModel.h>
#include <FECore/FESurface.h>
#include <FECore/writeplot.h>
#include <FECore/FEDomainParameter.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Store the nodal displacements
bool FEPlotDisplacement::Save(FEMesh& m, FEDataStream& a)
{
    // loop over all nodes
	writeNodalValues<vec3d>(m, a, [](const FENode& node) {
		return node.m_rt - node.m_r0;
	});
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal fluid velocity (only for CFD domains)
bool FEPlotNodalFluidVelocity::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
    int dofWX = fem->GetDOFIndex("wx");
    int dofWY = fem->GetDOFIndex("wy");
    int dofWZ = fem->GetDOFIndex("wz");
    int dofVX = fem->GetDOFIndex("vfx");
    int dofVY = fem->GetDOFIndex("vfy");
    int dofVZ = fem->GetDOFIndex("vfz");

	bool bcfd = false;
	if ((dofVX == -1) && (dofVY == -1) && (dofVZ == -1))
	{
		bcfd = true;
	}
    
    if (bcfd) {
        writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
            return node.get_vec3d(dofWX, dofWY, dofWZ);
        });
        return true;
    }
    else {
        writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
            return node.get_vec3d(dofVX, dofVY, dofVZ);
        });
        return true;
    }
}

//-----------------------------------------------------------------------------
//! Store the nodal relative fluid velocity
bool FEPlotNodalRelativeFluidVelocity::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	int dofWX = fem->GetDOFIndex("wx");
    int dofWY = fem->GetDOFIndex("wy");
    int dofWZ = fem->GetDOFIndex("wz");
    
	writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
		return node.get_vec3d(dofWX, dofWY, dofWZ);
	});
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal dilatations
bool FEPlotFluidDilatation::Save(FEMesh& m, FEDataStream& a)
{
    // get the dilatation dof index
    int dof_e = GetFEModel()->GetDOFIndex("ef");
	if (dof_e < 0) return false;

    // loop over all nodes
	writeNodalValues<double>(m, a, [=](const FENode& node) {
		return node.get(dof_e);
	});
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal effective fluid pressure
bool FEPlotFluidEffectivePressure::Save(FEDomain& dom, FEDataStream& a)
{
    // get the dilatation dof index
    int dof_e = GetFEModel()->GetDOFIndex("ef");
    if (dof_e < 0) return false;
    
    FEFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEFluid>();
    if (pfluid == 0) return false;
    
    // loop over all nodes
    writeNodalValues<double>(dom, a, [=, &dom](int i) {
        FENode& node = dom.Node(i);
        return pfluid->Pressure(node.get(dof_e));
    });
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal polar fluid angular velocity
bool FEPlotNodalPolarFluidAngularVelocity::Save(FEMesh& m, FEDataStream& a)
{
    FEModel* fem = GetFEModel();
    int dofGX = fem->GetDOFIndex("gx");
    int dofGY = fem->GetDOFIndex("gy");
    int dofGZ = fem->GetDOFIndex("gz");
    
    writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
        return node.get_vec3d(dofGX, dofGY, dofGZ);
    });
    return true;
}

//-----------------------------------------------------------------------------
//! Store the nodal temperatures
bool FEPlotNodalFluidTemperature::Save(FEMesh& m, FEDataStream& a)
{
    // get the dilatation dof index
    int dof_T = GetFEModel()->GetDOFIndex("T");
    if (dof_T < 0) return false;

    // loop over all nodes
    writeNodalValues<double>(m, a, [=](const FENode& node) {
        return node.get(dof_T);
    });
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
    
    int NF = pcs->Elements();
    vec3d fn(0,0,0);    // initialize
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
		FESurfaceElement& el = pcs->Element(j);

        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* pfluid = pm->ExtractProperty<FEFluidMaterial>();

            if (!pfluid) {
                pe = el.m_elem[1];
                if (pe) pfluid = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEFluidMaterial>();
            }

            // see if this is a fluid element
            if (pfluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                mat3d s(mat3dd(0));
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += pt.m_sf;
                    if (pfluid->GetViscousPolar())
                        s += pfluid->GetViscousPolar()->SkewStress(mp);
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
bool FEPlotFluidSurfaceMoment::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    int NF = pcs->Elements();
    vec3d mn(0,0,0);    // initialize
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
        }
        m_binit = false;
    }
    
    // calculate net fluid moment
    for (int j=0; j<NF; ++j)
    {
        FESurfaceElement& el = pcs->Element(j);
        
        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* pfluid = pm->ExtractProperty<FEFluidMaterial>();
            
            if (!pfluid) {
                pe = el.m_elem[1];
                if (pe) pfluid = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEFluidMaterial>();
            }
            
            // see if this is a fluid element
            if (pfluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                mat3d s(mat3dd(0));
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    if (pfluid->GetViscousPolar())
                        s += pfluid->GetViscousPolar()->CoupleStress(mp);
                }
                s /= nint;
                
                // Evaluate contribution to net moment on surface.
                // Negate the fluid couple vector since we want the couple vector on the surface,
                // which is the opposite of the traction on the fluid.
                mn -= s*m_area[j];
            }
        }
    }
    
    // save results
    a << mn;
    
    return true;
}

//-----------------------------------------------------------------------------
// Plot contact pressure
bool FEPlotFluidSurfacePressure::Save(FESurface &surf, FEDataStream& a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;
    
    const int dof_EF = GetFEModel()->GetDOFIndex("ef");
    const int dof_T = GetFEModel()->GetDOFIndex("T");
    
    writeElementValue<double>(surf, a, [=](int nface) {
        FESurfaceElement& el = pcs->Element(nface);
        double ef = pcs->Evaluate(nface, dof_EF);
        double T = pcs->Evaluate(nface, dof_T);
        FEElement* pe = el.m_elem[0];
        if (pe) {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
            if (!fluid) {
                pe = el.m_elem[1];
                if (pe) fluid = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEFluidMaterial>();
            }
            if (fluid) return fluid->Pressure(ef, T);
            else return 0.;
        }
        else return 0.;
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSurfaceTractionPower::Save(FESurface &surf, FEDataStream &a)
{
    FESurface* pcs = &surf;
    if (pcs == 0) return false;

    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
		FESurfaceElement& el = pcs->Element(j);

        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* pfluid = pm->ExtractProperty<FEFluidMaterial>();

            // see if this is a fluid element
            if (pfluid) {
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

    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    // initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
    if (m_binit) {
        m_area.resize(NF);
        for (int j=0; j<NF; ++j)
        {
            FESurfaceElement& el = pcs->Element(j);
            m_area[j] = pcs->SurfaceNormal(el,0,0)*pcs->FaceArea(el);
        }
        m_binit = false;
    }
    
    // calculate net fluid force
    for (int j=0; j<NF; ++j)
    {
		FESurfaceElement& el = pcs->Element(j);

        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* pfluid = pm->ExtractProperty<FEFluidMaterial>();

            // see if this is a fluid element
            if (pfluid) {
                // evaluate the average stress in this element
                int nint = pe->GaussPoints();
                double s = 0;
                for (int n=0; n<nint; ++n)
                {
                    FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                    FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
                    s += pfluid->EnergyDensity(mp)*(pt.m_vft*m_area[j]);
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

	FEModel* fem = GetFEModel();
    int dofWX = fem->GetDOFIndex("wx");
    int dofWY = fem->GetDOFIndex("wy");
    int dofWZ = fem->GetDOFIndex("wz");
    int dofEF = fem->GetDOFIndex("ef");
    
    int NF = pcs->Elements();
    double fn = 0;    // initialize
    
    FEMesh* m_pMesh = pcs->GetMesh();
    
    // calculate net fluid mass flow rate
    for (int j=0; j<NF; ++j)
    {
		FESurfaceElement& el = pcs->Element(j);

        // get the element this surface element belongs to
        FEElement* pe = el.m_elem[0];
        if (pe)
        {
            // get the material
            FEMaterial* pm = fem->GetMaterial(pe->GetMatID());
            
            // see if this is a fluid element
            FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
            if (fluid) {
                double rhor = fluid->m_rhor;
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

	int NF = pcs->Elements();
	double fn = 0;    // initialize

	// initialize on the first pass to calculate the vectorial area of each surface element and to identify solid element associated with this surface element
	if (m_binit) {
		m_area.resize(NF);
		for (int j = 0; j<NF; ++j)
		{
			FESurfaceElement& el = pcs->Element(j);
			m_area[j] = pcs->SurfaceNormal(el, 0, 0)*pcs->FaceArea(el);
		}
		m_binit = false;
	}

	// calculate net flow rate normal to this surface
	for (int j = 0; j<NF; ++j)
	{
		FESurfaceElement& el = pcs->Element(j);

		// get the element this surface element belongs to
		FEElement* pe = el.m_elem[0];
		if (pe)
		{
			// evaluate the average fluid flux in this element
			int nint = pe->GaussPoints();
			vec3d w(0, 0, 0);
			for (int n = 0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
				FEFluidMaterialPoint* ptf = mp.ExtractData<FEFluidMaterialPoint>();
				if (ptf) w += ptf->m_vft / (ptf->m_ef + 1);
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
bool FEPlotFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
		return (pt ? pt->m_pf : 0.0);
	});
	return true;
}

//-----------------------------------------------------------------------------
// TODO: This plots the exact same value as FEFluidPressure. Delete?
bool FEPlotElasticFluidPressure::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
	writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* pt = (mp.ExtractData<FEFluidMaterialPoint>());
		return (pt ? pt->m_pf : 0.0);
	});
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidTemperature::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		return pfluid->Temperature(const_cast<FEMaterialPoint&>(mp));
	});
	return true;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidVolumeRatio
{
public:
	FEFluidVolumeRatio(FEModel* fem, FESolidDomain& dom) : m_dom(dom), m_el(0)
	{
		dofEF = fem->GetDOFIndex("ef");
	}

	double operator()(const FEMaterialPoint& mp)
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
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue<double>(dom, a, FEFluidVolumeRatio(GetFEModel(), sd));
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidDensity
{
public:
	FEFluidDensity(FEModel* fem, FESolidDomain& dom, FEFluidMaterial* pm) : m_dom(dom), m_mat(pm), m_el(0)
	{
		dofEF = fem->GetDOFIndex("ef");
	}

	double operator()(const FEMaterialPoint& mp)
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
	FEFluidMaterial*	m_mat;
	int dofEF;

	double et[FEElement::MAX_NODES];
};

bool FEPlotFluidDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue<double>(dom, a, FEFluidDensity(GetFEModel(), sd, pfluid));
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// NOTE: This is not thread safe!
class FEFluidDensityRate
{
public:
	FEFluidDensityRate(FEModel* fem, FESolidDomain& dom, FEFluidMaterial* pm) : m_dom(dom), m_el(0), m_mat(pm)
	{
		dofVX = fem->GetDOFIndex("vx");
		dofVY = fem->GetDOFIndex("vy");
		dofVZ = fem->GetDOFIndex("vz");
		dofEF = fem->GetDOFIndex("ef");
		dofAEF = fem->GetDOFIndex("aef");
	}

	double operator()(const FEMaterialPoint& mp)
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
	FEFluidMaterial* m_mat;

	int dofVX, dofVY, dofVZ;
	int dofEF, dofAEF;
	vec3d vt[FEElement::MAX_NODES];
	double et[FEElement::MAX_NODES];
	double aet[FEElement::MAX_NODES];
};

bool FEPlotFluidDensityRate::Save(FEDomain &dom, FEDataStream& a)
{
	FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeAverageElementValue<double>(sd, a, FEFluidDensityRate(GetFEModel(), sd, pfluid));
		return true;
    }
    
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_vft : vec3d(0.));
	});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotRelativeFluidVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* fpt = mp.ExtractData<FEFluidMaterialPoint>();
        const FEElasticMaterialPoint* ept = mp.ExtractData<FEElasticMaterialPoint>();
        return (fpt ? fpt->m_vft - ept->m_v : vec3d(0.0));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotBFSIPorosity::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicFSI* bp = dom.GetMaterial()->ExtractProperty<FEBiphasicFSI>();
    if (bp == 0) return false;

    // write solid element data
    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        return bp->Porosity(const_cast<FEMaterialPoint&>(mp));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotBFSISolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicFSI* bp = dom.GetMaterial()->ExtractProperty<FEBiphasicFSI>();
    if (bp == 0) return false;
    
    // write solid element data
    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        return bp->SolidVolumeFrac(const_cast<FEMaterialPoint&>(mp));
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFSIFluidFlux::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFSIMaterialPoint* ppt = mp.ExtractData<FEFSIMaterialPoint>();
        return (ppt ? ppt->m_w : vec3d(0.0));
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPermeability::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicFSI* bp = dom.GetMaterial()->ExtractProperty<FEBiphasicFSI>();
    if (bp == 0) return false;
    
    writeAverageElementValue<mat3ds>(dom, a, [=](const FEMaterialPoint& mp) {
        return bp->Permeability(const_cast<FEMaterialPoint&>(mp));
    });
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotGradJ::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicFSI* pbfsi = dom.GetMaterial()->ExtractProperty<FEBiphasicFSI>();
    if (pbfsi == 0) return false;
    
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEBiphasicFSIMaterialPoint* bpt = mp.ExtractData<FEBiphasicFSIMaterialPoint>();
        return (bpt ? bpt->m_gradJ : vec3d(0.0));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotGradPhiF::Save(FEDomain &dom, FEDataStream& a)
{
    FEBiphasicFSI* bp = dom.GetMaterial()->ExtractProperty<FEBiphasicFSI>();
    if (bp == 0) return false;
    
    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [=](const FEMaterialPoint& mp) {
        return bp->gradPorosity(const_cast<FEMaterialPoint&>(mp));
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidAcceleration::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_aft : vec3d(0.));
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidVorticity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->Vorticity() : vec3d(0.));
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPolarFluidAngularVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEPolarFluidMaterialPoint* ppt = mp.ExtractData<FEPolarFluidMaterialPoint>();
        const FEFluidMaterialPoint* pt = mp.ExtractData<FEFluidMaterialPoint>();
        vec3d g(0,0,0);
        if (ppt) g = ppt->m_gf;
        else if (pt) g = pt->Vorticity()/2;
        return g;
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPolarFluidRelativeAngularVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* pt = mp.ExtractData<FEFluidMaterialPoint>();
        const FEPolarFluidMaterialPoint* ppt = mp.ExtractData<FEPolarFluidMaterialPoint>();
        return (ppt ? ppt->m_gf - pt->Vorticity()/2 : vec3d(0.));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotPolarFluidRegionalAngularVelocity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
        return (ppt ? ppt->Vorticity()/2 : vec3d(0.));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidHeatFlux::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;

    // write solid element data
    writeAverageElementValue<vec3d>(dom, a, [](const FEMaterialPoint& mp) {
        const FEThermoFluidMaterialPoint* ppt = mp.ExtractData<FEThermoFluidMaterialPoint>();
        return (ppt ? ppt->m_q : vec3d(0.));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotFluidStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->m_sf : mat3ds(0.));
	});
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotElementFluidRateOfDef::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<mat3ds>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? ppt->RateOfDeformation() : mat3ds(0.));
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidStressPowerDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
		const FEFluidMaterialPoint* ppt = mp.ExtractData<FEFluidMaterialPoint>();
		return (ppt ? (ppt->m_sf*ppt->m_Lf).trace() : 0.0);
	});

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidHeatSupplyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		FEFluidMaterialPoint* ppt = (mp_noconst.ExtractData<FEFluidMaterialPoint>());
		return (ppt ? -(pfluid->GetViscous()->Stress(mp_noconst)*ppt->m_Lf).trace() : 0.0);
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidShearViscosity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->GetViscous()->ShearViscosity(mp_noconst);
	});

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidStrainEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->StrainEnergyDensity(mp_noconst);
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidKineticEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->KineticEnergyDensity(mp_noconst);
	});
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidEnergyDensity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    // write solid element data
	writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		return pfluid->EnergyDensity(mp_noconst);
	});

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidBulkModulus::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;

    // write solid element data
    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->BulkModulus(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidElementStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(sd, a, [=](const FEMaterialPoint& mp) {
			FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
			return pfluid->StrainEnergyDensity(mp_noconst);
		});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidElementKineticEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<double>(sd, a, [=](const FEMaterialPoint& mp) {
			FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
			return pfluid->KineticEnergyDensity(mp_noconst);
		});
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// NOTE: Can't use helper functions due to division by m.
bool FEPlotFluidElementCenterOfMass::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	double dens = pfluid->m_rhor;
        
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
bool FEPlotFluidElementLinearMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<vec3d>(bd, a, [=](const FEMaterialPoint& mp) {
			double dens = pfluid->m_rhor;
			const FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
			return pt.m_vft*dens;
		});
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidElementAngularMomentum::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		FESolidDomain& bd = static_cast<FESolidDomain&>(dom);
		writeIntegratedElementValue<vec3d>(bd, a, [=](const FEMaterialPoint& mp) {
			double dens = pfluid->m_rhor;
			const FEFluidMaterialPoint& pt = *(mp.ExtractData<FEFluidMaterialPoint>());
			return pt.m_vft*dens;
		});
		return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificFreeEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificFreeEnergy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificEntropy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificEntropy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificInternalEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificInternalEnergy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificGageEnthalpy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificGageEnthalpy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificFreeEnthalpy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificFreeEnthalpy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidSpecificStrainEnergy::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->SpecificStrainEnergy(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidIsochoricSpecificHeatCapacity::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->IsochoricSpecificHeatCapacity(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidIsobaricSpecificHeatCapacity::Save(FEDomain &dom, FEDataStream& a)
{
    FEElasticFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEElasticFluid>();
    if (pfluid == 0) return false;

    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->IsobaricSpecificHeatCapacity(mp_noconst);
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidThermalConductivity::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidThermalConductivity* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidThermalConductivity>();
    if (pfluid == 0) return false;
    
    writeAverageElementValue<double>(dom, a, [=](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->ThermalConductivity(mp_noconst);
    });

    return true;
}

//-----------------------------------------------------------------------------
//! Store the average stresses for each element.
bool FEPlotFSISolidStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEFluid>();
    if (pfluid == 0) return false;
    
    if (dom.Class() == FE_DOMAIN_SOLID)
    {
        FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
        writeAverageElementValue<mat3ds>(sd, a, [](const FEMaterialPoint& mp) {
            const FEFSIMaterialPoint* pt = (mp.ExtractData<FEFSIMaterialPoint>());
            return (pt ? pt->m_ss : mat3ds(0.0));
        });
        return true;
    }
    
    return false;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidShearStressError::Save(FEDomain& dom, FEDataStream& a)
{
	FEFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEFluid>();
	if (pfluid == 0) return false;

	if (dom.Class() == FE_DOMAIN_SOLID)
	{
		writeRelativeError(dom, a, [](FEMaterialPoint& mp) {
			FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
			mat3ds s = fp.m_sf;
			double v = s.max_shear();
			return v;
		});
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! Store the average polar fluid stresses for each element.
bool FEPlotPolarFluidStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEViscousPolarFluid* vpfluid = dom.GetMaterial()->ExtractProperty<FEViscousPolarFluid>();
    if (vpfluid == 0) return false;
    FEViscousFluid* vfluid = dom.GetMaterial()->ExtractProperty<FEViscousFluid>();

    // write solid element data
    writeAverageElementValue<mat3d>(dom, a, [&](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return (vpfluid->SkewStress(mp_noconst) + vfluid->Stress(mp_noconst));
    });
    
    return true;
}

//-----------------------------------------------------------------------------
//! Store the average polar fluid couple stresses for each element.
bool FEPlotPolarFluidCoupleStress::Save(FEDomain& dom, FEDataStream& a)
{
    FEViscousPolarFluid* pfluid = dom.GetMaterial()->ExtractProperty<FEViscousPolarFluid>();
    if (pfluid == 0) return false;

    // write solid element data
    writeAverageElementValue<mat3d>(dom, a, [&](const FEMaterialPoint& mp) {
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        return pfluid->CoupleStress(mp_noconst);
    });
    
    return true;
}

//-----------------------------------------------------------------------------
bool FEPlotFluidRelativeReynoldsNumber::Save(FEDomain &dom, FEDataStream& a)
{
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    writeAverageElementValue<double>(dom, a, [&pfluid](const FEMaterialPoint& mp) {
        const FEFluidMaterialPoint* fpt = mp.ExtractData<FEFluidMaterialPoint>();
        const FEElasticMaterialPoint* ept = mp.ExtractData<FEElasticMaterialPoint>();
        FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
        double nu = pfluid->KinematicViscosity(mp_noconst);
        vec3d v(0,0,0);
        if (ept) v = ept->m_v;
        return (fpt->m_vft - v).Length()/nu;
    });
    
    return true;
}

//=================================================================================================
//-----------------------------------------------------------------------------
FEPlotFluidRelativePecletNumber::FEPlotFluidRelativePecletNumber(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
    DOFS& dofs = pfem->GetDOFS();
    int nsol = dofs.GetVariableSize("concentration");
    SetArraySize(nsol);
    
    // collect the names
    int ndata = pfem->GlobalDataItems();
    vector<string> s;
    for (int i = 0; i<ndata; ++i)
    {
        FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
        if (ps)
        {
            s.push_back(ps->GetName());
            m_sol.push_back(ps->GetID());
        }
    }
    assert(nsol == (int)s.size());
    SetArrayNames(s);
    SetUnits(UNIT_RECIPROCAL_LENGTH);
}

//-----------------------------------------------------------------------------
bool FEPlotFluidRelativePecletNumber::Save(FEDomain &dom, FEDataStream& a)
{
    FESoluteInterface* pm = dynamic_cast<FESoluteInterface*>(dom.GetMaterial());
    if (pm == 0) return false;
    
    FEFluidMaterial* pfluid = dom.GetMaterial()->ExtractProperty<FEFluidMaterial>();
    if (pfluid == 0) return false;
    
    // figure out the local solute IDs. This depends on the material
    int nsols = (int)m_sol.size();
    vector<int> lid(nsols, -1);
    int negs = 0;
    for (int i = 0; i<(int)m_sol.size(); ++i)
    {
        lid[i] = pm->FindLocalSoluteID(m_sol[i]);
        if (lid[i] < 0) negs++;
    }
    if (negs == nsols) return false;
    
    // loop over all elements
    int N = dom.Elements();
    for (int i = 0; i<N; ++i)
    {
        FEElement& el = dom.ElementRef(i);
        
        for (int k=0; k<nsols; ++k)
        {
            int nsid = lid[k];
            if (nsid == -1) a << 0.f;
            else
            {
                // calculate average relative Peclet number
                double ew = 0;
                for (int j = 0; j<el.GaussPoints(); ++j)
                {
                    FEMaterialPoint& mp = *el.GetMaterialPoint(j);
                    const FEFluidMaterialPoint* fpt = mp.ExtractData<FEFluidMaterialPoint>();
                    const FEElasticMaterialPoint* ept = mp.ExtractData<FEElasticMaterialPoint>();
                    vec3d v(0,0,0);
                    if (ept) v = ept->m_v;
                    ew += (fpt->m_vft - v).Length()/pm->GetFreeDiffusivity(mp, nsid);
                }
                ew /= el.GaussPoints();
                a << ew;
            }
        }
        
    }
    return true;
}

