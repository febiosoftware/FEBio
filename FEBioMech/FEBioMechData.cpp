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
#include "FEBioMechData.h"
#include "FEElasticMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FEDamageMaterialPoint.h"
#include "FEReactivePlasticityMaterialPoint.h"
#include "FEReactivePlasticDamageMaterialPoint.h"
#include "FEElasticMixture.h"
#include "FEElasticMultigeneration.h"
#include "FERigidMaterial.h"
#include "FESolidSolver.h"
#include "FESolidSolver2.h"
#include "FERigidBody.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FERigidConnector.h"
#include "FEVolumeConstraint.h"
#include "FEContactSurface.h"
#include "FEDiscreteElasticMaterial.h"

//-----------------------------------------------------------------------------
double FENodeXPos::value(const FENode& node)
{
	return node.m_rt.x; 
}

//-----------------------------------------------------------------------------
double FENodeYPos::value(const FENode& node)
{
	return node.m_rt.y; 
}

//-----------------------------------------------------------------------------
double FENodeZPos::value(const FENode& node)
{
	return node.m_rt.z; 
}

//-----------------------------------------------------------------------------
double FENodeXDisp::value(const FENode& node)
{
	const int dof_X = GetFEModel()->GetDOFIndex("x");
	return node.get(dof_X); 
}

//-----------------------------------------------------------------------------
double FENodeYDisp::value(const FENode& node)
{
	const int dof_Y = GetFEModel()->GetDOFIndex("y");
	return node.get(dof_Y); 
}

//-----------------------------------------------------------------------------
double FENodeZDisp::value(const FENode& node)
{
	const int dof_Z = GetFEModel()->GetDOFIndex("z");
	return node.get(dof_Z); 
}

//-----------------------------------------------------------------------------
double FENodeXVel::value(const FENode& node)
{
	const int dof_VX = GetFEModel()->GetDOFIndex("vx");
	return node.get(dof_VX);
}

//-----------------------------------------------------------------------------
double FENodeYVel::value(const FENode& node)
{
	const int dof_VY = GetFEModel()->GetDOFIndex("vy");
	return node.get(dof_VY); 
}

//-----------------------------------------------------------------------------
double FENodeZVel::value(const FENode& node)
{
	const int dof_VZ = GetFEModel()->GetDOFIndex("vz");
	return node.get(dof_VZ);
}

//-----------------------------------------------------------------------------
double FENodeXAcc::value(const FENode& node)
{
	return node.m_at.x; 
}

//-----------------------------------------------------------------------------
double FENodeYAcc::value(const FENode& node)
{
	return node.m_at.y; 
}

//-----------------------------------------------------------------------------
double FENodeZAcc::value(const FENode& node)
{
	return node.m_at.z; 
}

//-----------------------------------------------------------------------------
double FENodeForceX::value(const FENode& node)
{
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		const vector<int>& id = node.m_ID;
		return (-id[0] - 2 >= 0 ? Fr[-id[0] - 2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceY::value(const FENode& node)
{
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		const vector<int>& id = node.m_ID;
		return (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceZ::value(const FENode& node)
{
	FESolidSolver2* psolid_solver = dynamic_cast<FESolidSolver2*>(GetFEModel()->GetCurrentStep()->GetFESolver());
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		const vector<int>& id = node.m_ID;
		return (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FELogContactGap::value(FESurfaceElement& el)
{
	double g = 0.0;
	for (int i = 0; i < el.GaussPoints(); ++i)
	{
		FEMaterialPoint* mp = el.GetMaterialPoint(i);
		FEContactMaterialPoint* cp = dynamic_cast<FEContactMaterialPoint*>(mp);
		if (cp)
		{
			g += cp->m_gap;
		}
	}
	g /= (double)el.GaussPoints();

	return g;
}

//-----------------------------------------------------------------------------
double FELogContactPressure::value(FESurfaceElement& el)
{
	double Lm = 0.0;
	for (int i = 0; i < el.GaussPoints(); ++i)
	{
		FEMaterialPoint* mp = el.GetMaterialPoint(i);
		FEContactMaterialPoint* cp = dynamic_cast<FEContactMaterialPoint*>(mp);
		if (cp)
		{
			Lm += cp->m_Ln;
		}
	}
	Lm /= (double)el.GaussPoints();

	return Lm;
}

//-----------------------------------------------------------------------------
double FELogContactTractionX::value(FESurfaceElement& el)
{
	FEContactSurface* ps = dynamic_cast<FEContactSurface*>(el.GetMeshPartition());
	if (ps == nullptr) return 0.0;

	FEContactInterface* pci = ps->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		vec3d tn;
		ps->GetSurfaceTraction(el.m_lid, tn);
		return tn.x;
	}
	return 0.0;
}

//-----------------------------------------------------------------------------
double FELogContactTractionY::value(FESurfaceElement& el)
{
	FEContactSurface* ps = dynamic_cast<FEContactSurface*>(el.GetMeshPartition());
	if (ps == nullptr) return 0.0;

	FEContactInterface* pci = ps->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		vec3d tn;
		ps->GetSurfaceTraction(el.m_lid, tn);
		return tn.y;
	}
	return 0.0;
}

//-----------------------------------------------------------------------------
double FELogContactTractionZ::value(FESurfaceElement& el)
{
	FEContactSurface* ps = dynamic_cast<FEContactSurface*>(el.GetMeshPartition());
	if (ps == nullptr) return 0.0;

	FEContactInterface* pci = ps->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		vec3d tn;
		ps->GetSurfaceTraction(el.m_lid, tn);
		return tn.z;
	}
	return 0.0;
}

//-----------------------------------------------------------------------------
double FELogElemPosX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
		val += pt.m_rt.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
		val += pt.m_rt.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
		val += pt.m_rt.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemJacobian::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_J;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainEffective::value(FEElement& el)
{
	int nint = el.GaussPoints();
	mat3ds Eavg; Eavg.zero();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

		mat3ds C = ep.LeftCauchyGreen();
		mat3dd I(1.0);
		mat3ds E = (C - I)*0.5;

		Eavg += E;
	}
	Eavg /= (double)nint;
	double val = Eavg.effective_norm();

	return val;
}

//-----------------------------------------------------------------------------
double FELogElemStrain2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xx();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.yy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.zz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.yz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemInfStrainXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds e = pt.SmallStrain();
		val += e.xz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.xx();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.yy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.zz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.xy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.yz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        val += U.xz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretch1::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        U.exact_eigen(l);
        val += l[0];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretch2::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        U.exact_eigen(l);
        val += l[1];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretch3::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds U = pt.RightStretch();
        U.exact_eigen(l);
        val += l[2];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightStretchEffective::value(FEElement& el)
{
    int nint = el.GaussPoints();
    mat3ds Uavg; Uavg.zero();
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
        
        mat3ds U = ep.RightStretch();
        
        Uavg += U;
    }
    Uavg /= (double)nint;
    double val = Uavg.effective_norm();
    
    return val;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.xx();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.yy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.zz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.xy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.yz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        val += V.xz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretch1::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        V.exact_eigen(l);
        val += l[0];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretch2::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        V.exact_eigen(l);
        val += l[1];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretch3::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds V = pt.LeftStretch();
        V.exact_eigen(l);
        val += l[2];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftStretchEffective::value(FEElement& el)
{
    int nint = el.GaussPoints();
    mat3ds Vavg; Vavg.zero();
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
        
        mat3ds V = ep.LeftStretch();
        
        Vavg += V;
    }
    Vavg /= (double)nint;
    double val = Vavg.effective_norm();
    
    return val;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.xx();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.yy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.zz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.xy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.yz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        val += H.xz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHencky1::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        H.exact_eigen(l);
        val += l[0];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHencky2::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        H.exact_eigen(l);
        val += l[1];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHencky3::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds H = pt.RightHencky();
        H.exact_eigen(l);
        val += l[2];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemRightHenckyEffective::value(FEElement& el)
{
    int nint = el.GaussPoints();
    mat3ds Havg; Havg.zero();
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
        
        mat3ds H = ep.RightHencky();
        
        Havg += H;
    }
    Havg /= (double)nint;
    double val = Havg.effective_norm();
    
    return val;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.xx();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.yy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.zz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.xy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.yz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        val += h.xz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHencky1::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        h.exact_eigen(l);
        val += l[0];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHencky2::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        h.exact_eigen(l);
        val += l[1];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHencky3::value(FEElement& el)
{
    double l[3];
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        mat3ds h = pt.LeftHencky();
        h.exact_eigen(l);
        val += l[2];
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemLeftHenckyEffective::value(FEElement& el)
{
    int nint = el.GaussPoints();
    mat3ds havg; havg.zero();
    for (int n = 0; n < nint; ++n)
    {
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
        
        mat3ds h = ep.LeftHencky();
        
        havg += h;
    }
    havg /= (double)nint;
    double val = havg.effective_norm();
    
    return val;
}

//-----------------------------------------------------------------------------
double FELogElemStressX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressEffective::value(FEElement& el)
{
	int nint = el.GaussPoints();
	mat3ds savg; savg.zero();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

		savg += ep.m_s;
	}
	savg /= (double)nint;
	double val = savg.effective_norm();

	return val;
}


//-----------------------------------------------------------------------------
double FELogElemStress1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.xx();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.yy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.zz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.xy();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.yz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemPK2StressXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		mat3ds S = pt.pull_back(pt.m_s);
		val += S.xz();
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressEigenVector::value(FEElement& el)
{
	assert(m_eigenVector >= 0);
	assert(m_component >= 0);

	// calculate average stress
	mat3ds s(0.0);
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		s += pt.m_s;
	}
	s /= (double)nint;

	// get the eigen vectors
	vec3d e[3];
	double l[3];
	s.eigen(l, e);

	// extract component
	double v = 0.0;
	switch (m_component)
	{
	case 0: v = e[m_eigenVector].x; break;
	case 1: v = e[m_eigenVector].y; break;
	case 2: v = e[m_eigenVector].z; break;
	}

	return v;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemElasticity_::value(FEElement& el, int n)
{
    FEElasticMaterial* pme = GetFEModel()->GetMaterial(el.GetMatID())->ExtractProperty<FEElasticMaterial>();
    if ((pme == 0) || pme->IsRigid()) return 0;

    tens4dmm c;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        c = pme->SolidTangent(pt);
		val += c.d[n];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainEnergyDensity::value(FEElement& el)
{
	FEElasticMaterial* pme = GetFEModel()->GetMaterial(el.GetMatID())->ExtractProperty<FEElasticMaterial>();
	if ((pme == 0) || pme->IsRigid()) return 0;
    
    double sed;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        sed = pme->StrainEnergyDensity(pt);
		val += sed;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDevStrainEnergyDensity::value(FEElement& el)
{
	FEElasticMaterial* pme = GetFEModel()->GetMaterial(el.GetMatID())->ExtractProperty<FEElasticMaterial>();
	FEUncoupledMaterial* pmu = dynamic_cast<FEUncoupledMaterial*>(pme);
    if ((pme == 0) || pme->IsRigid() || (pmu == 0)) return 0;
    
    double sed;
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(i);
        sed = pmu->DevStrainEnergyDensity(pt);
		val += sed;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFiberStretch::value(FEElement& el)
{
	int matID = el.GetMatID();
	FEMaterial* mat = GetFEModel()->GetMaterial(matID);

	int n = el.GaussPoints();
	double l = 0.0;
	for (int j=0; j<n; ++j)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(j);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3d Q = mat->GetLocalCS(mp);
		vec3d ri = Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.norm();
	}
	l /= (double) n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorX::value(FEElement& el)
{
	int matID = el.GetMatID();
	FEMaterial* mat = GetFEModel()->GetMaterial(matID);

	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(j);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3d Q = mat->GetLocalCS(mp);

		vec3d ri = Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.x;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorY::value(FEElement& el)
{
	int matID = el.GetMatID();
	FEMaterial* mat = GetFEModel()->GetMaterial(matID);

	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(j);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3d Q = mat->GetLocalCS(mp);

		vec3d ri = Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.y;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogElemFiberVectorZ::value(FEElement& el)
{
	int matID = el.GetMatID();
	FEMaterial* mat = GetFEModel()->GetMaterial(matID);

	int n = el.GaussPoints();
	double l = 0.0;
	for (int j = 0; j<n; ++j)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(j);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		mat3d Q = mat->GetLocalCS(mp);

		vec3d ri = Q.col(0);
		vec3d r = pt.m_F*ri;

		l += r.z;
	}
	l /= (double)n;
	return l;
}

//-----------------------------------------------------------------------------
double FELogDamage::value(FEElement& el)
{
    int nint = el.GaussPoints();
    double D = 0;
    for (int j=0; j<nint; ++j)
    {
        FEMaterialPoint& pt = *el.GetMaterialPoint(j);
        FEDamageMaterialPoint* ppd = pt.ExtractData<FEDamageMaterialPoint>();
        FEElasticMixtureMaterialPoint* pem = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        FEMultigenerationMaterialPoint* pmg = pt.ExtractData<FEMultigenerationMaterialPoint>();
        if (ppd) D += (float) ppd->m_D;
        else if (pem) {
            for (int k=0; k<pem->Components(); ++k)
            {
                FEDamageMaterialPoint* ppd = pem->GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
                if (ppd) D += (float) ppd->m_D;
            }
        }
        else if (pmg) {
            for (int k=0; k<pmg->Components(); ++k)
            {
                FEDamageMaterialPoint* ppd = pmg->GetPointData(k)->ExtractData<FEDamageMaterialPoint>();
                FEElasticMixtureMaterialPoint* pem = pmg->GetPointData(k)->ExtractData<FEElasticMixtureMaterialPoint>();
                if (ppd) D += (float) ppd->m_D;
                else if (pem)
                {
                    for (int l=0; l<pem->Components(); ++l)
                    {
                        FEDamageMaterialPoint* ppd = pem->GetPointData(l)->ExtractData<FEDamageMaterialPoint>();
                        if (ppd) D += (float) ppd->m_D;
                    }
                }
            }
        }
    }
    D /= (double) nint;
    return D;
}

//-----------------------------------------------------------------------------
double FELogOctahedralPlasticStrain::value(FEElement& el)
{
    int nint = el.GaussPoints();
    double D = 0;
    for (int j=0; j<nint; ++j)
    {
        FEMaterialPoint& pt = *el.GetMaterialPoint(j);
        FEReactivePlasticityMaterialPoint* prp = pt.ExtractData<FEReactivePlasticityMaterialPoint>();
        FEReactivePlasticDamageMaterialPoint* prd = pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
        FEElasticMixtureMaterialPoint* pem = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        FEMultigenerationMaterialPoint* pmg = pt.ExtractData<FEMultigenerationMaterialPoint>();
        if (prp) D += (float) prp->m_gp[0];
        else if (prd) D += (float) prd->m_gp[0];
        else if (pem) {
            for (int k=0; k<pem->Components(); ++k)
            {
                FEReactivePlasticityMaterialPoint* prp = pt.ExtractData<FEReactivePlasticityMaterialPoint>();
                FEReactivePlasticDamageMaterialPoint* prd = pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
                if (prp) D += (float) prp->m_gp[0];
                else if (prd) D += (float) prd->m_gp[0];
            }
        }
        else if (pmg) {
            for (int k=0; k<pmg->Components(); ++k)
            {
                FEReactivePlasticityMaterialPoint* prp = pt.ExtractData<FEReactivePlasticityMaterialPoint>();
                FEReactivePlasticDamageMaterialPoint* prd = pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
                FEElasticMixtureMaterialPoint* pem = pmg->GetPointData(k)->ExtractData<FEElasticMixtureMaterialPoint>();
                if (prp) D += (float) prp->m_gp[0];
                else if (prd) D += (float) prd->m_gp[0];
                else if (pem)
                {
                    for (int l=0; l<pem->Components(); ++l)
                    {
                        FEReactivePlasticityMaterialPoint* prp = pt.ExtractData<FEReactivePlasticityMaterialPoint>();
                        FEReactivePlasticDamageMaterialPoint* prd = pt.ExtractData<FEReactivePlasticDamageMaterialPoint>();
                        if (prp) D += (float) prp->m_gp[0];
                        else if (prd) D += (float) prd->m_gp[0];
                    }
                }
            }
        }
    }
    D /= (double) nint;
    return D;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementStretch::value(FEElement& el)
{
	if (dynamic_cast<FEDiscreteElement*>(&el) == nullptr) return 0.0;

	FEDiscreteElement& del = dynamic_cast<FEDiscreteElement&>(el);

	FEMesh& mesh = GetFEModel()->GetMesh();

	vec3d ra0 = mesh.Node(del.m_node[0]).m_r0;
	vec3d ra1 = mesh.Node(del.m_node[0]).m_rt;
	vec3d rb0 = mesh.Node(del.m_node[1]).m_r0;
	vec3d rb1 = mesh.Node(del.m_node[1]).m_rt;

	double L0 = (rb0 - ra0).norm();
	double Lt = (rb1 - ra1).norm();

	double l = Lt / L0;

	return l;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementElongation::value(FEElement& el)
{
	if (dynamic_cast<FEDiscreteElement*>(&el) == nullptr) return 0.0;

	FEDiscreteElement& del = dynamic_cast<FEDiscreteElement&>(el);

	FEMesh& mesh = GetFEModel()->GetMesh();

	vec3d ra0 = mesh.Node(del.m_node[0]).m_r0;
	vec3d ra1 = mesh.Node(del.m_node[0]).m_rt;
	vec3d rb0 = mesh.Node(del.m_node[1]).m_r0;
	vec3d rb1 = mesh.Node(del.m_node[1]).m_rt;

	double L0 = (rb0 - ra0).norm();
	double Lt = (rb1 - ra1).norm();

	double Dl = Lt - L0;

	return Dl;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementForce::value(FEElement& el)
{
	if (dynamic_cast<FEDiscreteElement*>(&el) == nullptr) return 0.0;

	FEDiscreteElement& del = dynamic_cast<FEDiscreteElement&>(el);
	FEMesh& mesh = GetFEModel()->GetMesh();

	// get the (one) material point data
	FEDiscreteElasticMaterialPoint& mp = dynamic_cast<FEDiscreteElasticMaterialPoint&>(*el.GetMaterialPoint(0));

	vec3d ra1 = mesh.Node(del.m_node[0]).m_rt;
	vec3d rb1 = mesh.Node(del.m_node[1]).m_rt;
	vec3d e = rb1 - ra1; e.unit();

	vec3d F = mp.m_Ft;

	double Fm = F * e;

	return Fm;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementForceX::value(FEElement& el)
{
	FEDiscreteElasticMaterialPoint* mp = dynamic_cast<FEDiscreteElasticMaterialPoint*>(el.GetMaterialPoint(0));
	if (mp) return mp->m_Ft.x;
	else return 0.0;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementForceY::value(FEElement& el)
{
	FEDiscreteElasticMaterialPoint* mp = dynamic_cast<FEDiscreteElasticMaterialPoint*>(el.GetMaterialPoint(0));
	if (mp) return mp->m_Ft.y;
	else return 0.0;
}

//-----------------------------------------------------------------------------
double FELogDiscreteElementForceZ::value(FEElement& el)
{
	FEDiscreteElasticMaterialPoint* mp = dynamic_cast<FEDiscreteElasticMaterialPoint*>(el.GetMaterialPoint(0));
	if (mp) return mp->m_Ft.z;
	else return 0.0;
}

//-----------------------------------------------------------------------------
double FELogElementMixtureStress::value(FEElement& el)
{
	if (m_comp < 0) return 0.0;

	double s = 0.0;
	for (int n = 0; n < el.GaussPoints(); ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMixtureMaterialPoint* mmp = mp.ExtractData< FEElasticMixtureMaterialPoint>();
		if (mmp)
		{
			if (m_comp < mmp->Components())
			{
				FEElasticMaterialPoint& ep = *mmp->GetPointData(m_comp)->ExtractData<FEElasticMaterialPoint>();

				switch (m_metric)
				{
				case 0: s += ep.m_s.xx(); break;
				case 1: s += ep.m_s.xy(); break;
				case 2: s += ep.m_s.yy(); break;
				case 3: s += ep.m_s.xz(); break;
				case 4: s += ep.m_s.yz(); break;
				case 5: s += ep.m_s.zz(); break;
				}
			}
		}
	}
	s /= (double)el.GaussPoints();

	return s;
}

//-----------------------------------------------------------------------------
double FELogRigidBodyR11::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0,0)); }
double FELogRigidBodyR12::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0, 1)); }
double FELogRigidBodyR13::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(0, 2)); }
double FELogRigidBodyR21::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 0)); }
double FELogRigidBodyR22::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 1)); }
double FELogRigidBodyR23::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(1, 2)); }
double FELogRigidBodyR31::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 0)); }
double FELogRigidBodyR32::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 1)); }
double FELogRigidBodyR33::value(FERigidBody& rb) { return (rb.GetRotation().RotationMatrix()(2, 2)); }

//-----------------------------------------------------------------------------
double FELogRigidBodyPosX::value(FERigidBody& rb) { return rb.m_rt.x; }
double FELogRigidBodyPosY::value(FERigidBody& rb) { return rb.m_rt.y; }
double FELogRigidBodyPosZ::value(FERigidBody& rb) { return rb.m_rt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyVelX::value(FERigidBody& rb) { return rb.m_vt.x; }
double FELogRigidBodyVelY::value(FERigidBody& rb) { return rb.m_vt.y; }
double FELogRigidBodyVelZ::value(FERigidBody& rb) { return rb.m_vt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAccX::value(FERigidBody& rb) { return rb.m_at.x; }
double FELogRigidBodyAccY::value(FERigidBody& rb) { return rb.m_at.y; }
double FELogRigidBodyAccZ::value(FERigidBody& rb) { return rb.m_at.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngPosX::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).x*rb.GetRotation().GetAngle()); }
double FELogRigidBodyAngPosY::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).y*rb.GetRotation().GetAngle()); }
double FELogRigidBodyAngPosZ::value(FERigidBody& rb) { return ((rb.GetRotation().GetVector()).z*rb.GetRotation().GetAngle()); }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngVelX::value(FERigidBody& rb) { return rb.m_wt.x; }
double FELogRigidBodyAngVelY::value(FERigidBody& rb) { return rb.m_wt.y; }
double FELogRigidBodyAngVelZ::value(FERigidBody& rb) { return rb.m_wt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyAngAccX::value(FERigidBody& rb) { return rb.m_alt.x; }
double FELogRigidBodyAngAccY::value(FERigidBody& rb) { return rb.m_alt.y; }
double FELogRigidBodyAngAccZ::value(FERigidBody& rb) { return rb.m_alt.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyQuatX::value(FERigidBody& rb) { return rb.GetRotation().x; }
double FELogRigidBodyQuatY::value(FERigidBody& rb) { return rb.GetRotation().y; }
double FELogRigidBodyQuatZ::value(FERigidBody& rb) { return rb.GetRotation().z; }
double FELogRigidBodyQuatW::value(FERigidBody& rb) { return rb.GetRotation().w; }

//-----------------------------------------------------------------------------
double FELogRigidBodyForceX::value(FERigidBody& rb) { return rb.m_Fr.x; }
double FELogRigidBodyForceY::value(FERigidBody& rb) { return rb.m_Fr.y; }
double FELogRigidBodyForceZ::value(FERigidBody& rb) { return rb.m_Fr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyTorqueX::value(FERigidBody& rb) { return rb.m_Mr.x; }
double FELogRigidBodyTorqueY::value(FERigidBody& rb) { return rb.m_Mr.y; }
double FELogRigidBodyTorqueZ::value(FERigidBody& rb) { return rb.m_Mr.z; }

//-----------------------------------------------------------------------------
double FELogRigidBodyKineticEnergy::value(FERigidBody& rb) {
    FERigidBody&rbl = static_cast<FERigidBody&>(rb);
    return (rbl.m_mass*(rbl.m_vt*rbl.m_vt) + rbl.m_wt*(rbl.m_moi*rbl.m_wt))/2;
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorForceX::value(FENLConstraint& rc) 
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.x : 0); 
}

double FELogRigidConnectorForceY::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.y : 0);
}

double FELogRigidConnectorForceZ::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_F.z : 0);
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorMomentX::value(FENLConstraint& rc)
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.x : 0);
}

double FELogRigidConnectorMomentY::value(FENLConstraint& rc)
{ 
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.y : 0);
}

double FELogRigidConnectorMomentZ::value(FENLConstraint& rc)
{
	FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
	return (prc ? prc->m_M.z : 0);
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorTranslationX::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeTranslation().x : 0);
}

double FELogRigidConnectorTranslationY::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeTranslation().y : 0);
}

double FELogRigidConnectorTranslationZ::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeTranslation().z : 0);
}

//-----------------------------------------------------------------------------
double FELogRigidConnectorRotationX::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeRotation().x : 0);
}

double FELogRigidConnectorRotationY::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeRotation().y : 0);
}

double FELogRigidConnectorRotationZ::value(FENLConstraint& rc)
{
    FERigidConnector* prc = dynamic_cast<FERigidConnector*>(&rc);
    return (prc ? prc->RelativeRotation().z : 0);
}

//-----------------------------------------------------------------------------
double FELogVolumeConstraint::value(FENLConstraint& rc)
{
    FEVolumeConstraint* prc = dynamic_cast<FEVolumeConstraint*>(&rc);
    return (prc ? prc->EnclosedVolume() : 0);
}

//-----------------------------------------------------------------------------
double FELogVolumePressure::value(FENLConstraint& rc)
{
    FEVolumeConstraint* prc = dynamic_cast<FEVolumeConstraint*>(&rc);
    return (prc ? prc->Pressure() : 0);
}

//=============================================================================
double FELogContactArea::value(FESurface& surface)
{
	FEContactSurface* pcs = dynamic_cast<FEContactSurface*>(&surface);
	if (pcs == 0) return 0.0;

	// make sure the corresponding contact interface is active
	// (in case the parent was not set, we'll proceed regardless)
	FEContactInterface* pci = pcs->GetContactInterface(); assert(pci);
	if ((pci == 0) || pci->IsActive())
	{
		double area = pcs->GetContactArea();
		return area;
	}
	return 0.0;
}
