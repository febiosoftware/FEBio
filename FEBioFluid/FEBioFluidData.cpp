#include "stdafx.h"
#include "FEBioFluidData.h"
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
double FENodeFluidXVel::value(int nnode)
{
    const int dof_VFX = GetFEModel()->GetDOFIndex("wx");
    FEMesh& mesh = GetFEModel()->GetMesh();
    FENode& node = mesh.Node(nnode);
    return node.get(dof_VFX);
}

//-----------------------------------------------------------------------------
double FENodeFluidYVel::value(int nnode)
{
    const int dof_VFY = GetFEModel()->GetDOFIndex("wy");
    FEMesh& mesh = GetFEModel()->GetMesh();
    FENode& node = mesh.Node(nnode);
    return node.get(dof_VFY);
}

//-----------------------------------------------------------------------------
double FENodeFluidZVel::value(int nnode)
{
    const int dof_VFZ = GetFEModel()->GetDOFIndex("wz");
    FEMesh& mesh = GetFEModel()->GetMesh();
    FENode& node = mesh.Node(nnode);
    return node.get(dof_VFZ);
}

//-----------------------------------------------------------------------------
double FELogElemFluidPosX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* pt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        FEElasticMaterialPoint* ept = el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        if (pt) val += pt->m_r0.x;
        else if (ept) val += ept->m_rt.x;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidPosY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* pt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        FEElasticMaterialPoint* ept = el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        if (pt) val += pt->m_r0.y;
        else if (ept) val += ept->m_rt.y;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidPosZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* pt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        FEElasticMaterialPoint* ept = el.GetMaterialPoint(i)->ExtractData<FEElasticMaterialPoint>();
        if (pt) val += pt->m_r0.z;
        else if (ept) val += ept->m_rt.z;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElasticFluidPressure::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
		if (ppt) val += ppt->m_pf;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVolumeRatio::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_Jf;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidDensity::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    FEMaterial* pmat = GetFEModel()->GetMaterial(el.GetMatID());
    FEFluid* pfmat = dynamic_cast<FEFluid*>(pmat);
    if (pfmat) {
        for (int i=0; i<nint; ++i)
        {
            FEMaterialPoint* pt = el.GetMaterialPoint(i);
            val += pfmat->Density(*pt);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogFluidStressPower::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += (ppt->m_sf*ppt->m_Lf).trace();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVelocityX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_vft.x;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVelocityY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_vft.y;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVelocityZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_vft.z;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidAccelerationX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_aft.x;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidAccelerationY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_aft.y;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidAccelerationZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_aft.z;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVorticityX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->Vorticity().x;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVorticityY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->Vorticity().y;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidVorticityZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->Vorticity().z;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressXX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.xx();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressYY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.yy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressZZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.zz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.xy();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.yz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidStressXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_sf.xz();
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefXX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.xx();
        }
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefYY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.yy();
        }
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefZZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.zz();
        }
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefXY::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.xy();
        }
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefYZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.yz();
        }
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidRateOfDefXZ::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) {
            mat3ds D = ppt->RateOfDeformation();
            val += D.xz();
        }
    }
    return val / (double) nint;
}

