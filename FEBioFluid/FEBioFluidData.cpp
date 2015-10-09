#include "stdafx.h"
#include "FEBioFluidData.h"
#include "FEFluid.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
double FELogElemFluidPosX::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        val += pt.m_r0.x;
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
        FEFluidMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        val += pt.m_r0.y;
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
        FEFluidMaterialPoint& pt = *el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        val += pt.m_r0.z;
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
		if (ppt) val += ppt->m_p;
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
        if (ppt) val += ppt->m_J;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFluidDensity::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    FEMaterial* pmat = m_pfem->GetMaterial(el.GetMatID());
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
        if (ppt) val += (ppt->m_s*ppt->m_L).trace();
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
        if (ppt) val += ppt->m_vt.x;
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
        if (ppt) val += ppt->m_vt.y;
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
        if (ppt) val += ppt->m_vt.z;
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
        if (ppt) val += ppt->m_at.x;
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
        if (ppt) val += ppt->m_at.y;
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
        if (ppt) val += ppt->m_at.z;
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
        if (ppt) val += ppt->m_s.xx();
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
        if (ppt) val += ppt->m_s.yy();
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
        if (ppt) val += ppt->m_s.zz();
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
        if (ppt) val += ppt->m_s.xy();
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
        if (ppt) val += ppt->m_s.yz();
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
        if (ppt) val += ppt->m_s.xz();
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

