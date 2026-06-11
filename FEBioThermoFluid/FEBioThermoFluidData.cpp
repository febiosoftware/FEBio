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



#include "FEBioThermoFluidData.h"
#include "FEThermoFluid.h"
#include "FEThermoElasticFluid.h"
#include <FEBioFluid/FEFluid.h>
#include <FECore/FEModel.h>

//=============================================================================
// N O D E  D A T A
//=============================================================================

//-----------------------------------------------------------------------------
double FENodeThermoFluidTemperature::value(const FENode& node)
{
    const int dof_T = GetFEModel()->GetDOFIndex("T");
    return node.get(dof_T);
}

//=============================================================================
// E L E M E N T   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
double FELogThermoElasticFluidPressure::value(FEElement& el)
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
double FELogThermoFluidVolumeRatio::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    for (int i=0; i<nint; ++i)
    {
        FEFluidMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEFluidMaterialPoint>();
        if (ppt) val += ppt->m_ef + 1;
    }
    return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidDensity::value(FEElement& el)
{
    double val = 0.0;
    int nint = el.GaussPoints();
    FEMaterial* pmat = GetFEModel()->GetMaterial(el.GetMatID());
    FEFluid* pfmat = dynamic_cast<FEFluid*>(pmat);
    FEThermoFluid* ptfmat = dynamic_cast<FEThermoFluid*>(pmat);
    if (pfmat) {
        for (int i=0; i<nint; ++i)
        {
            FEMaterialPoint* pt = el.GetMaterialPoint(i);
            val += pfmat->Density(*pt);
        }
        return val / (double) nint;
    }
    else if (ptfmat) {
        for (int i=0; i<nint; ++i)
        {
            FEMaterialPoint* pt = el.GetMaterialPoint(i);
            val += ptfmat->Density(*pt);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidSpecificFreeEnergy::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->SpecificFreeEnergy(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->SpecificFreeEnergy(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidSpecificEntropy::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->SpecificEntropy(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->SpecificEntropy(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidSpecificInternalEnergy::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->SpecificInternalEnergy(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->SpecificInternalEnergy(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidSpecificEnthalpy::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->SpecificEnthalpy(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->SpecificGaugeEnthalpy(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidSpecificStrainEnergy::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->SpecificStrainEnergy(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->SpecificStrainEnergy(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidIsochoricSpecificHeatCapacity::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->IsochoricSpecificHeatCapacity(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->IsochoricSpecificHeatCapacity(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidIsobaricSpecificHeatCapacity::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->IsobaricSpecificHeatCapacity(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->IsobaricSpecificHeatCapacity(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

//-----------------------------------------------------------------------------
double FELogThermoFluidTemperature::value(FEElement& el)
{
    double val = 0.0;
    FEModel* fem = GetFEModel();
    FEMaterial* pmat = fem->GetMaterial(el.GetMatID());
    FEElasticFluid* pfluid = pmat->ExtractProperty<FEElasticFluid>();
    FEThermoElasticFluid* tfluid = pmat->ExtractProperty<FEThermoElasticFluid>();
    if (pfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += pfluid->Temperature(mp);
        }
        return val / (double) nint;
    }
    else if (tfluid) {
        int nint = el.GaussPoints();
        for (int i=0; i<nint; ++i) {
            FEMaterialPoint& mp = *el.GetMaterialPoint(i);
            val += tfluid->Temperature(mp);
        }
        return val / (double) nint;
    }
    else
        return 0;
}

