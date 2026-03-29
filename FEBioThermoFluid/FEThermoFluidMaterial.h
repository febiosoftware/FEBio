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
#pragma once
#include <FEBioFluid/FEFluidMaterial.h>
#include "FEThermoViscousFluid.h"
#include "FEFluidHeatFlux.h"
#include "febiothermofluid_api.h"

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEBIOTHERMOFLUID_API FEThermoFluidMaterial : public FEFluidMaterial
{
public:
    FEThermoFluidMaterial(FEModel* pfem);
    virtual ~FEThermoFluidMaterial() {}

public:
    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override;
    
    //! return thermal conductivity
    FEFluidHeatFlux* GetHeatFlux() { return m_pHeatFlux; }

    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp)  { return m_pViscous->Tangent_RateOfDeformation(mp); }
    
    //! heat flux
    vec3d HeatFlux(FEMaterialPoint& mp);
    
private: // material properties
    FEThermoViscousFluid*   m_pViscous; //!< pointer to thermo viscous part of fluid material
    FEFluidHeatFlux*        m_pHeatFlux; //!< pointer to thermal conductivity material
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
