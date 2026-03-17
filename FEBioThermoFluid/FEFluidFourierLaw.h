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
#include "FEFluidHeatFlux.h"
#include "febiothermofluid_api.h"
#include "FEFluidThermalConductivity.h"

//-----------------------------------------------------------------------------
// This class evaluates the heat flux from Fourier's law for a fluid

class FEBIOTHERMOFLUID_API FEFluidFourierLaw : public FEFluidHeatFlux
{
public:
    //! constructor
    FEFluidFourierLaw(FEModel* pfem);
    
    //! initialize
    bool Init() override;
    
    //! heat flux
    vec3d HeatFlux(FEMaterialPoint& pt) override;

    double Conductivity(FEMaterialPoint& mp) override { return m_Kr*m_Khat->NormalizedConductivity(mp); }
    
    double Tangent_Conductivity_Temperature(FEMaterialPoint& mp) override { return m_Kr*m_Khat->Tangent_NormalizedConductivity_Temperature(mp); }
    
    double Tangent_Conductivity_Strain(FEMaterialPoint& mp) override { return m_Kr*m_Khat->Tangent_NormalizedConductivity_Strain(mp); }

public:
    double                      m_Kr;   // referential thermal conductivity
    FEFluidThermalConductivity* m_Khat; // normalized thermal conductivity

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
