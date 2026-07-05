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
#include "FEThermoViscousFluid.h"
#include "febiothermofluid_api.h"
#include "FEThermoFluidProperty.h"
#include <FEBioFluid/FENewtonianFluid.h>
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Newtonian fluid

class FEBIOTHERMOFLUID_API FENewtonianThermoFluid :	public FEThermoViscousFluid
{
public:
    //! constructor
    FENewtonianThermoFluid(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! viscous stress
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to temperature
    mat3ds Tangent_Temperature(FEMaterialPoint& mp) override;
    
    //! dynamic shear viscosity
    double ShearViscosity(FEMaterialPoint& mp) override;
    double TangentShearViscosityTemperature(FEMaterialPoint& mp);
    double TangentShearViscosityStrain(FEMaterialPoint& mp);

    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp) override;
    double TangentBulkViscosityTemperature(FEMaterialPoint& mp);
    double TangentBulkViscosityStrain(FEMaterialPoint& mp);

public:
    // viscosities
    double              m_kappa;    //!< referential bulk viscosity
    double              m_mu;       //!< referential shear viscosity
    FEThermoFluidProperty*  m_kappahat; //!< normalized bulk viscosity
    FEThermoFluidProperty*  m_muhat;    //!< normalized shear viscosity

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
