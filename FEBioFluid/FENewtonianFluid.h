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
#include "FEViscousFluid.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Newtonian fluid

class FEBIOFLUID_API FENewtonianFluid :	public FEViscousFluid
{
public:
    //! constructor
    FENewtonianFluid(FEModel* pfem);
    
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
    
    //! dynamic viscosity
    double ShearViscosity(FEMaterialPoint& mp) override;
    double TangentShearViscosityTemperature(FEMaterialPoint& mp);

    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp) override;
    double TangentBulkViscosityTemperature(FEMaterialPoint& mp);

public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;		//!< shear viscosity
    double  m_Tr;       //!< referential temperature
    FEFunction1D*   m_kappahat; //!< normalized bulk viscosity vs normalized temperature
    FEFunction1D*   m_muhat;    //!< normalized shear viscosity vs normalized temperature

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
