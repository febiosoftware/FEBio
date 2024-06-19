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
// This class evaluates the viscous stress in a Newtonian viscous real vapor

class FEBIOFLUID_API FENewtonianRealVapor :	public FEViscousFluid
{
public:
    enum { MAX_NVC = 4 };
    
public:
    //! constructor
    FENewtonianRealVapor(FEModel* pfem);
    
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
    double TangentShearViscosityStrain(FEMaterialPoint& mp);

    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp) override;
    double TangentBulkViscosityTemperature(FEMaterialPoint& mp);
    double TangentBulkViscosityStrain(FEMaterialPoint& mp);

public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;		//!< shear viscosity
    double  m_Tr;       //!< referential temperature
    double          m_Tc;           //!< normalized critical temperature (Tc/Tr)
    double          m_alpha;        //!< exponent alpha used for calculating temperature map
    int             m_nvc;          //!< number of virial coefficients for isochoric specific heat capacity
    FEFunction1D*   m_esat;         //!< dilatation on saturation curve
    FEFunction1D*   m_musat;        //!< normalized shear viscosity vs normalized temperature on saturation curve
    FEFunction1D*   m_C[MAX_NVC];   //!< non-dimensional pressure coefficient (multiply by m_Pr to get actual value)

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
