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
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Ideal gas under isentropic conditions.

class FEBIOFLUID_API FEIdealGasIsentropic : public FEFluid
{
public:
    FEIdealGasIsentropic(FEModel* pfem);
    
public:
    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! elastic pressure
    double Pressure(FEMaterialPoint& mp) override;
    double Pressure(const double e, const double T = 0) override;
    
    //! tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of elastic pressure with respect to strain J
    double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp) override;
    
    //! invert pressure-dilatation relation
    bool Dilatation(const double T, const double p, double& e) override;
    
    //! evaluate temperature
    double Temperature(FEMaterialPoint& mp) override;
    
public:
    double      m_gamma;    //!< ratio of specific heats (constant pressure/constant volume)
    double      m_M;        //!< molar mass
    double      m_Pr;       //!< ambient pressure
    double      m_Tr;       //!< ambient temperature
    double      m_R;        //!< universal gas constant
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
