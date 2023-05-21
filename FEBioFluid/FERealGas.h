/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEElasticFluid.h"
#include "FEThermoFluid.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
//! Base class for the viscous part of the fluid response.
//! These materials provide the viscous stress and its tangents.
//!
class FEBIOFLUID_API FERealGas : public FEElasticFluid
{
public:
    enum { MAX_NVA = 7, MAX_NVC = 4 };
    
public:
    FERealGas(FEModel* pfem);
    ~FERealGas() {}
    
    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! gage pressure
    double Pressure(FEMaterialPoint& pt) override;
    
    //! tangent of pressure with respect to strain J
    double Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of pressure with respect to strain J
    double Tangent_Strain_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of pressure with respect to temperature T
    double Tangent_Temperature(FEMaterialPoint& mp) override;
    
    //! 2nd tangent of pressure with respect to temperature T
    double Tangent_Temperature_Temperature(FEMaterialPoint& mp) override;
    
    //! tangent of pressure with respect to strain J and temperature T
    double Tangent_Strain_Temperature(FEMaterialPoint& mp) override;
    
    //! specific free energy
    double SpecificFreeEnergy(FEMaterialPoint& mp) override;
    
    //! specific entropy
    double SpecificEntropy(FEMaterialPoint& mp) override;
    
    //! specific strain energy
    double SpecificStrainEnergy(FEMaterialPoint& mp) override;
    
    //! isochoric specific heat capacity
    double IsochoricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! tangent of isochoric specific heat capacity with respect to strain J
    double Tangent_cv_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of isochoric specific heat capacity with respect to temperature T
    double Tangent_cv_Temperature(FEMaterialPoint& mp) override;
    
    //! isobaric specific heat capacity
    double IsobaricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! dilatation from temperature and pressure
    bool Dilatation(const double T, const double p, double& e) override;

public:
    double      m_R;        //!< universal gas constant
    double      m_Pr;       //!< referential absolute pressure
    double      m_Tr;       //!< referential absolute temperature
    double      m_rhor;     //!< referential mass density
    int             m_nva;          //!< number of virial coefficients for pressure constitutive relation
    int             m_nvc;          //!< number of virial coefficients for isochoric specific heat capacity
    FEFunction1D*   m_a0;           //!< specific free energy at zero pressure
    FEFunction1D*   m_A[MAX_NVA];   //!< non-dimensional virial coefficients for pressure constitutive relation
    FEFunction1D*   m_C[MAX_NVC];   //!< non-dimensional virial coefficients for isochoric specific heat capacity

    FEThermoFluid*  m_pMat; //!< parent thermo-fluid material
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
    
};
