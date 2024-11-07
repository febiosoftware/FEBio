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
//! Elastic fluid model for real vapor
//!
class FEBIOFLUID_API FERealVapor : public FEElasticFluid
{
public:
    enum { MAX_NVP = 7, MAX_NVC = 4 };
    
public:
    FERealVapor(FEModel* pfem);
    ~FERealVapor() {}
    
    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;

    //! gauge pressure
    double Pressure(FEMaterialPoint& pt) override;
    
    //! specific free energy
    double SpecificFreeEnergy(FEMaterialPoint& mp) override;
    
    //! specific entropy
    double SpecificEntropy(FEMaterialPoint& mp) override;
    
    //! specific strain energy
    double SpecificStrainEnergy(FEMaterialPoint& mp) override;
    
    //! isochoric specific heat capacity
    double IsochoricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! isobaric specific heat capacity
    double IsobaricSpecificHeatCapacity(FEMaterialPoint& mp) override;
    
    //! dilatation from temperature and pressure
    bool Dilatation(const double T, const double p, double& e) override;

    //! tangent of pressure with respect to strain J
    double Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of pressure with respect to temperature T
    double Tangent_Temperature(FEMaterialPoint& mp) override;
    
public:
    double          m_Pr;           //!< referential absolute pressure
    double          m_Tr;           //!< referential absolute temperature
    double          m_Tc;           //!< normalized critical temperature (Tc/Tr)
    double          m_alpha;        //!< exponent alpha used for calculating temperature map
    double          m_rhor;         //!< referential mass density
    int             m_nvp;          //!< number of virial coefficients for pressure
    int             m_nvc;          //!< number of virial coefficients for isochoric specific heat capacity
    FEFunction1D*   m_D[MAX_NVP];   //!< non-dimensional pressure coefficient (multiply by m_Pr to get actual value)
    FEFunction1D*   m_psat;         //!< normalized gauge pressure on saturation curve (multiply by m_Pr to get actual value)
    FEFunction1D*   m_asat;         //!< normalized specific free energy on saturation curve (multiply by m_Pr/m_rhor to get actual value)
    FEFunction1D*   m_ssat;         //!< normalized specific entropy on saturation curve (multiply by m_Pr/m_rhor*m_Tr to get actual value)
    FEFunction1D*   m_esat;         //!< dilatation on saturation curve
    FEFunction1D*   m_cvsat;        //!< non-dimensional isochoric specific heat capacity on saturation curve (multiply by m_Pr/m_rhor*m_Tr to get actual value)
    FEFunction1D*   m_C[MAX_NVC];   //!< non-dimensional virial coefficients for isochoric specific heat capacity

    FEThermoFluid*  m_pMat; //!< parent thermo-fluid material
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
    
};
