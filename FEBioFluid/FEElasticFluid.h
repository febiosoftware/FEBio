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
#include <FECore/FEMaterial.h>
#include <FECore/tens4d.h>
#include "febiofluid_api.h"

class FEFluidMaterial;

//-----------------------------------------------------------------------------
//! Base class for the viscous part of the fluid response.
//! These materials provide the viscous stress and its tangents.
//!
class FEBIOFLUID_API FEElasticFluid : public FEMaterialProperty
{
public:
    FEElasticFluid(FEModel* pfem) : FEMaterialProperty(pfem) {}
    virtual ~FEElasticFluid() {}

	void SetParentFluid(FEFluidMaterial* pf) { m_pFluid = pf; }
    
    //! gauge pressure
    virtual double Pressure(FEMaterialPoint& pt) = 0;
    
    //! tangent of pressure with respect to strain J
    virtual double Tangent_Strain(FEMaterialPoint& mp);
    
    //! 2nd tangent of pressure with respect to strain J
    virtual double Tangent_Strain_Strain(FEMaterialPoint& mp);
    
    //! specific free energy
    virtual double SpecificFreeEnergy(FEMaterialPoint& mp) = 0;
    
    //! calculate dilatation for given (effective) pressure and temperature
    virtual bool Dilatation(const double T, const double p, double& e) = 0;
    
    //! calculate fluid pressure and its derivatives from state variables
    double Pressure(const double ef, const double T);
    double Tangent_Strain(const double ef, const double T);

protected:
	FEFluidMaterial* m_pFluid = nullptr;	//!< pointer to parent fluid material

    FECORE_BASE_CLASS(FEElasticFluid)
};
