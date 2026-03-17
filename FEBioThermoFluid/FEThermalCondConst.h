/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FEFluidThermalConductivity.h"
#include <FECore/FEModel.h>
#include "febiothermofluid_api.h"

//-----------------------------------------------------------------------------
// This class implements a thermal conductivity which is constant

class FEBIOTHERMOFLUID_API FEThermalCondConst :	public FEFluidThermalConductivity
{
public:
	//! constructor
    FEThermalCondConst(FEModel* pfem);
		
    //! viscosity
    double NormalizedConductivity(FEMaterialPoint& pt) override { return 1.; };
    
    //! tangent of normalized viscosity with respect to temperature
    double Tangent_NormalizedConductivity_Temperature(FEMaterialPoint& mp) override { return 0; }

    //! tangent of normalized viscosity with respect to volumetric strain (or J)
    double Tangent_NormalizedConductivity_Strain(FEMaterialPoint& mp) override { return 0; }
};
