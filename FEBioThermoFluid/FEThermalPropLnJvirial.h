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
#include "FEThermoFluidProperty.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
// This class implements a thermofluid property which is a virial expansion in lnJ

class FEBIOTHERMOFLUID_API FEThermalPropLnJvirial :	public FEThermoFluidProperty
{
public:
    enum { MAX_COEF = 7 };
    
public:
	//! constructor
    FEThermalPropLnJvirial(FEModel* pfem);
    
    //! initialize
    bool Init() override;
		
    //! viscosity
    double NormalizedProperty(FEMaterialPoint& pt) override;
    
    //! tangent of normalized viscosity with respect to temperature
    double Tangent_NormalizedProperty_Temperature(FEMaterialPoint& mp) override;

    //! tangent of normalized viscosity with respect to volumetric strain (or J)
    double Tangent_NormalizedProperty_Strain(FEMaterialPoint& mp) override;
    
public:
	FEFunction1D*	m_coef[MAX_COEF];			//!< virial coefficients
    double          m_Pr;           //!< referential absolute pressure
    double          m_Tr;           //!< referential absolute temperature

private:
    int             m_ncoef;                    //!< number of coefficients in the virial expansion
		
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
