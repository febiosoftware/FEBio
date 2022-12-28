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
#include "FEOsmoticCoefficient.h"
#include <FECore/FEFunction1D.h>

//-----------------------------------------------------------------------------
// This class implements a material that has an osmotic coefficient behaving
// according to the Wells-Manning theory.  The Wells correction is provided
// by a loadcurve.

class FEBIOMIX_API FEOsmCoefManning : public FEOsmoticCoefficient
{
public:
    //! constructor
    FEOsmCoefManning(FEModel* pfem);
    
    //! destructor
    ~FEOsmCoefManning() {}
    
    //! Initialization
    bool Init() override;
    
    //! osmotic coefficient
    double OsmoticCoefficient(FEMaterialPoint& pt) override;
    
    //! Tangent of osmotic coefficient with respect to strain (J=detF)
    double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp) override;
    
    //! Tangent of osmotic coefficient with respect to concentration
    double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp, const int isol) override;

    //! Manning response
    double OsmoticCoefficient_Manning(FEMaterialPoint& pt);
    double Tangent_OsmoticCoefficient_Strain_Manning(FEMaterialPoint& mp);
    double Tangent_OsmoticCoefficient_Concentration_Manning(FEMaterialPoint& mp, const int isol);

    //! Wells response
    double OsmoticCoefficient_Wells(FEMaterialPoint& pt);
    double Tangent_OsmoticCoefficient_Strain_Wells(FEMaterialPoint& mp);
    double Tangent_OsmoticCoefficient_Concentration_Wells(FEMaterialPoint& mp, const int isol);
    
public:
    double			m_ksi;	//!< Manning parameter
    int				m_sol;	//!< global id of co-ion
    int				m_lsol;	//!< local id of co-ion
    FEFunction1D*	m_osmc;	//!< osmotic coefficient for Wells correction (mobile ion - mobile interaction)

    // declare parameter list
    DECLARE_FECORE_CLASS();
};
