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
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a porosity and concentration
// dependent diffusivity which is isotropic, according to the constitutive relation
// of Albro et al (CMBE 2009)

class FEBIOMIX_API FEDiffAlbroIso : public FESoluteDiffusivity
{
public:
    //! constructor
    FEDiffAlbroIso(FEModel* pfem);
    
    //! free diffusivity
    double Free_Diffusivity(FEMaterialPoint& pt) override;
    
    //! Tangent of free diffusivity with respect to concentration
    double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& pt, const int isol) override;
    
    //! diffusivity
    mat3ds Diffusivity(FEMaterialPoint& pt) override;
    
    //! Tangent of diffusivity with respect to strain
    tens4dmm Tangent_Diffusivity_Strain(FEMaterialPoint& mp) override;
    
    //! Tangent of diffusivity with respect to concentration
    mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol=0) override;
    
    //! data initialization and checking
    bool Init() override;
    
public:
    double	m_diff0;        //!< free diffusivity
    double	m_cdinv;		//!< inverse of characteristic concentration c_D
    double	m_alphad;		//!< non-dimensional coefficient for porosity term
    int     m_lsol;         //!< local ID of solute
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
