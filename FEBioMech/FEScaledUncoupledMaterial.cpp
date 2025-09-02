/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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


#include "FEScaledUncoupledMaterial.h"

BEGIN_FECORE_CLASS(FEScaledUncoupledMaterial, FEUncoupledMaterial)
	ADD_PROPERTY(m_pBase, "solid");
    ADD_PARAMETER(m_scale, FE_RANGE_CLOSED(0,1), "scale");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! stress function
mat3ds FEScaledUncoupledMaterial::DevStress(FEMaterialPoint& pt)
{
    // get the elastic material point data
    FEElasticMaterialPoint& mp = *pt.ExtractData<FEElasticMaterialPoint>();
    
    double scale = m_scale(pt);
    return m_pBase->DevStress(pt)*scale;
}

//-----------------------------------------------------------------------------
//! tangent function
tens4ds FEScaledUncoupledMaterial::DevTangent(FEMaterialPoint& pt)
{
    // get the elastic material point data
    FEElasticMaterialPoint& mp = *pt.ExtractData<FEElasticMaterialPoint>();
    
    double scale = m_scale(pt);
    return m_pBase->DevTangent(pt)*scale;
}

//-----------------------------------------------------------------------------
//! strain energy density function
double FEScaledUncoupledMaterial::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    // get the elastic material point data
    FEElasticMaterialPoint& mp = *pt.ExtractData<FEElasticMaterialPoint>();
    
    double scale = m_scale(pt);
    return m_pBase->DevStrainEnergyDensity(pt)*scale;
}

FEMaterialPointData* FEScaledUncoupledMaterial::CreateMaterialPointData()
{
    return m_pBase->CreateMaterialPointData();
}