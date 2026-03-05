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
#include "FENonlinearElasticFluid.h"
#include "FEFluidMaterialPoint.h"
#include "FEFluid.h"

//-----------------------------------------------------------------------------
FENonlinearElasticFluid::FENonlinearElasticFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
}

//-----------------------------------------------------------------------------
//! initialization
bool FENonlinearElasticFluid::Init()
{
    return true;
}

//-----------------------------------------------------------------------------
//! gauge pressure
double FENonlinearElasticFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
    double p =k*(1/(1+fp.m_ef)-1);
    
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of pressure with respect to strain J
double FENonlinearElasticFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    return -m_pFluid->m_k /pow(1+fp.m_ef,2);
}

//-----------------------------------------------------------------------------
//! 2nd tangent of pressure with respect to strain J
double FENonlinearElasticFluid::Tangent_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    return 2* m_pFluid->m_k /pow(1+fp.m_ef,3);
}

//-----------------------------------------------------------------------------
//! specific free energy
double FENonlinearElasticFluid::SpecificFreeEnergy(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
	double k = m_pFluid->m_k;
	double rhor = m_pFluid->m_rhor;
    double a = k/rhor*(fp.m_ef-log(1+fp.m_ef));
    return a;
}

//-----------------------------------------------------------------------------
//! dilatation from temperature and pressure
bool FENonlinearElasticFluid::Dilatation(const double T, const double p, double& e)
{
    e = 1.0/(1+p/ m_pFluid->m_k)-1.0;
    return true;
}
