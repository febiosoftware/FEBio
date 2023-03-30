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

#include "FECore/FEMaterial.h"
#include "febiomech_api.h"

//---------------------------------------------------------------------------
// Base class for fiber density distribution functions
//
class FEBIOMECH_API FEFiberDensityDistribution : public FEMaterialProperty
{
public:
    FEFiberDensityDistribution(FEModel* pfem) : FEMaterialProperty(pfem) {}
    
    // Evaluation of fiber density along n0
    virtual double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) = 0;

    FECORE_BASE_CLASS(FEFiberDensityDistribution)
};

//---------------------------------------------------------------------------
// Spherical fiber density distribution
//
class FESphericalFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
    FESphericalFiberDensityDistribution(FEModel* pfem) : FEFiberDensityDistribution(pfem) {}
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override { return 1.0; }
};

//---------------------------------------------------------------------------
// Ellipsoidal fiber density distribution
//
class FEEllipsoidalFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
	FEEllipsoidalFiberDensityDistribution(FEModel* pfem);
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    FEParamVec3 m_spa;      // semi-principal axes of ellipsoid
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------
// 3D axisymmetric von Mises fiber density distribution
//
class FEVonMises3DFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
    FEVonMises3DFiberDensityDistribution(FEModel* pfem) : FEFiberDensityDistribution(pfem) { m_b = 0; }
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    FEParamDouble m_b;         // concentration parameter
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------
// 3D 2-fiber families axisymmetric von Mises fiber density distribution
//
class FEVonMises3DTwoFDDAxisymmetric : public FEFiberDensityDistribution
{
public:
	FEVonMises3DTwoFDDAxisymmetric(FEModel* pfem);
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    FEParamDouble	m_b;		// concentration parameter
    FEParamDouble	m_c;         // cosine of ±angle offset of fiber families
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------
// Circular fiber density distribution (2d)
//
class FECircularFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
    FECircularFiberDensityDistribution(FEModel* pfem) : FEFiberDensityDistribution(pfem) {}
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override { return 1.0; }
};

//---------------------------------------------------------------------------
// Elliptical fiber density distribution (2d)
//
class FEEllipticalFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
    FEEllipticalFiberDensityDistribution(FEModel* pfem) : FEFiberDensityDistribution(pfem) { m_spa[0] = 1; m_spa[1] = 1; }
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    FEParamDouble m_spa[2];    // semi-principal axes of ellipse
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------
// 2D planar von Mises fiber density distribution
//
class FEVonMises2DFiberDensityDistribution : public FEFiberDensityDistribution
{
public:
    FEVonMises2DFiberDensityDistribution(FEModel* pfem) : FEFiberDensityDistribution(pfem) { m_b = 0; }
    
    double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    FEParamDouble m_b;         // concentration parameter
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------
class FEStructureTensorDistribution : public FEFiberDensityDistribution
{
public:
	FEStructureTensorDistribution(FEModel* fem);

	double FiberDensity(FEMaterialPoint& mp, const vec3d& n0) override;

public:
	FEParamMat3ds	m_SPD;

	DECLARE_FECORE_CLASS();
};
