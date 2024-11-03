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
#include "febiomech_api.h"
#include "FENonLocalKernel.h"
#include <FECore/FEMaterial.h>
#include <FECore/FEMeshTopo.h>

//-----------------------------------------------------------------------------
// Virtual base class for nonlocal averaging method
class FEBIOMECH_API FENonLocalAveraging : public FEMaterialProperty
{
public:
    FENonLocalAveraging(FEModel* pfem) : FEMaterialProperty(pfem) { m_krnl = nullptr; m_binit = false; }
    virtual ~FENonLocalAveraging() {}
    
    virtual bool Init() override;
    
    //! damage criterion average
    virtual double DamageCriterionAverage(FEMaterialPoint& pt) = 0;
    
    virtual mat3ds DamageCriterionTangentAverage(FEMaterialPoint& pt) = 0;
        
    //! plasticity criterion average
    virtual void PlasticityCriterionAverage(FEMaterialPoint& pt, FEMaterialPoint& rp) = 0;
    
public:
    FENonLocalKernel* m_krnl;   //! kernel function
    bool    m_binit;            //!< initialization flag

public:
    std::vector<std::vector<int>>    m_EPL; //!< list of element proximity lists
    FEMeshTopo      m_topo;                 //!< mesh topology;
    
    FECORE_BASE_CLASS(FENonLocalAveraging);
};

//-----------------------------------------------------------------------------
// Use Bazant non-local averaging scheme
class FEBIOMECH_API FENLABazant : public FENonLocalAveraging
{
public:
    FENLABazant(FEModel* pfem) : FENonLocalAveraging(pfem) {}
    
    //! criterion average
    double DamageCriterionAverage(FEMaterialPoint& pt) override;
    
    //! criterion stress tangent average
    mat3ds DamageCriterionTangentAverage(FEMaterialPoint& pt) override;

    //! plasticity criterion average
    void PlasticityCriterionAverage(FEMaterialPoint& pt, FEMaterialPoint& rp) override {};
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Use Borino-Polizzotti non-local averaging scheme
class FEBIOMECH_API FENLABorino : public FENonLocalAveraging
{
public:
    FENLABorino(FEModel* pfem) : FENonLocalAveraging(pfem) {}
    
    //! criterion average
    double DamageCriterionAverage(FEMaterialPoint& pt) override;
    
    //! criterion stress tangent average
    mat3ds DamageCriterionTangentAverage(FEMaterialPoint& pt) override;
    
    //! plasticity criterion average
    void PlasticityCriterionAverage(FEMaterialPoint& pt, FEMaterialPoint& rp) override {};

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Average over the element
class FEBIOMECH_API FENLAElement : public FENonLocalAveraging
{
public:
    FENLAElement(FEModel* pfem) : FENonLocalAveraging(pfem) {}
    
    bool Init() override;
    
    //! criterion average
    double DamageCriterionAverage(FEMaterialPoint& pt) override;
    
    //! criterion stress tangent average
    mat3ds DamageCriterionTangentAverage(FEMaterialPoint& pt) override;
    
    //! plasticity criterion average
    void PlasticityCriterionAverage(FEMaterialPoint& pt, FEMaterialPoint& rp) override;

};
