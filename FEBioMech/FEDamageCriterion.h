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
#include "FENonLocalAveraging.h"
#include <FECore/FEMaterial.h>
#include <FECore/FEMeshTopo.h>

//-----------------------------------------------------------------------------
// Virtual base class for damage criterion

class FEBIOMECH_API FEDamageCriterion : public FEMaterialProperty
{
public:
    FEDamageCriterion(FEModel* pfem) : FEMaterialProperty(pfem) { m_nla = nullptr; }
    virtual ~FEDamageCriterion() {}
    
	//! damage
    virtual double DamageCriterion(FEMaterialPoint& pt) { if (m_nla) return m_nla->DamageCriterionAverage(pt); else return DCpt(pt); }
    
    //! criterion tangent with respect to stress
    virtual mat3ds CriterionStressTangent(FEMaterialPoint& pt) { if (m_nla) return m_nla->DamageCriterionTangentAverage(pt); else return CSTpt(pt); }
    
    //! damage at single point
    virtual double DCpt(FEMaterialPoint& pt) = 0;

    //! criterion tangent with respect to stress at single pont
    virtual mat3ds CSTpt(FEMaterialPoint& pt) { return mat3ds(0); }
    
public:
    FENonLocalAveraging*    m_nla;  //! optional nonlocal averaging scheme
    
    FECORE_BASE_CLASS(FEDamageCriterion);
};

//-----------------------------------------------------------------------------
// Simo's damage criterion

class FEDamageCriterionSimo : public FEDamageCriterion
{
public:
    FEDamageCriterionSimo(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Strain energy density as damage criterion

class FEDamageCriterionSED : public FEDamageCriterion
{
public:
	FEDamageCriterionSED(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Specific strain energy as damage criterion

class FEDamageCriterionSSE : public FEDamageCriterion
{
public:
    FEDamageCriterionSSE(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// von Mises stress as damage criterion

class FEDamageCriterionVMS : public FEDamageCriterion
{
public:
	FEDamageCriterionVMS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    //! criterion tangent with respect to stress at single pont
    mat3ds CSTpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Drucker yield criterion

class FEDamageCriterionDrucker : public FEDamageCriterion
{
public:
    FEDamageCriterionDrucker(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    //! criterion tangent with respect to stress at single pont
    mat3ds CSTpt(FEMaterialPoint& pt) override;

public:
    FEParamDouble   m_c;    //!< Drucker material parameter
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// max shear stress as damage criterion

class FEDamageCriterionMSS : public FEDamageCriterion
{
public:
	FEDamageCriterionMSS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    //! criterion tangent with respect to stress at single pont
    mat3ds CSTpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// max normal stress as damage criterion

class FEDamageCriterionMNS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    //! criterion tangent with respect to stress at single pont
    mat3ds CSTpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// max normal Lagrange strain as damage criterion

class FEDamageCriterionMNLS : public FEDamageCriterion
{
public:
	FEDamageCriterionMNLS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// max octahedral strain as damage criterion

class FEDamageCriterionOSS : public FEDamageCriterion
{
public:
    FEDamageCriterionOSS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// max octahedral natural strain as damage criterion

class FEDamageCriterionONS : public FEDamageCriterion
{
public:
    FEDamageCriterionONS(FEModel* pfem) : FEDamageCriterion(pfem) {}
    
    //! damage at single point
    double DCpt(FEMaterialPoint& pt) override;

    DECLARE_FECORE_CLASS();
};
