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
#include <FECore/FEPrescribedBC.h>
#include <FECore/FEModelParam.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Prescribe the actual fluid pressure in a multiphasic mixture
//!
class FEBIOMIX_API FEMultiphasicFluidPressureBC : public FEPrescribedSurface
{
public:
    //! constructor
    FEMultiphasicFluidPressureBC(FEModel* pfem);

    //! set the dilatation
    void Update() override;

    //! initialize
    bool Init() override;

    //! serialization
    void Serialize(DumpStream& ar) override;

public:
    // return the value for node i, dof j
    void GetNodalValues(int nodelid, std::vector<double>& val) override;

    // copy data from another class
    void CopyFrom(FEBoundaryCondition* pbc) override;

private:
    FEParamDouble   m_p;        //!< prescribed fluid pressure
    vector<double>  m_pe;       //!< effective fluid pressure

private:
    double      m_Rgas;
    double      m_Tabs;

    int         m_dofP;
    int         m_dofC;

public:
    bool        m_bshellb;      //!< flag for prescribing pressure on shell bottom

    DECLARE_FECORE_CLASS();
};
