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
#include "febiofluid_api.h"
#include <FECore/FEPrescribedBC.h>
#include <FECore/FESurface.h>
#include <FECore/FENodeNodeList.h>

//-----------------------------------------------------------------------------
//! FEThermoFluidTemperatureBC prescribes a temperature BC on an inflow/outflow boundary
//! by enforcing a simplified version of the energy balance as an essential BC
//!
class FEBIOFLUID_API FEThermoFluidTemperatureBC : public FEPrescribedSurface
{
public:
    //! constructor
    FEThermoFluidTemperatureBC(FEModel* pfem);

    //! set the dilatation
    void Update() override;

    //! initialize
    bool Init() override;

    //! serialization
    void Serialize(DumpStream& ar) override;

    //! evaluate flow rate
    void MarkBackFlow();
    
public:
    // return the value for node i, dof j
    void GetNodalValues(int nodelid, std::vector<double>& val) override;

    // copy data from another class
    void CopyFrom(FEBoundaryCondition* pbc) override;

private:
    vector<double>  m_T;
    FESurface*      m_surf;
    FENodeNodeList  m_nnlist;
    vector<bool>    m_backflow;

    FEDofList       m_dofW;
    int             m_dofT;
};
