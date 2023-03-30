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
#include "FEFluidMaterial.h"

//-----------------------------------------------------------------------------
//! FEFluidRCBC is a fluid surface that has an RC-equivalent circuit for outflow conditions
//!
class FEBIOFLUID_API FEFluidRCBC : public FEPrescribedSurface
{
public:
    //! constructor
    FEFluidRCBC(FEModel* pfem);
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;

    //! Update 
    void Update() override;
	void UpdateModel() override;

    //! serialization
    void Serialize(DumpStream& ar) override;

public:
	void PrepStep(std::vector<double>& ui, bool brel);

    // return the value for node i, dof j
    void GetNodalValues(int nodelid, std::vector<double>& val) override;

    // copy data from another class
    void CopyFrom(FEBoundaryCondition* pbc) override;

private:
	//! set the dilatation
	void UpdateDilatation();

private:
    double          m_R;        //!< flow resistance
    double          m_p0;       //!< initial fluid pressure
    double          m_C;        //!< capacitance
    bool            m_Bern;     //!< Use Bernoulli's Relation (Q*|Q|)

    double          m_qt;       //!< flow rate at current time step
    double          m_dqt;      //!< flow rate time derivative at current time step
    double          m_pt;       //!< pressure at current time step
    double          m_dpt;      //!< pressure derivative at current time step
    double          m_qp;       //!< flow rate at previous time step
    double          m_dqp;      //!< flow rate time derivative at previous time step
    double          m_pp;       //!< pressure at previoust time step
    double          m_dpp;      //!< pressure derivative at previoust time step
    double          m_e;
    
private:
    FEFluidMaterial*    m_pfluid;   //!< pointer to fluid
    FESurface*          m_psurf;    //!< pointer to surface

    FEDofList   m_dofW;
    int         m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
