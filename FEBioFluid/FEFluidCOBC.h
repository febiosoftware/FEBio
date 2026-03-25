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
#include <FECore/FEFunction1D.h>
#include "FEFluidMaterial.h"

//-----------------------------------------------------------------------------
//! FEFluidCOBC is a fluid surface load that implements a coronary outflow boundary condition
//!
class FEBIOFLUID_API FEFluidCOBC : public FEPrescribedSurface
{
public:
    //! constructor
    FEFluidCOBC(FEModel* pfem);
    
    //! set the dilatation
    void Update() override;
    void UpdateModel() override;
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;

public:
	void PrepStep(std::vector<double>& ui, bool brel) override;

    // return the value for node i, dof j
    void GetNodalValues(int nodelid, std::vector<double>& val) override;

    // copy data from another class
    void CopyFrom(FEBoundaryCondition* pbc) override;

private:
	//! set the dilatation
	void UpdateDilatation();

private:
    double          m_Ra;       //!< arterial resistance
    double          m_Ca;       //!< arterial compliance
    double          m_Ram;      //!< microvascular arterial resistance
    double          m_Cim;      //!< myocardial compliance
    double          m_Rv;       //!< ventricular resistance
    double          m_Rvm;      //!< microvascular ventricular resistance
    double          m_p0;       //!< initial fluid pressure
    FEFunction1D*   m_Pmi;      //!< myocardial pressure
    FEFunction1D*   m_PRA;      //!< right atrium pressure
    double          m_pd;       //!< downstream pressure
    
private:
    double              m_pn;   //!< fluid pressure at current time point
    double              m_pp;   //!< fluid pressure at previous time point
    double              m_ppp;  //!< fluid pressure at 2nd-previous time point
    double              m_q;    //!< flow rate at current time point
    double              m_qp;   //!< flow rate at previous time point
    double              m_qpp;  //!< flow rate at 2nd-previous time point
    double              m_pdn;  //!< downstream fluid pressure at current time point
    double              m_pdp;  //!< downstream fluid pressure at previous time point
    double              m_tp;   //!< previous time
    double              m_tpp;  //!< 2nd-previous time
    double              m_e;
    FEFluidMaterial*    m_pfluid;   //!< pointer to fluid
    FESurface*          m_psurf;    //!< pointer to surface
    
    FEDofList   m_dofW;
    int         m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
