/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! This surface load represents the traction applied on the solid at the
//! interface between a fluid and solid in an FSI analysis.
class FEBIOFLUID_API FEFluidFSITraction : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidFSITraction(FEModel* pfem);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! calculate residual
    void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! initialization
    bool Init() override;
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const FETimeInfo& tp, const int iel);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const FETimeInfo& tp, const int iel);
    
protected:
	vector<double>      m_K;        //!< fluid bulk modulus
	vector<double>      m_s;        //!< scale factor
	vector<bool>        m_bself;    //!< flag if fluid pressure is applied on its own FSI mesh
	vector<FEElement*>  m_elem;     //!< list of fluid-FSI elements

    // degrees of freedom
    int		m_dofX, m_dofY, m_dofZ;
    int     m_dofSX, m_dofSY, m_dofSZ;
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofEF, m_dofEFP;
};
