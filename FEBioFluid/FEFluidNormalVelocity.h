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
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FEFluidNormalVelocity is a fluid surface that has a normal
//! velocity prescribed on it.  This routine simultaneously prescribes a
//! surface load and nodal prescribed velocities
class FEBIOFLUID_API FEFluidNormalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalVelocity(FEModel* pfem);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS) override {}
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R) override;
    
	//! serializatsion
	void Serialize(DumpStream& ar) override;
    
    //! set the velocity
    void Update() override;
    
    //! initialization
    bool Init() override;
    
    //! activate
    void Activate() override;
    
    //! parabolic velocity profile
    bool SetParabolicVelocity();

    //! rim pressure
    bool SetRimPressure();
    
    double ScalarLoad(FESurfaceMaterialPoint& mp) override;

private:
    double NormalVelocity(FESurfaceMaterialPoint& mp);

private:
    FEParamDouble	m_velocity;	//!< average velocity
    bool            m_brim;     //!< flag for setting rim pressure

    // obsolete parameters
    bool            m_bpv;      //!< flag for prescribing nodal values
    bool            m_bpar;

private:
    vector<vec3d>   m_nu;       //!< nodal normals
    vector<double>  m_VN;       //!< nodal values for velocity
    vector<int>     m_rim;      //!< list of nodes on the rim

	FEDofList	m_dofW;
    FEDofList   m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
