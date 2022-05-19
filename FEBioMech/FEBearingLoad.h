/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEModelParam.h>
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! The bearing load is a surface domain that sustains a parabolic or sinusoidal pressure boundary
//! condition that simulates the contact pressure in a bearing.
//!
class FEBearingLoad : public FESurfaceLoad
{
public:
    enum P_PROFILE {
        P_SINE,         // sinusoidal pressure distribution
        P_PARA          // parabolic pressure distribution
    };

public:
    //! constructor
    FEBearingLoad(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    //! update
    void Update() override;
    
    //! evaluate bearing pressure
    double ScalarLoad(FESurfaceMaterialPoint& mp) override;

    //! serialization
    void Serialize(DumpStream& ar) override;

public:
    //! calculate residual
    void LoadVector(FEGlobalVector& R) override;
    
    //! calculate stiffness
    void StiffnessMatrix(FELinearSystem& LS) override;
    
protected:
    double          m_scale;        //!< scale factor for bearing load
    vec3d           m_force;        //!< bearing load
    bool            m_bsymm;        //!< use symmetric formulation
    bool            m_blinear;      //!< is the load linear (i.e. it will be calculated in the reference frame and assummed deformation independent)
    bool            m_bshellb;      //!< flag for prescribing pressure on shell bottom
    int             m_profile;      //!< profile of pressure distribution

private:
    vec3d           m_er;           //!< unit vector along radial direction of bearing surface
    FESurfaceMap*   m_pc;            //!< pressure value

    DECLARE_FECORE_CLASS();
};
