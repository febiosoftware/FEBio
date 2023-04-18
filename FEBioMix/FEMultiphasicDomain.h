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
#include "FECore/FESolidDomain.h"
#include "FEMultiphasic.h"
#include "FEBioMech/FEElasticDomain.h"

//-----------------------------------------------------------------------------
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FESolver;

//-----------------------------------------------------------------------------
//! Abstract interface class for multiphasic domains.

//! A biphasic domain is used by the biphasic solver.
//! This interface defines the functions that have to be implemented by a
//! biphasic domain. There are basically two categories: residual functions
//! that contribute to the global residual vector. And stiffness matrix
//! function that calculate contributions to the global stiffness matrix.
class FEBIOMIX_API FEMultiphasicDomain : public FEElasticDomain
{
public:
    FEMultiphasicDomain(FEModel* pfem);
    virtual ~FEMultiphasicDomain(){}
    
    // --- R E S I D U A L ---
    
    //! internal work for steady-state case
    virtual void InternalForcesSS(FEGlobalVector& R) = 0;
    
    // --- S T I F F N E S S   M A T R I X ---
    
    //! calculates the global stiffness matrix for this domain
    virtual void StiffnessMatrix(FELinearSystem& LS, bool bsymm) = 0;
    
    //! calculates the global stiffness matrix (steady-state case)
    virtual void StiffnessMatrixSS(FELinearSystem& LS, bool bsymm) = 0;
    
protected:
    FEMultiphasic*      m_pMat;
    int                 m_dofP;		//!< pressure dof index
    int                 m_dofQ;     //!< shell extra pressure dof index
    int                 m_dofC;		//!< concentration dof index
    int                 m_dofD;		//!< shell extra concentration dof index
    int                 m_dofVX;
    int                 m_dofVY;
    int                 m_dofVZ;
    
protected:
    bool                m_breset;   //! flag for calling reset
};
