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
#include "FEElasticShellDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! shell element. Results indicate that using this class produces poorer convergence
//! with shells than the standard FEElasticShellDomain.  This class is included
//! only for development purposes.
class FEBIOMECH_API FE3FieldElasticShellDomain : public FEElasticShellDomain
{
protected:
    struct ELEM_DATA
    {
        double    eJ;        // average element jacobian
        double    ep;        // average pressure
        double    Lk;        // Lagrangian multiplier

		void Serialize(DumpStream& ar);
    };
    
public:
    //! constructor
	FE3FieldElasticShellDomain(FEModel* pfem);
    
    //! \todo Is this really used?
	FE3FieldElasticShellDomain& operator = (FE3FieldElasticShellDomain& d);
    
    //! initialize class
	bool Init() override;
    
    //! Reset data
    void Reset() override;
    
    //! augmentation
    bool Augment(int naug) override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
public: // overridden from FEElasticDomain
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // calculate stiffness matrix
    void StiffnessMatrix(FELinearSystem& LS) override;
    
protected:
    //! Dilatational stiffness component for nearly-incompressible materials
    void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke);
    
    //! material and geometrical stiffness components
    void ElementStiffness(int iel, matrix& ke);
    
    //! update the stress of an element
    void UpdateElementStress(int iel);

public:
	bool DoAugmentations() const;
    
protected:
    vector<ELEM_DATA>    m_Data;

	bool	m_blaugon;		//!< augmented lagrangian flag
	double	m_augtol;		//!< augmented lagrangian tolerance
	int		m_naugmin;		//!< minimum number of augmentations
	int		m_naugmax;		//!< max number of augmentations

	DECLARE_FECORE_CLASS();
};
