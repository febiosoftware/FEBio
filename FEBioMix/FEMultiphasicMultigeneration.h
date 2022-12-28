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
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
// forward declaration of material class
class FEMultiphasicMultigeneration;

//-----------------------------------------------------------------------------
//! Multigenerational SBM material point.
//! This material point stores the inverse of the relative deformation gradient,
//! and the increment in the mass of solid-bound molecular species in
//! multiple generations.

class FEBIOMIX_API FEMultigenSBMMaterialPoint : public FEMaterialPointData
{
public:
	FEMultigenSBMMaterialPoint(FEMultiphasicMultigeneration* pm, FEMaterialPointData* pt) : m_pmat(pm), FEMaterialPointData(pt) { m_tgen = 0.0; }
    
	FEMaterialPointData* Copy();
    
	void Serialize(DumpStream& ar);
    
	void Init();
	void Update(const FETimeInfo& timeInfo);
    
public:
	// multigenerational material data
    int                         m_ngen;     //!< number of generations
	vector <mat3d>              m_Fi;       //!< inverse of relative deformation gradient
	vector <double>             m_Ji;       //!< determinant of Fi (store for efficiency)
    int                         m_nsbm;     //!< number of solid-bound molecules
    vector< vector<double> >    m_gsbmr;    //!< sbmr content at each generation
    vector< vector<double> >    m_gsbmrp;   //!< gsbmr at previous time point
    vector<double>              m_lsbmr;    //!< last generation sbmr values
	double                      m_tgen;     //!< last generation time
    
private:
	FEMultiphasicMultigeneration*	m_pmat;
};

//-----------------------------------------------------------------------------
//! Multigeneration multiphasic material.

class FEMultiphasicMultigeneration : public FEMultiphasic
{
public:
	//! constructor
	FEMultiphasicMultigeneration(FEModel* pfem);
    
    //! returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
	
    //! Update solid bound molecules
    void UpdateSolidBoundMolecules(FEMaterialPoint& mp) override;

	int CheckGeneration(const double t);
    double GetGenerationTime(const int igen);
    
public:
	double	m_gtime;	//!< time duration of each generation
    
    DECLARE_FECORE_CLASS();

};
