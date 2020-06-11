/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEBiphasic.h"
#include <map>

typedef std::map<int,double> idmap;     //!< map integer id with double value
typedef std::map<int,double>::iterator itridmap;

//-----------------------------------------------------------------------------
// This class implements a material that has a solvent supply following
// Starling's equation

class FEBIOMIX_API FESolventSupplyStarling :	public FESolventSupply
{
public:
	//! constructor
	FESolventSupplyStarling(FEModel* pfem);
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt) override;
	
	//! Tangent of supply with respect to strain
	mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to pressure
	double Tangent_Supply_Pressure(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp, const int isol);
	
    //! set parameter attribute for indexed solute parameters
	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval) override;
    
	//! set value of indexed parameters
	void SetIndexedParameter(idmap& p, int id, double val) { p.insert(std::pair<int, double>(id, val)); }
    
public:
	double		m_kp;				//!< coefficient of pressure drop
	double		m_pv;				//!< prescribed (e.g., vascular) pressure
	vector<double>		m_qc;       //!< coefficients of concentration drops
	vector<double>		m_cv;       //!< prescribed (e.g., vascular) concentrations
    double  m_qctmp;                //!< helper variable for reading in m_qc
    idmap	m_qcinp;                //!< m_qc for each solute (input)
    double  m_cvtmp;                //!< helper variable for reading in m_cv
    idmap	m_cvinp;                //!< m_cv for each solute (input)
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
