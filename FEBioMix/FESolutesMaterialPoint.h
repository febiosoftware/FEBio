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
#include <FECore/FEMaterialPoint.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Class for storing material point data for solute materials

class FEBIOMIX_API FESolutesMaterialPoint : public FEMaterialPoint
{
public:
	//! Constructor
	FESolutesMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	//! Create a shallow copy
	FEMaterialPoint* Copy();
	
	//! serialize data
	void Serialize(DumpStream& ar);
    
	//! Initialize material point data
	void Init();
	
public:
	// solutes material data
	int				m_nsol;		//!< number of solutes
	vector<double>	m_c;		//!< effective solute concentration
	vector<vec3d>	m_gradc;	//!< spatial gradient of solute concentration
	vector<vec3d>	m_j;		//!< solute molar flux
	vector<double>	m_ca;		//!< actual solute concentration
    vector<double>  m_crp;      //!< referential actual solute concentration at previous time step
	double			m_psi;		//!< electric potential
	vec3d			m_Ie;		//!< current density
	double			m_cF;		//!< fixed charge density in current configuration
	int				m_nsbm;		//!< number of solid-bound molecules
	double			m_rhor;		//!< current referential mass density
	vector<double>	m_sbmr;		//!< referential mass concentration of solid-bound molecules
	vector<double>	m_sbmrp;	//!< m_sbmr at previoust time step
	vector<double>	m_sbmrhat;	//!< referential mass supply of solid-bound molecules
    vector<double>  m_sbmrhatp; //!< referential mass supply of solid-bound molecules at previous time step
	vector<double>	m_sbmrmin;	//!< minimum value of m_sbmr
	vector<double>	m_sbmrmax;	//!< maximum value of m_sbmr
	vector<double>	m_k;		//!< solute partition coefficient
	vector<double>	m_dkdJ;		//!< 1st deriv of m_k with strain (J)
	vector<double>	m_dkdJJ;	//!< 2nd deriv of m_k with strain (J)
	vector< vector<double> >	m_dkdc;			//!< 1st deriv of m_k with effective concentration
	vector< vector<double> >	m_dkdJc;		//!< cross deriv of m_k with J and c
	vector< vector< vector<double> > > m_dkdcc;	// 2nd deriv of m_k with c
	vector< vector<double> >	m_dkdr;			//!< 1st deriv of m_k with m_sbmr
	vector< vector<double> >	m_dkdJr;		//!< cross deriv of m_k with J and m_sbmr
	vector< vector< vector<double> > > m_dkdrc;	//!< cross deriv of m_k with m_sbmr and c
    vector<int>     m_cri;      //!< optional integer data needed for chemical reactions
    vector<double>  m_crd;      //!< optional double data needed for chemical reactions
    double          m_strain;   //!< areal strain
    double          m_pe;       //!< effective fluid pressure on external side
    double          m_pi;       //!< effective fluid pressure on internal side
    vector<double>  m_ce;       //!< effective solute concentration on external side
    vector<double>  m_ci;       //!< effective solute concentration on internal side
    vector<int>     m_ide;      //!< solute IDs on external side
    vector<int>     m_idi;      //!< solute IDs on internal side
};

