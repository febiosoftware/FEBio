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



#include "stdafx.h"
#include "FESolutesMaterialPoint.h"
#include "FECore/DumpStream.h"
using namespace std;

//=============================================================================
//   FESolutesMaterialPoint
//=============================================================================


//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FESolutesMaterialPoint::Copy()
{
	FESolutesMaterialPoint* pt = new FESolutesMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FESolutesMaterialPoint::Init()
{
	m_nsol = m_nsbm = 0;
	m_psi = m_cF = 0;
	m_Ie = vec3d(0,0,0);
	m_rhor = 0;
    m_c.clear();
    m_gradc.clear();
    m_j.clear();
    m_ca.clear();
    m_crp.clear();
    m_sbmr.clear();
    m_sbmrp.clear();
    m_sbmrhat.clear();
    m_sbmrhatp.clear();
    m_sbmrmin.clear();
    m_sbmrmax.clear();
    m_k.clear();
    m_dkdJ.clear();
    m_dkdJJ.clear();
    m_dkdc.clear();
    m_dkdJc.clear();
    m_dkdcc.clear();
    m_dkdr.clear();
    m_dkdJr.clear();
    m_dkdrc.clear();
    m_cri.clear();
    m_crd.clear();
    m_strain = 0;
    m_pe = m_pi = 0;
    m_ce.clear();
    m_ci.clear();
    m_ide.clear();
    m_idi.clear();
    m_bsb.clear();
    
	// don't forget to initialize the base class
	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Serialize material point data to the archive
void FESolutesMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_nsol & m_psi & m_cF & m_Ie & m_nsbm;
	ar & m_c & m_gradc & m_j & m_ca & m_crp & m_k & m_dkdJ;
	ar & m_dkdc;
	ar & m_sbmr & m_sbmrp & m_sbmrhat & m_sbmrhatp;
	ar & m_cri;
	ar & m_crd;
	ar & m_strain & m_pe & m_pi;
	ar & m_ce & m_ide;
	ar & m_ci & m_idi;
    ar & m_bsb;
}

//-----------------------------------------------------------------------------
double FESolutesMaterialPoint::Osmolarity() const
{
    double ew = 0.0;
    for (int isol = 0; isol < (int)m_ca.size(); ++isol)
    {
        // exclude solid-bound 'solutes'
        if (!m_bsb[isol]) ew += m_ca[isol];
    }
    return ew;
}
