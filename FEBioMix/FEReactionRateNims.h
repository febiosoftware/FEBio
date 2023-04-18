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
#include "FEChemicalReaction.h"

//-----------------------------------------------------------------------------
//! Concentration-history-dependent reaction rate.
//! Reaction rate depends on concentration of a solute (e.g., growth factor)
//! and whether solute has been released (removed) at some release time.
//! Before release, reaction rate varies linearly with history of maximum solute
//! concentration cmax, from k0 at cmax=0 to kc at cmax=cc, then holds constant at kc.
//! After release, reaction rate increases from k0 at cmax=0 to kr at cmax=cr,
//! then holds constant at kr. Release time is trel.

class FEBIOMIX_API FEReactionRateNims : public FEReactionRate
{
public:
	//! constructor
	FEReactionRateNims(FEModel* pfem);
	
	//! data initialization and checking
	bool Init() override;
	
	//! reaction rate at material point
	double ReactionRate(FEMaterialPoint& pt) override;
	
	//! tangent of reaction rate with strain at material point
	mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override;
	
	//! tangent of reaction rate with effective fluid pressure at material point
	double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override;
	
    //! reset, initialize and update chemical reaction data in the FESolutesMaterialPoint
    void ResetElementData(FEMaterialPoint& mp) override;
    void InitializeElementData(FEMaterialPoint& mp) override;
    void UpdateElementData(FEMaterialPoint& mp) override;
    
public:
    int             m_sol;                  //!< solute id (1-based)
    int             m_lid;                  //!< local id of solute (zero-based)
    FEParamDouble   m_k0;                   //!< reaction rate at zero concentration
    FEParamDouble   m_cc;                   //!< concentration cc
    FEParamDouble   m_kc;                   //!< reaction rate at cc;
    FEParamDouble   m_cr;                   //!< concentration cr;
    FEParamDouble   m_kr;                   //!< reaction rate at cr;
    FEParamDouble   m_trel;                 //!< release time;
    int             m_cmax;                 //!< index of entry in m_crd
	
	DECLARE_FECORE_CLASS();
};
