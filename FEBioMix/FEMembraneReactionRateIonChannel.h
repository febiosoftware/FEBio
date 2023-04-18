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
#include "FEMembraneReaction.h"

class FEBIOMIX_API FEMembraneReactionRateIonChannel : public FEMembraneReactionRate
{
public:
    //! constructor
    FEMembraneReactionRateIonChannel(FEModel* pfem);
    
    // initialization
    bool Init() override;
    
    //! reaction rate at material point
    double ReactionRate(FEMaterialPoint& pt) override;
    
    //! tangent of reaction rate with strain at material point
    double Tangent_ReactionRate_Strain(FEMaterialPoint& pt) override;
    
    //! tangent of reaction rate with effective fluid pressure at material point
    double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) override {return 0; }
    double Tangent_ReactionRate_Pe(FEMaterialPoint& pt) override { return 0; }
    double Tangent_ReactionRate_Pi(FEMaterialPoint& pt) override { return 0; }
    
    //! tangent of reaction rate with effective solute concentration at material point
    double Tangent_ReactionRate_Concentration(FEMaterialPoint& pt, const int isol) override {return 0; }
    double Tangent_ReactionRate_Ce(FEMaterialPoint& pt, const int isol) override;
    double Tangent_ReactionRate_Ci(FEMaterialPoint& pt, const int isol) override;
    
public:
    int     m_sol;      //!< solute id (1-based) for ion
    int     m_lid;      //!< local id of solute (zero-based)
    int     m_z;        //!< charge number of channel ion
    int     m_sbm;      //!< sbm id (1-based) for channel protein
    double  m_g;        //!< channel conductance
    
    
    DECLARE_FECORE_CLASS();
};
