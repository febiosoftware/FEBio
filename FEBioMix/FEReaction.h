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
#include <FECore/FEMaterial.h>
#include "febiomix_api.h"
#include <map>

//-----------------------------------------------------------------------------
class FESoluteInterface;

//-----------------------------------------------------------------------------
//! Base class for reactions.

typedef std::map<int,int> intmap;
typedef std::map<int,int>::iterator itrmap;

class FEBIOMIX_API FEReaction : public FEMaterialProperty
{
public:
    //! constructor
    FEReaction(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
public:
    //! set stoichiometric coefficient
    void SetStoichiometricCoefficient(intmap& RP, int id, int v) { RP.insert(std::pair<int, int>(id, v)); }

public: //TODO: Make this protected again
    FESoluteInterface* m_psm;   //!< solute interface to parent class
};


//-----------------------------------------------------------------------------
class FEBIOMIX_API FEReactionSpeciesRef : public FEMaterialProperty
{
public:
    enum SpeciesType { UnknownSpecies, SoluteSpecies, SBMSpecies };

public:
    FEReactionSpeciesRef(FEModel* fem);

    bool Init() override;

    int GetSpeciesType() const;

    bool IsSolute() const;
    bool IsSBM() const;

public:
    int     m_speciesID;        // the species ID
    int     m_v;                // stoichiometric coefficient

    // these parameters are mostly for parsing older files that used "sol" or "sbm"
    int     m_solId;
    int     m_sbmId;

private:
    int     m_speciesType;  // solute or sbm?

    DECLARE_FECORE_CLASS();
};

class FEBIOMIX_API FEReactantSpeciesRef : public FEReactionSpeciesRef
{
public: FEReactantSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEReactantSpeciesRef)
};

class FEBIOMIX_API FEProductSpeciesRef : public FEReactionSpeciesRef {
public: FEProductSpeciesRef(FEModel* fem) : FEReactionSpeciesRef(fem) {}
      DECLARE_FECORE_CLASS();
      FECORE_BASE_CLASS(FEProductSpeciesRef)
};
