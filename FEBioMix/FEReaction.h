//
//  FEReaction.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/4/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEReaction_hpp
#define FEReaction_hpp

#include "FECore/FEMaterial.h"
#include "FEBioMix/FESolutesMaterialPoint.h"
#include <map>

//-----------------------------------------------------------------------------
class FEMultiphasic;

//-----------------------------------------------------------------------------
//! Base class for reactions.

typedef std::map<int,int> intmap;
typedef std::map<int,int>::iterator itrmap;

class FEReaction : public FEMaterial
{
public:
    //! constructor
    FEReaction(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
public:
    //! set stoichiometric coefficient
    void SetStoichiometricCoefficient(intmap& RP, int id, int v) { RP.insert(std::pair<int, int>(id, v)); }
    
public:
    FEMultiphasic*    m_pMP;        //!< pointer to multiphasic material where reaction occurs
};

#endif /* FEReaction_hpp */
