#pragma once
#include "FECore/FEMaterial.h"
#include "FEBioMix/FESolutesMaterialPoint.h"
#include <map>

//-----------------------------------------------------------------------------
class FEMultiphasic;

//-----------------------------------------------------------------------------
//! Base class for reactions.

typedef std::map<int,int> intmap;
typedef std::map<int,int>::iterator itrmap;

class FEBIOMIX_API FEReaction : public FEMaterial
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
