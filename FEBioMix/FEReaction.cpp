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
#include "FEReaction.h"
#include <FECore/log.h>
#include "FESoluteInterface.h"
#include "FESolute.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEReaction::FEReaction(FEModel* pfem) : FEMaterialProperty(pfem)
{
    m_psm = nullptr;
}

//-----------------------------------------------------------------------------
bool FEReaction::Init()
{
    // make sure the parent class is set
    m_psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
    assert(m_psm);
    if (m_psm == 0) {
        feLogError("Parent class not set or of incorrect type");
        return false;
    }
    
    // now call base class
    return FEMaterialProperty::Init();
}


//=============================================================================
BEGIN_FECORE_CLASS(FEReactionSpeciesRef, FEMaterialProperty)
	ADD_PARAMETER(m_speciesID, "species", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_solId    , "sol")->SetFlags(FE_PARAM_ATTRIBUTE | FE_PARAM_HIDDEN);
	ADD_PARAMETER(m_sbmId    , "sbm")->SetFlags(FE_PARAM_ATTRIBUTE | FE_PARAM_HIDDEN);
END_FECORE_CLASS();

FEReactionSpeciesRef::FEReactionSpeciesRef(FEModel* fem) : FEMaterialProperty(fem) 
{
    m_speciesID = -1;
    m_solId = -1;
    m_sbmId = -1;
    
    m_v = 0;

    m_speciesType = UnknownSpecies;
}

int FEReactionSpeciesRef::GetSpeciesType() const
{
    assert(m_speciesType != UnknownSpecies);
    return m_speciesType;
}

bool FEReactionSpeciesRef::IsSolute() const
{
    return (m_speciesType == SoluteSpecies);
}

bool FEReactionSpeciesRef::IsSBM() const
{
    return (m_speciesType == SBMSpecies);
}

bool FEReactionSpeciesRef::Init()
{
    // make sure the species ID is valid
    if ((m_speciesID == -1) && (m_solId == -1) && (m_sbmId == -1)) return false;

    // figure out if this is a solute or sbm
    if (m_speciesType == UnknownSpecies)
    {
        FEModel& fem = *GetFEModel();
        if (m_speciesID != -1)
        {
            int id = m_speciesID - 1;
            int n = 0;
            for (int i = 0; i < fem.GlobalDataItems(); ++i)
            {
                FESoluteData* sol = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
                FESBMData* sbm = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
                if (sol || sbm)
                {
                    if (id == n)
                    {
                        if (sol) { m_speciesType = SoluteSpecies; m_speciesID = sol->GetID(); }
                        if (sbm) { m_speciesType = SBMSpecies; m_speciesID = sbm->GetID(); }

                        break;
                    }
                    else n++;
                }
            }
        }
        else if (m_solId != -1)
        {
            for (int i = 0; i < fem.GlobalDataItems(); ++i)
            {
                FESoluteData* sol = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
                if (sol && (sol->GetID() == m_solId))
                {
                    m_speciesType = SoluteSpecies;
                    m_speciesID = sol->GetID(); 
                }
            }
        }
        else if (m_sbmId != -1)
        {
            for (int i = 0; i < fem.GlobalDataItems(); ++i)
            {
                FESBMData* sbm = dynamic_cast<FESBMData*>(fem.GetGlobalData(i));
                if (sbm && (sbm->GetID() == m_sbmId))
                {
                    m_speciesType = SBMSpecies;
                    m_speciesID = sbm->GetID();
                }
            }
        }

        assert(m_speciesType != UnknownSpecies);
        if (m_speciesType == UnknownSpecies) return false;
    }

    return FEMaterialProperty::Init();
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEReactantSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vR");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEProductSpeciesRef, FEReactionSpeciesRef)
    ADD_PARAMETER(m_v, "vP");
END_FECORE_CLASS();
