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
#include "FEReactionERD.h"
#include <FECore/log.h>
#include "FEElasticReactionDiffusionInterface.h"
#include <FEBioMix/FESolute.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEReactionERD::FEReactionERD(FEModel* pfem) : FEMaterialProperty(pfem)
{
    m_psm = nullptr;
}

//-----------------------------------------------------------------------------
bool FEReactionERD::Init()
{
    // make sure the parent class is set
    m_psm = dynamic_cast<FEElasticReactionDiffusionInterface*>(GetAncestor());
    assert(m_psm);
    if (m_psm == 0) {
        feLogError("Parent class not set or of incorrect type");
        return false;
    }

    // now call base class
    return FEMaterialProperty::Init();
}


//=============================================================================
BEGIN_FECORE_CLASS(FEReactionSpeciesRefERD, FEMaterialProperty)
ADD_PARAMETER(m_speciesID, "species", FE_PARAM_ATTRIBUTE, "$(species)");
ADD_PARAMETER(m_solId, "sol")->SetFlags(FE_PARAM_ATTRIBUTE | FE_PARAM_HIDDEN);
END_FECORE_CLASS();

FEReactionSpeciesRefERD::FEReactionSpeciesRefERD(FEModel* fem) : FEMaterialProperty(fem)
{
    m_speciesID = -1;
    m_solId = -1;

    m_v = 0;
}

bool FEReactionSpeciesRefERD::Init()
{
    // make sure the species ID is valid
    if ((m_speciesID == -1) && (m_solId == -1)) return false;

    // figure out if this is a solute or sbm
    FEModel& fem = *GetFEModel();
    if (m_speciesID != -1)
    {
        int id = m_speciesID - 1;
        int n = 0;
        for (int i = 0; i < fem.GlobalDataItems(); ++i)
        {
            FESoluteData* sol = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
            if (sol)
            {
                if (id == n)
                {
                    if (sol) { m_speciesID = sol->GetID(); }

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
                m_speciesID = sol->GetID();
            }
        }
    }

    return FEMaterialProperty::Init();
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEReactantSpeciesRefERD, FEReactionSpeciesRefERD)
ADD_PARAMETER(m_v, "vR");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEProductSpeciesRefERD, FEReactionSpeciesRefERD)
ADD_PARAMETER(m_v, "vP");
END_FECORE_CLASS();
