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
#include "FEERDDomainFactory.h"
#include <FEBioMix/FEBiphasic.h>
#include <FEBioMix/FEBiphasicSolute.h>
#include <FEBioMix/FETriphasic.h>
#include <FEBioMix/FEMultiphasic.h>
#include "FEElasticReactionDiffusion.h"
#include <FEBioMix/FEBiphasicSolidDomain.h>
#include <FEBioMix/FEBiphasicSoluteDomain.h>
#include <FEBioMix/FETriphasicDomain.h>
#include <FEBioMix/FEMultiphasicDomain.h>
#include "FEElasticReactionDiffusionDomain.h"
#include <FECore/FEShellDomain.h>

//-----------------------------------------------------------------------------
FEDomain* FEERDDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;

	FEDomain* pd = nullptr;

	if (eclass == FE_ELEM_SOLID)
	{
		const char* sztype = 0;
		if (dynamic_cast<FEElasticReactionDiffusion*>(pmat)) sztype = "elastic-reaction-diffusion-solid";

		if (sztype) pd = fecore_new<FESolidDomain>(sztype, pfem);
	}

	if (pd) pd->SetMaterial(pmat);
	return pd;
}
