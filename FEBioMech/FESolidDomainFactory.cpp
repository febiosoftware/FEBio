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
#include "FESolidDomainFactory.h"
#include "FERigidMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FELinearTrussDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FERemodelingElasticDomain.h"
#include "FEUDGHexDomain.h"
#include "FEUT4Domain.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEDiscreteElasticDomain.h"
#include "FERemodelingElasticMaterial.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FEDiscreteElementMaterial.h"
#include "FESRIElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FESolidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	FE_Element_Shape eshape = spec.eshape;
	FE_Element_Type etype = spec.etype;

	// this will store the domain we are going to allocate
	FEDomain* pd = nullptr;

	// try to allocate the domain
	if (eclass == FE_ELEM_SOLID)
	{
		if (dynamic_cast<FESolidMaterial*>(pmat) == nullptr) return nullptr;

		const char* sztype = nullptr;
		if      (dynamic_cast<FERigidMaterial*            >(pmat)) sztype = "rigid-solid";
		else if (dynamic_cast<FERemodelingElasticMaterial*>(pmat)) sztype = "remodeling-solid";
		else if (eshape == ET_HEX8)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else
			{
				if (etype == FE_HEX8G1) sztype = "udg-hex";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_TET10) || (eshape == ET_TET15) || (eshape == ET_TET20))
		{
			if ((etype == FE_TET10G8RI4) || (etype == FE_TET10G4RI1))
			{
				sztype = "sri-solid";
			}
			else
			{
				if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_HEX20) || (eshape == ET_HEX27) || (eshape == ET_PYRA13))
		{
//			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			sztype = "elastic-solid";
		}
		else if (eshape == ET_TET4)
		{
			if (spec.m_but4) sztype = "ut4-solid";
			else sztype = "elastic-solid";
		}
		else if (eshape == ET_TET5)
		{
			sztype = "elastic-solid";
		}
		else if (eshape == ET_PENTA6)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else sztype = "elastic-solid";
		}
		else if (eshape == ET_PENTA15)
		{
			// three-field implementation for uncoupled materials
//          if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			sztype = "elastic-solid";
		}
		else if (eshape == ET_PYRA5)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else sztype = "elastic-solid";
		}

		if (sztype) pd = fecore_new<FESolidDomain>(sztype, pfem);
	}
	else if (eclass == FE_ELEM_SHELL)
	{
		if (dynamic_cast<FESolidMaterial*>(pmat) == nullptr) return nullptr;

		const char* sztype = nullptr;
		if (dynamic_cast<FERigidMaterial*>(pmat))
		{
			if (spec.m_shell_formulation == OLD_SHELL) sztype = "rigid-shell-old";
			else sztype = "rigid-shell";
		}
		else if (dynamic_cast<FESolidMaterial*>(pmat))
		{
			if ((eshape == ET_QUAD4) || (eshape == ET_TRI3) || (eshape == ET_QUAD8) || (eshape == ET_TRI6))
			{
				switch (spec.m_shell_formulation)
				{
				case NEW_SHELL:
					// three-field implementation for uncoupled materials
					if (dynamic_cast<FEUncoupledMaterial*>(pmat)) {
						if (spec.m_bthree_field) sztype = "three-field-shell";
						else if ((eshape == ET_QUAD4) && spec.m_bthree_field) sztype = "three-field-shell";
						else if ((eshape == ET_TRI3) && spec.m_bthree_field) sztype = "three-field-shell";
						else sztype = "elastic-shell";
					}
					else sztype = "elastic-shell";
					break;
				case OLD_SHELL: sztype = "elastic-shell-old"; break;
				case EAS_SHELL: sztype = "elastic-shell-eas"; break;
				case ANS_SHELL: sztype = "elastic-shell-ans"; break;
				default:
					return 0;
				}
			}
		}

		if (sztype) pd = fecore_new<FEShellDomain>(sztype, pfem);
	}
	else if (eclass == FE_ELEM_TRUSS)
	{
		const char* sztype = nullptr;
		if (dynamic_cast<FETrussMaterial*>(pmat))
			if (eshape == ET_TRUSS2) sztype = "linear-truss";

		if (sztype) pd = fecore_new<FETrussDomain>(sztype, pfem);
	}
	else if (eclass == FE_ELEM_DISCRETE)
	{
		const char* sztype = nullptr;
		if (dynamic_cast<FEDiscreteMaterial*>(pmat))
		{
	//		if      (eclass == FE_ELEM_WIRE    ) sztype = "deformable-spring";
			if      (eclass == FE_ELEM_WIRE    ) sztype = "deformable-spring2";
			else if (eclass == FE_ELEM_DISCRETE)
			{
				if (dynamic_cast<FEDiscreteElasticMaterial*>(pmat)) sztype = "discrete";
			}
			else return 0;
		}

		if (sztype) pd = fecore_new<FEDiscreteDomain>(sztype, pfem);
	}

	if (pd) pd->SetMaterial(pmat);
	return pd;
}
