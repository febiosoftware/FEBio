#include "stdafx.h"
#include "FESolidDomainFactory.h"
#include "FELinearSolidDomain.h"
#include "FERigidMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEElasticTrussDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FERemodelingElasticDomain.h"
#include "FEUDGHexDomain.h"
#include "FEUT4Domain.h"
#include "FELinearElastic.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEDiscreteSpringDomain.h"
#include "FERemodelingElasticMaterial.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FEElasticDomain2O.h"
#include "FESRIElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FESolidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FEDomain* pd = 0;
	FE_Element_Shape eshape = spec.eshape;
	FE_Element_Type etype = spec.etype;
	if (dynamic_cast<FERigidMaterial*>(pmat))
	{
		// rigid elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) pd = new FERigidSolidDomain(pfem);
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3)) pd = new FERigidShellDomain(pfem);
		else return 0;
	}
	else if (dynamic_cast<FERemodelingElasticMaterial*>(pmat)) pd = new FERemodelingElasticDomain(pfem);
	else if (dynamic_cast<FEMicroMaterial2O*          >(pmat)) pd = new FEElasticDomain2O        (pfem);
	else if (dynamic_cast<FESolidMaterial*>(pmat))
	{
		// structural elements
		if (eshape == ET_HEX8)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) pd = new FE3FieldElasticSolidDomain(pfem);
			else
			{
				if (etype == FE_HEX8G1) pd = new FEUDGHexDomain(pfem);
				else pd = new FEElasticSolidDomain(pfem);
			}
		}
		else if ((eshape == ET_TET10) || (eshape == ET_TET15))
		{
			if ((etype == FE_TET10G8RI4)||(etype == FE_TET10G4RI1))
			{
				pd = new FESRIElasticSolidDomain(pfem);
			}
			else
			{
				if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_tet)) pd = new FE3FieldElasticSolidDomain(pfem);
				pd = new FEElasticSolidDomain(pfem);
			}
		}
		else if ((eshape == ET_HEX20) || (eshape == ET_HEX27))
		{
//			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) pd = new FE3FieldElasticSolidDomain(pfem);
			pd = new FEElasticSolidDomain(pfem);
		}
		else if (eshape == ET_TET4)
		{
			if (spec.m_but4) pd = new FEUT4Domain(pfem);
			else pd = new FEElasticSolidDomain(pfem);
		}
		else if (eshape == ET_PENTA6) 
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) pd = new FE3FieldElasticSolidDomain(pfem);
			else pd = new FEElasticSolidDomain(pfem);
		}
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3)) pd = new FEElasticShellDomain(pfem);
		else if ((eshape == ET_TRUSS2)) pd = new FEElasticTrussDomain(pfem);
		else return 0;
	}
	else if (dynamic_cast<FEDiscreteMaterial*>(pmat))
	{
		// TODO: check the element type
		pd = new FEDiscreteSpringDomain(pfem);
	}

	if (pd) pd->SetMaterial(pmat);
	return pd;
}
