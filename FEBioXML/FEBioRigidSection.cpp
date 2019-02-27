#include "stdafx.h"
#include "FEBioRigidSection.h"
#include <FECore/FEModel.h>
#include <FEBioMech/FERigidSystem.h>
#include <FECore/FECoreKernel.h>
#include <FEBioMech/FERigidSurface.h>
#include <FEBioMech/FEMechModel.h>

void FEBioRigidSection::Parse(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& RS = *fem.GetRigidSystem();

	++tag;
	do
	{
		if (tag == "rigid")
		{
			const char* sztype = tag.AttributeValue("type");
			const char* szname = tag.AttributeValue("name");
			FERigidSurface* rs = fecore_new<FERigidSurface>(sztype, &fem);
			if (rs == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			rs->SetName(szname);
			RS.AddRigidSurface(rs);

			FEParameterList& pl = rs->GetParameterList();
			ReadParameterList(tag, pl);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
