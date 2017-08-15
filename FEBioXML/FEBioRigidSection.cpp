#include "stdafx.h"
#include "FEBioRigidSection.h"
#include <FECore/FEModel.h>
#include <FECore/FERigidSystem.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FERigidSurface.h>

void FEBioRigidSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FERigidSystem& RS = *fem.GetRigidSystem();

	++tag;
	do
	{
		if (tag == "rigid")
		{
			const char* sztype = tag.AttributeValue("type");
			const char* szname = tag.AttributeValue("name");
			FERigidSurface* rs = fecore_new<FERigidSurface>(FERIGIDOBJECT_ID, sztype, &fem);
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
