#include "stdafx.h"
#include "FEBioCodeSection.h"
#include "FECore/FECallBack.h"
#include "FECore/FECoreKernel.h"

void FEBioCodeSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if (tag == "callback")
		{
			const char* szname = tag.AttributeValue("name");
			FECallBack* pcb = fecore_new<FECallBack>(szname, GetFEModel());

			// TODO: The constructor of FECallBack already registered the callback class, so
			// we don't need to do anything else here. Of course, the question is who
			// is going to cleanup pcb?
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

