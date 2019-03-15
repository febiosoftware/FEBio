#include "stdafx.h"
#include "FEBioMeshAdaptorSection.h"
#include <FECore/FEMeshAdaptor.h>

void FEBioMeshAdaptorSection::Parse(XMLTag& tag)
{
	if (tag.isempty()) return;

	++tag;
	do
	{
		if (tag == "mesh_adaptor")
		{
			ParseMeshAdaptor(tag);
		}
		++tag;
	}
	while (!tag.isend());
}

void FEBioMeshAdaptorSection::ParseMeshAdaptor(XMLTag& tag)
{
	const char* sztype = tag.AttributeValue("type");

	FEModel* fem = GetFEModel();

	FEMeshAdaptor* meshAdaptor = fecore_new<FEMeshAdaptor>(sztype, fem);
	if (meshAdaptor == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	fem->AddMeshAdaptor(meshAdaptor);

	ReadParameterList(tag, meshAdaptor);
}
