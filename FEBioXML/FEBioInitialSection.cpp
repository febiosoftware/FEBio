#include "stdafx.h"
#include "FEBioInitialSection.h"

//-----------------------------------------------------------------------------
//! Read the Initial from the FEBio input file
//!
void FEBioInitialSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// make sure we've read the nodes section
	if (mesh.Nodes() == 0) throw XMLReader::InvalidTag(tag);

	for (int i=0; i<mesh.Nodes(); ++i) mesh.Node(i).m_v0 = vec3d(0,0,0);

	// read nodal data
	++tag;
	do
	{
		if (tag == "velocity")
		{
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					vec3d v;
					tag.value(v);
					mesh.Node(nid).m_v0 += v;
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "fluid_pressure")
		{
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					double p;
					tag.value(p);
					mesh.Node(nid).m_p0 += p;
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "concentration")
		{
			int isol = 0;
			const char* sz = tag.AttributeValue("sol", true);
			if (sz) isol = atoi(sz) - 1;
			if ((isol < 0) || (isol >= MAX_CDOFS))
				throw XMLReader::InvalidAttributeValue(tag, "sol", sz);
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					double c;
					tag.value(c);
					mesh.Node(nid).m_c0[isol] += c;
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
