#include "stdafx.h"
#include "FEBioInitialSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
//! Read the Initial from the FEBio input file
//!
void FEBioInitialSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

    // get number of DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetCDOFS();
    
	// make sure we've read the nodes section
	if (mesh.Nodes() == 0) throw XMLReader::InvalidTag(tag);

	// read nodal data
	++tag;
	do
	{
		if (tag == "velocity")
		{
			FEInitialVelocity* pic = dynamic_cast<FEInitialVelocity*>(fecore_new<FEInitialCondition>(FEIC_ID, "velocity", &fem));
			fem.AddInitialCondition(pic);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					vec3d v;
					m_pim->value(tag, v);

					pic->Add(nid, v);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "fluid_pressure")
		{
			FEInitialPressure* pic = dynamic_cast<FEInitialPressure*>(fecore_new<FEInitialCondition>(FEIC_ID, "fluid_pressure", &fem));
			fem.AddInitialCondition(pic);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					double p;
					m_pim->value(tag, p);
					pic->Add(nid, p);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "concentration")
		{
			FEInitialConcentration* pic = dynamic_cast<FEInitialConcentration*>(fecore_new<FEInitialCondition>(FEIC_ID, "concentration", &fem));
			fem.AddInitialCondition(pic);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			int isol = 0;
			const char* sz = tag.AttributeValue("sol", true);
			if (sz) isol = atoi(sz) - 1;
			if ((isol < 0) || (isol >= MAX_CDOFS))
				throw XMLReader::InvalidAttributeValue(tag, "sol", sz);

			pic->SetSoluteID(isol);

			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					double c;
					m_pim->value(tag, c);
					pic->Add(nid, c);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "temperature")
		{
			FEInitialTemperature* pic = dynamic_cast<FEInitialTemperature*>(fecore_new<FEInitialCondition>(FEIC_ID, "temperature", &fem));
			fem.AddInitialCondition(pic);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					double T;
					m_pim->value(tag, T);
					pic->Add(nid, T);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "ic")
		{
			// get the type attribute
			const char* sztype = tag.AttributeValue("type");

			// allocate a new initial condition of this type
			FEInitialCondition* pic = fecore_new<FEInitialCondition>(FEIC_ID, sztype, &fem);
			if (pic == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// add the initial condition
			fem.AddInitialCondition(pic);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			if (!tag.isempty())
			{
				// TODO: read the parameter list
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
