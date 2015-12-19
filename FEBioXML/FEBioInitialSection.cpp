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
	DOFS& dofs = fem.GetDOFS();

	// make sure we've read the nodes section
	if (mesh.Nodes() == 0) throw XMLReader::InvalidTag(tag);

	// read nodal data
	++tag;
	do
	{
		if (tag == "velocity")
		{
			const int dof_VX = dofs.GetDOF("vx");
			const int dof_VY = dofs.GetDOF("vy");
			const int dof_VZ = dofs.GetDOF("vz");
			FEInitialBCVec3D* pic = dynamic_cast<FEInitialBCVec3D*>(fecore_new<FEInitialCondition>(FEIC_ID, "init_bc_vec3d", &fem));
			pic->SetDOF(dof_VX, dof_VY, dof_VZ);
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
		else
		{
			// Get the degree of freedom
			int ndof = -1;
			if      (tag == "temperature"   ) ndof = dofs.GetDOF("t");
			else if (tag == "fluid_pressure") ndof = dofs.GetDOF("p");
			else if (tag == "concentration" )
			{
				// TODO: Add a check to make sure that a solute with this ID exists
				int isol = 0;
				const char* sz = tag.AttributeValue("sol", true);
				if (sz) isol = atoi(sz) - 1;

				ndof = dofs.GetDOF("concentration", isol);
			}
			else if (tag == "init")
			{
				// TODO: In the future I want to make this the preferred way of defining initial conditions
				// get the bc attribute
				const char* szbc = tag.AttributeValue("bc");
				ndof = dofs.GetDOF(szbc);
			}
			else throw XMLReader::InvalidTag(tag);
			if (ndof == -1) throw XMLReader::InvalidTag(tag);

			// allocate initial condition
			FEInitialBC* pic = dynamic_cast<FEInitialBC*>(fecore_new<FEInitialCondition>(FEIC_ID, "init_bc", &fem));
			pic->SetDOF(ndof);

			// add this boundary condition to the current step
			fem.AddInitialCondition(pic);
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddModelComponent(pic);
				pic->Deactivate();
			}

			// read the node list and values
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
		++tag;
	}
	while (!tag.isend());
}
