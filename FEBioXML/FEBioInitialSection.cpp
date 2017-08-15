#include "stdafx.h"
#include "FEBioInitialSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"
#include <FECore/FEInitialCondition.h>
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
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

			// add it to the model
			GetBuilder()->AddInitialCondition(pic);

			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = ReadNodeID(tag);
					vec3d v;
					value(tag, v);

					pic->Add(nid, v);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "ic")
		{
			const char* sztype = tag.AttributeValue("type");
			FEInitialCondition* pic = dynamic_cast<FEInitialCondition*>(fecore_new<FEInitialCondition>(FEIC_ID, sztype, &fem));

			if (tag.isleaf() == false)
			{
				FEParameterList& pl = pic->GetParameterList();
				++tag;
				do
				{
					if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}
			
			// add it to the model
			GetBuilder()->AddInitialCondition(pic);
		}
		else
		{
			// Get the degree of freedom
			int ndof = -1;
			if      (tag == "temperature"   ) ndof = dofs.GetDOF("T");
			else if (tag == "fluid_pressure") ndof = dofs.GetDOF("p");
            else if (tag == "fluid_shell_pressure") ndof = dofs.GetDOF("q");
            else if (tag == "dilatation"    ) ndof = dofs.GetDOF("e");
			else if (tag == "concentration" )
			{
				// TODO: Add a check to make sure that a solute with this ID exists
				int isol = 0;
				const char* sz = tag.AttributeValue("sol", true);
				if (sz) isol = atoi(sz) - 1;

				ndof = dofs.GetDOF("concentration", isol);
			}
            else if (tag == "shell_concentration" )
            {
                // TODO: Add a check to make sure that a solute with this ID exists
                int isol = 0;
                const char* sz = tag.AttributeValue("sol", true);
                if (sz) isol = atoi(sz) - 1;
                
                ndof = dofs.GetDOF("shell concentration", isol);
            }
			else throw XMLReader::InvalidTag(tag);
			if (ndof == -1) throw XMLReader::InvalidTag(tag);

			// allocate initial condition
			FEInitialBC* pic = dynamic_cast<FEInitialBC*>(fecore_new<FEInitialCondition>(FEIC_ID, "init_bc", &fem));
			pic->SetDOF(ndof);

			// add it to the model
			GetBuilder()->AddInitialCondition(pic);

			// read the node list and values
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = ReadNodeID(tag);
					double p;
					value(tag, p);
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

//-----------------------------------------------------------------------------
void FEBioInitialSection25::Parse(XMLTag& tag)
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
		if (tag == "init")
		{
			// get the BC
			const char* sz = tag.AttributeValue("bc");
			int ndof = dofs.GetDOF(sz);
			if (ndof == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			// get the node set
			const char* szset = tag.AttributeValue("node_set");
			FENodeSet* pns = mesh.FindNodeSet(szset);
			if (pns == 0) throw XMLReader::InvalidTag(tag);

			// allocate initial condition
			FEInitialBC* pic = dynamic_cast<FEInitialBC*>(fecore_new<FEInitialCondition>(FEIC_ID, "init_bc", &fem));
			pic->SetDOF(ndof);
			pic->SetNodes(*pns);

			// add it to the model
			GetBuilder()->AddInitialCondition(pic);

			// read parameters
			ReadParameterList(tag, pic);
		}
		else if (tag == "ic")
		{
			const char* sztype = tag.AttributeValue("type");
			FEInitialCondition* pic = dynamic_cast<FEInitialCondition*>(fecore_new<FEInitialCondition>(FEIC_ID, sztype, &fem));

			ReadParameterList(tag, pic);
			
			// add it to the model
			GetBuilder()->AddInitialCondition(pic);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
