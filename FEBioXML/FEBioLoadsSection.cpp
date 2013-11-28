#include "stdafx.h"
#include "FEBioLoadsSection.h"
#include "FEBioMech/FEPointBodyForce.h"
#include "FEBioMech/FEPressureLoad.h"
#include "FEBioMech/FETractionLoad.h"
#include "FEBioMix/FEPoroTraction.h"
#include "FEBioMix/FEFluidFlux.h"
#include "FEBioMix/FESoluteFlux.h"
#include "FEBioHeat/FEHeatSource.h"
#include "FEBioHeat/FEHeatFlux.h"
#include "FEBioHeat/FEConvectiveHeatFlux.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//!  Parses the loads section from the xml file (version 1.2 or up)
//!
void FEBioLoadsSection::Parse(XMLTag& tag)
{
	assert(m_pim->Version() >= 0x0102);
	
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "force"              ) ParseBCForce             (tag);
		else if (tag == "pressure"           ) ParseBCPressure          (tag);
		else if (tag == "traction"           ) ParseBCTraction          (tag);
		else if (tag == "normal_traction"    ) ParseBCPoroNormalTraction(tag);
		else if (tag == "fluidflux"          ) ParseBCFluidFlux         (tag);
		else if (tag == "soluteflux"         ) ParseBCSoluteFlux        (tag);
		else if (tag == "heatflux"           ) ParseBCHeatFlux          (tag);
		else if (tag == "convective_heatflux") ParseBCConvectiveHeatFlux(tag);
		else if (tag == "body_force"         ) ParseBodyForce           (tag);
		else if (tag == "heat_source"        ) ParseHeatSource          (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// NOTE: note that this section used to be in the Globals section (version 1.1)
void FEBioLoadsSection::ParseBodyForce(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) szt = "const";

	if (strcmp(szt, "point") == 0)
	{
		FEPointBodyForce* pf = new FEPointBodyForce(&fem);
		FEParameterList& pl = pf->GetParameterList();
		++tag;
		do
		{
			if (tag == "a")
			{
				const char* szlc = tag.AttributeValue("lc");
//						pf->lc[0] = pf->lc[1] = pf->lc[2] = atoi(szlc);
				tag.value(pf->m_a);
			}
			else if (tag == "node")
			{
				tag.value(pf->m_inode); 
				pf->m_inode -= 1;
			}
			else if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());

		fem.AddBodyLoad(pf);
	}
	else
	{
		// see if the kernel knows this force
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEBodyForce* pf = febio.Create<FEBodyForce>(szt, &fem);
		if (pf)
		{
			if (!tag.isleaf())
			{
				FEParameterList& pl = pf->GetParameterList();
				++tag;
				do
				{
					if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}

			fem.AddBodyLoad(pf);
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseHeatSource(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEHeatSource* phs = new  FEHeatSource(&fem);
	FEParameterList& PL = phs->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, PL) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
	fem.AddBodyLoad(phs);
}


//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCForce(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	int nversion = m_pim->Version();
	if (nversion >= 0x0200)
	{
		// count how many nodal forces there are
		int ncnf = tag.children();

		// get the bc
		int bc = -1;
		const char* sz = tag.AttributeValue("bc");

		if      (strcmp(sz, "x") == 0) bc = DOF_X;
		else if (strcmp(sz, "y") == 0) bc = DOF_Y;
		else if (strcmp(sz, "z") == 0) bc = DOF_Z;
		else if (strcmp(sz, "p") == 0) bc = DOF_P;
		else if (strcmp(sz, "t") == 0) bc = DOF_T;
		else if (strcmp(sz, "c") == 0) bc = DOF_C;
		else if (strcmp(sz, "c1") == 0) bc = DOF_C;
		else if (strcmp(sz, "c2") == 0) bc = DOF_C+1;
		else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

		// read the prescribed data
		++tag;
		for (int i=0; i<ncnf; ++i)
		{
			// get the nodal ID
			int n = atoi(tag.AttributeValue("id"))-1;

			// get the load curve
			sz = tag.AttributeValue("lc");
			int lc = atoi(sz)-1;

			// create new nodal force
			FENodalForce* pfc = new FENodalForce;
			pfc->node = n;
			pfc->bc = bc;
			pfc->lc = lc;
			tag.value(pfc->s);
			fem.AddNodalLoad(pfc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(pfc);
				pfc->Deactivate();
			}

			++tag;
		}
	}
	else
	{
		// count how many nodal forces there are
		int ncnf = tag.children();

		// read the prescribed data
		++tag;
		for (int i=0; i<ncnf; ++i)
		{
			int n = atoi(tag.AttributeValue("id"))-1, bc;
			const char* sz = tag.AttributeValue("bc");

			if      (strcmp(sz, "x") == 0) bc = 0;
			else if (strcmp(sz, "y") == 0) bc = 1;
			else if (strcmp(sz, "z") == 0) bc = 2;
			else if (strcmp(sz, "p") == 0) bc = 6;
			else if (strcmp(sz, "t") == 0) bc = 10;
			else if (strcmp(sz, "c") == 0) bc = 11;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

			sz = tag.AttributeValue("lc");
			int lc = atoi(sz) - 1;

			FENodalForce* pfc = new FENodalForce;
			pfc->node = n;
			pfc->bc = bc;
			pfc->lc = lc;
			tag.value(pfc->s);
			fem.AddNodalLoad(pfc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddBoundaryCondition(pfc);
				pfc->Deactivate();
			}

			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCPressure(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}

	// count how many pressure cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate pressure data
	FEPressureLoad* ps = new FEPressureLoad(psurf, blinear);
	ps->create(npr);
	fem.AddSurfaceLoad(ps);

	// read the pressure data
	++tag;
	int nf[FEElement::MAX_NODES ], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEPressureLoad::LOAD& pc = ps->PressureLoad(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		pc.lc = atoi(sz)-1;

		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ps);
		ps->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCTraction(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;

	// count how many traction cards there are
	int ntc = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(ntc);
	fem.GetMesh().AddSurface(psurf);

	// allocate traction data
	FETractionLoad* pt = new FETractionLoad(psurf);
	fem.AddSurfaceLoad(pt);
	pt->create(ntc);

	// read the traction data
	++tag;
	int nf[FEElement::MAX_NODES], N;
	vec3d s;
	for (int i=0; i<ntc; ++i)
	{
		FETractionLoad::LOAD& tc = pt->TractionLoad(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		tc.lc = atoi(sz)-1;

		s.x  = atof(tag.AttributeValue("tx"));
		s.y  = atof(tag.AttributeValue("ty"));
		s.z  = atof(tag.AttributeValue("tz"));

		tc.s[0] = tc.s[1] = tc.s[2] = tc.s[3] = s;
		tc.s[4] = tc.s[5] = tc.s[6] = tc.s[7] = s;

		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pt);
		pt->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCPoroNormalTraction(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	
	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	bool beffective = false;
	sz = tag.AttributeValue("traction", true);
	if (sz)
	{
		if (strcmp(sz, "effective") == 0) beffective = true;
		else if ((strcmp(sz, "total") == 0) || (strcmp(sz, "mixture") == 0)) beffective = false;
		else throw XMLReader::InvalidAttributeValue(tag, "traction", sz);
	}
	
	// count how many normal traction cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate normal traction data
	FEPoroNormalTraction* ps = new FEPoroNormalTraction(psurf, blinear, beffective);
	ps->create(npr);
	fem.AddSurfaceLoad(ps);
	
	// read the normal traction data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEPoroNormalTraction::LOAD& pc = ps->NormalTraction(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		pc.lc = atoi(sz)-1;
		
		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;
		
		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);
		
		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;
		
		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ps);
		ps->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCFluidFlux(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	bool bmixture = false;
	sz = tag.AttributeValue("flux", true);
	if (sz)
	{
		if (strcmp(sz, "mixture") == 0) bmixture = true;
		else if (strcmp(sz, "fluid") == 0) bmixture = false;
		else throw XMLReader::InvalidAttributeValue(tag, "flux", sz);
	}
	
	// count how many fluid flux cards there are
	int nfr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(nfr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate fluid flux data
	FEFluidFlux* pfs = new FEFluidFlux(psurf, blinear, bmixture);
	pfs->create(nfr);
	fem.AddSurfaceLoad(pfs);
	
	// read the fluid flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<nfr; ++i)
	{
		FEFluidFlux::LOAD& fc = pfs->FluidFlux(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		fc.lc = atoi(sz) - 1;
		
		s  = atof(tag.AttributeValue("scale"));
		fc.s[0] = fc.s[1] = fc.s[2] = fc.s[3] = s;
		fc.s[4] = fc.s[5] = fc.s[6] = fc.s[7] = s;
		
		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);
		
		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;
		
		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pfs);
		pfs->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCSoluteFlux(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	
	const char* sz;
	bool blinear = false;
	sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if (strcmp(sz, "linear") == 0) blinear = true;
		else if (strcmp(sz, "nonlinear") == 0) blinear = false;
		else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	
	// count how many fluid flux cards there are
	int nfr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(nfr);
	fem.GetMesh().AddSurface(psurf);
	
	// allocate fluid flux data
	FESoluteFlux* pfs = new FESoluteFlux(psurf, blinear);
	pfs->create(nfr);
	fem.AddSurfaceLoad(pfs);
	
	// read the fluid flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<nfr; ++i)
	{
		FESoluteFlux::LOAD& fc = pfs->SoluteFlux(i);
		FESurfaceElement& el = psurf->Element(i);
		
		sz = tag.AttributeValue("lc");
		fc.lc = atoi(sz) - 1;
		
		s  = atof(tag.AttributeValue("scale"));
		fc.s[0] = fc.s[1] = fc.s[2] = fc.s[3] = s;
		fc.s[4] = fc.s[5] = fc.s[6] = fc.s[7] = s;
		
		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);
		
		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;
		
		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(pfs);
		pfs->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCHeatFlux(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many heatflux cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate flux data
	FEHeatFlux* ph = new FEHeatFlux(psurf);
	ph->create(npr);
	fem.AddSurfaceLoad(ph);

	const char* sz;

	// read the flux data
	++tag;
	int nf[4], N;
	double s;
	for (int i=0; i<npr; ++i)
	{
		FEHeatFlux::LOAD& pc = ph->HeatFlux(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		if (sz) pc.lc = atoi(sz) - 1;

		s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ph);
		ph->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection::ParseBCConvectiveHeatFlux(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many heatflux cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = new FESurface(&fem.GetMesh());
	psurf->create(npr);
	fem.GetMesh().AddSurface(psurf);

	// allocate flux data
	FEConvectiveHeatFlux* ph = new FEConvectiveHeatFlux(psurf);
	ph->create(npr);
	fem.AddSurfaceLoad(ph);

	const char* sz;

	// read the flux data
	++tag;
	int nf[8], N;
	for (int i=0; i<npr; ++i)
	{
		FEConvectiveHeatFlux::LOAD& pc = ph->HeatFlux(i);
		FESurfaceElement& el = psurf->Element(i);

		sz = tag.AttributeValue("lc");
		if (sz) pc.lc = atoi(sz) - 1;

		double s  = atof(tag.AttributeValue("scale"));
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;

		// heat transfer coefficient
		double hc = atof(tag.AttributeValue("hc"));
		pc.hc = hc;

		// read the element
		if      (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
//		else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
//		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

		++tag;
	}

	// add this boundary condition to the current step
	if (m_pim->m_nsteps > 0)
	{
		GetStep()->AddBoundaryCondition(ph);
		ph->Deactivate();
	}
}
