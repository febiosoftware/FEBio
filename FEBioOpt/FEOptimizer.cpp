#include "FEOptimizer.h"
#include "FELMOptimizeMethod.h"
#include "FECore/Logfile.h"

//=============================================================================
InvalidVariableName::InvalidVariableName(const char* sz)
{ 
	strcpy(szname, sz); 
}

//=============================================================================
// FEOptimizeInput
//=============================================================================

//-----------------------------------------------------------------------------
//! Read the data from the xml input file
//!
bool FEOptimizeInput::Input(const char* szfile, FEOptimizeData* pOpt)
{
	// try to open the file
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;

	// find the root element
	XMLTag tag;
	if (xml.FindTag("febio_optimize", tag) == false) return false;

	FEOptimizeData& opt = *pOpt;

	// parse the file
	try
	{
		++tag;
		bool bret = true;
		do
		{
			if		(tag == "Model"     ) bret = ParseModel     (tag, opt);
			else if (tag == "Options"   ) bret = ParseOptions   (tag, opt);
			else if (tag == "Function"  ) bret = ParseObjective (tag, opt);
			else if (tag == "Parameters") bret = ParseParameters(tag, opt);
			else if (tag == "LoadData"  ) bret = ParseLoadData  (tag, opt);
			else throw XMLReader::InvalidTag(tag);

			if (bret == false) return false;

			// go to the next tag
			++tag;
		}
		while (!tag.isend());
	}
	catch (InvalidVariableName e)
	{
		fprintf(stderr, "FATAL ERROR: the variable %s is not recognized.\n\n", e.szname);
		return false;
	}
	catch (NothingToOptimize)
	{
		fprintf(stderr, "FATAL ERROR: there is nothing to optimize.\n\n");
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		fprintf(stderr, "FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		fprintf(stderr, "FATAL ERROR: The element %s has an invalid value (%s).\n\n", e.tag.Name(), e.tag.szvalue());
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: an exception occured in the optimize routine.\n\n");
		return false;
	}
	xml.Close();

	return true;
}

//-----------------------------------------------------------------------------
//! Read the Model section of the input file
bool FEOptimizeInput::ParseModel(XMLTag& tag, FEOptimizeData& opt)
{
	FEModel& fem = opt.GetFEM();

	// get the input file
	char szfile[256];
	tag.value(szfile);

	// read input data
	if (fem.Input(szfile) == false) return false;

	// initialize data 
	if (fem.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Read the optimizer section of the input file
bool FEOptimizeInput::ParseOptions(XMLTag& tag, FEOptimizeData& opt)
{
	FELMOptimizeMethod* popt = 0;
	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) popt = new FELMOptimizeMethod;
	else
	{
		if (strcmp(szt, "marquardt") == 0) popt = new FELMOptimizeMethod;
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
	}

	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if      (tag == "obj_tol"     ) tag.value(popt->m_objtol);
			else if (tag == "f_diff_scale") tag.value(popt->m_fdiff );
			else if (tag == "max_iter"    ) tag.value(popt->m_nmax  );
//			else if (tag == "print_cov"   ) tag.value(popt->m_bcov  );
            else if (tag == "log_level"   )
            {
                char szval[256];
                tag.value(szval);
                if		(strcmp(szval, "LOG_DEFAULT"        ) == 0) {} // don't change the plot level
                else if (strcmp(szval, "LOG_NEVER"          ) == 0) popt->m_loglevel = Logfile::NEVER;
                else if (strcmp(szval, "LOG_FILE_ONLY"      ) == 0) popt->m_loglevel = Logfile::FILE_ONLY;
                else if (strcmp(szval, "LOG_SCREEN_ONLY"    ) == 0) popt->m_loglevel = Logfile::SCREEN_ONLY;
                else if (strcmp(szval, "LOG_FILE_AND_SCREEN") == 0) popt->m_loglevel = Logfile::FILE_AND_SCREEN;
                else throw XMLReader::InvalidValue(tag);
            }
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}

	opt.SetSolver(popt);

	return true;
}

//-----------------------------------------------------------------------------
//! Read the objectives section of the input file
bool FEOptimizeInput::ParseObjective(XMLTag &tag, FEOptimizeData& opt)
{
	FEModel& fem = opt.GetFEM();

	OPT_OBJECTIVE obj;

	++tag;
	do
	{
		if (tag == "fnc")
		{
			tag.value(obj.m_szname);

			// get the loadcurve for this objective function
			tag.AttributeValue("lc", obj.m_nlc);
			obj.m_nlc--;

			// find the variable
			obj.m_pd = fem.FindParameter(obj.m_szname);
			if (obj.m_pd == 0) throw InvalidVariableName(obj.m_szname);

			opt.SetObjective(obj);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Read the variables section of the input file
bool FEOptimizeInput::ParseParameters(XMLTag& tag, FEOptimizeData& opt)
{
	FEModel& fem = opt.GetFEM();

	// read the parameters
	OPT_VARIABLE var;
	++tag;
	do
	{
		if (tag == "param")
		{
			// get the variable name
			const char* sz = tag.AttributeValue("name");
			if (sz == 0) throw InvalidVariableName("[Unknown]");

			strcpy(var.m_szname, sz);

			// find the variable
			double* pd = fem.FindParameter(sz);
			if (pd == 0) throw InvalidVariableName(sz);
			
			var.m_pd = pd;

			// set initial values and bounds
			double d[4] = {0, 0, 0, 1};
			tag.value(d, 4);
			var.m_val = d[0];
			var.m_min = d[1];
			var.m_max = d[2];
			var.m_sf  = d[3];

			// add the variable
			opt.AddVariable(var);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Read the load data section of the input file
bool FEOptimizeInput::ParseLoadData(XMLTag &tag, FEOptimizeData& opt)
{
	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			FELoadCurve::INTFUNC ntype = FELoadCurve::LINEAR;
			FELoadCurve::EXTMODE nextm = FELoadCurve::CONSTANT;

			// get the (optional) type
			const char* szt = tag.AttributeValue("type", true);
			if (szt)
			{
				if      (strcmp(szt, "step"  ) == 0) ntype = FELoadCurve::STEP;
				else if (strcmp(szt, "linear") == 0) ntype = FELoadCurve::LINEAR;
				else if (strcmp(szt, "smooth") == 0) ntype = FELoadCurve::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
			}

			// get the optional extend mode
			const char* szm = tag.AttributeValue("extend", true);
			if (szm)
			{
				if      (strcmp(szm, "constant"     ) == 0) nextm = FELoadCurve::CONSTANT;
				else if (strcmp(szm, "extrapolate"  ) == 0) nextm = FELoadCurve::EXTRAPOLATE;
				else if (strcmp(szm, "repeat"       ) == 0) nextm = FELoadCurve::REPEAT;
				else if (strcmp(szm, "repeat offset") == 0) nextm = FELoadCurve::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", szt);
			}

			// count how many points we have
			XMLTag t(tag); ++t;
			int nlp = 0;
			while (!t.isend()) { ++nlp; ++t; }

			// create the loadcurve
			FELoadCurve* plc = new FELoadCurve;
			plc->Create(nlp);
			plc->SetInterpolation(ntype);
			plc->SetExtendMode(nextm);
			opt.AddLoadCurve(plc);

			// read the points
			double d[2];
			++tag;
			for (int i=0; i<nlp; ++i)
			{
				tag.value(d, 2);
				plc->LoadPoint(i).time  = d[0];
				plc->LoadPoint(i).value = d[1];

				++tag;
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
// FEOptimizeData
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
FEOptimizeData::FEOptimizeData(FEModel& fem) : m_fem(fem)
{
	m_pSolver = 0;
}

//-----------------------------------------------------------------------------
FEOptimizeData::~FEOptimizeData(void)
{
	delete m_pSolver;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Init()
{
	if (m_pSolver == 0) m_pSolver = new FELMOptimizeMethod;

	return true;
}

//-----------------------------------------------------------------------------
bool FEOptimizeData::Solve()
{
	return m_pSolver->Solve(this);
}

//-----------------------------------------------------------------------------
//! Read the data from the input file
//!
bool FEOptimizeData::Input(const char *szfile)
{
	FEOptimizeInput in;
	if (in.Input(szfile, this) == false) return false;
	return true;
}
