// fem.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FEM_H_07012006_
#define _FEM_H_07012006_

#include "FEBioLib/FEBioModel.h"
#include "FECore/DumpFile.h"
#include "FECore/FEContactInterface.h"
#include "FEBioLib/DataStore.h"
#include "FEBioPlot/PlotFile.h"
#include "FEBioLib/FEElasticMixture.h"
#include "FEBioLib/FEUncoupledElasticMixture.h"
#include "FEBioLib/FESolute.h"

#include <list>
using namespace std;

//-----------------------------------------------------------------------------
//! The Finite Element Model class. 

//! This class stores solver parameters, geometry data, material data, and 
//! other data that is needed to solve the FE problem.
//! FEBio is designed to solve finite element problems. All the finite element
//! data is collected here in this class. This class also provides
//! routines to initalize, input, output and update the FE data. Although this
//! class provides the main solve routine it does not really solve anything.
//! The actual solving is done by one of the classes derived from the FESolver class.

class FEM : public FEBioModel
{
public:
	//! constructor - sets default variables
	FEM();

	//! destructor
	virtual ~FEM();

	//! read the configuration file
	bool Configure(const char* szfile);

	//! Restart from restart point
	bool Restart(const char* szfile);

	//! Serialize the current state to/from restart file
	bool Serialize(DumpFile& ar);

	//! return a pointer to the named variable
	double* FindParameter(const char* szname);

	//! return a pointer to the parameter variable
	double* ReturnParameter(FEParam* pp, const int index);
	
	//! return a pointer to the named variable in a solid mixture
	double* FindSolidMixtureParameter(const char* szvar, const int index, FEElasticMixture* pme);
	double* FindUncoupledSolidMixtureParameter(const char* szvar, const int index, FEUncoupledElasticMixture* pme);
	
	//! find a boundary condition from the ID
	FEBoundaryCondition* FindBC(int nid);

	//! Evaluate parameter list
	void EvaluateParameterList(FEParameterList& pl);
	void EvaluateMaterialParameters(FEMaterial* pm);

	//! check for user-interruption
	virtual void CheckInterruption();

public:
	virtual void PushState();
	virtual void PopState ();

protected:
	void ShallowCopy(FEM& fem);

protected:
	void SerializeMaterials   (DumpFile& ar);
	void SerializeAnalysisData(DumpFile& ar);
	void SerializeGeometry    (DumpFile& ar);
	void SerializeMesh        (DumpFile& ar);
	void SerializeContactData (DumpFile& ar);
	void SerializeBoundaryData(DumpFile& ar);
	void SerializeIOData      (DumpFile& ar);
	void SerializeLoadData    (DumpFile& ar);
	void SerializeConstants   (DumpFile& ar);
	void SerializeDataStore   (DumpFile& ar);

public:
	bool	m_bInterruptable;	//!< true if this model can be interrupted with ctrl+c
};

#endif // _FEM_H_07012006_
