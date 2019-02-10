#pragma once
#include "cmdoptions.h"

class FEBioModel;

class FEBioApp
{
public:
	FEBioApp();

	int Init(int nargs, char* argv[]);

	int Run();

	void Finish();

	void ProcessCommands();

	CMDOPTIONS& CommandOptions();

	bool ParseCmdLine(int argc, char* argv[]);

	// run an febio model
	int RunModel();

public:
	// get the current model
	FEBioModel* GetCurrentModel();

protected:
	// show FEBio prompt
	int prompt();

	// set the currently active model
	void SetCurrentModel(FEBioModel* fem);

public:
	static FEBioApp* GetInstance();

private:
	CMDOPTIONS	m_ops;				// command line options

	FEBioModel*		m_fem;			// current model (or null if not model is running)

	static FEBioApp*	m_This;
};
