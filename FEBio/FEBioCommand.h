#pragma once
#include "Command.h"
#include "CommandManager.h"

class FEM;

//-----------------------------------------------------------------------------
//! Base class of FEBio commands

class FEBioCommand : public Command
{
public:
	FEBioCommand();
	virtual ~FEBioCommand(void);

	static void SetFEM(FEM* pfem);

protected:
	static FEM*	m_pfem;
};

//-----------------------------------------------------------------------------

class FERegisterCmd
{
public:
	FERegisterCmd(Command* pcmd, const char* szname, const char* szdesc) 
	{ 
		pcmd->SetName(szname);
		pcmd->SetDescription(szdesc);
		CommandManager* pCM = CommandManager::GetInstance();
		pCM->AddCommand(pcmd); 
	}
};

#define DECLARE_COMMAND(theCmd) public:	static FERegisterCmd m_##theCmd##_rc
#define REGISTER_COMMAND(theClass, theName, theDesc) FERegisterCmd theClass::m_##theClass##_rc(new theClass(), theName, theDesc)

//-----------------------------------------------------------------------------

class FEBioCmd_Help : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Help);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Quit : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Quit);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Cont : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Cont);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Conv : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Conv);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Debug : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Debug);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Dtmin : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Dtmin);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Fail : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Fail);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Plot : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Plot);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Print : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Print);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Restart : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Restart);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Version : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Version);
};

//-----------------------------------------------------------------------------

class FEBioCmd_Time : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Time);
};
