/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "Command.h"
#include "CommandManager.h"

//-----------------------------------------------------------------------------
class FEBioModel;

//-----------------------------------------------------------------------------
//! Base class of FEBio commands

class FEBioCommand : public Command
{
public:
	FEBioCommand();
	virtual ~FEBioCommand(void);

protected:
	FEBioModel* GetFEM();
};

//-----------------------------------------------------------------------------
// helper class for registering FEBio commands
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
class FEBioCmd_Run : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Run);
};

//-----------------------------------------------------------------------------
class FEBioCmd_Restart : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Restart);
};

//-----------------------------------------------------------------------------
class FEBioCmd_LoadPlugin : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_LoadPlugin);
};

//-----------------------------------------------------------------------------
class FEBioCmd_UnLoadPlugin : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_UnLoadPlugin);
};

//-----------------------------------------------------------------------------
class FEBioCmd_Config : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Config);
};

//-----------------------------------------------------------------------------
class FEBioCmd_Plugins : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Plugins);
};

//-----------------------------------------------------------------------------
class FEBioCmd_Help : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Help);
};

//-----------------------------------------------------------------------------
class FEBioCmd_Events : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_Events);
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

//-----------------------------------------------------------------------------
class FEBioCmd_svg : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_svg);
};

//-----------------------------------------------------------------------------
class FEBioCmd_out : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_out);
};

//-----------------------------------------------------------------------------
class FEBioCmd_where : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_where);
};

//-----------------------------------------------------------------------------
class FEBioCmd_break : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_break);
};

//-----------------------------------------------------------------------------
class FEBioCmd_breaks : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_breaks);
};

//-----------------------------------------------------------------------------
class FEBioCmd_clear_breaks : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_clear_breaks);
};

//-----------------------------------------------------------------------------
class FEBioCmd_list : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_list);
};

//-----------------------------------------------------------------------------
class FEBioCmd_hist: public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_hist);
};

//-----------------------------------------------------------------------------
class FEBioCmd_set : public FEBioCommand
{
public:
	int run(int nargs, char** argv);
	DECLARE_COMMAND(FEBioCmd_set);
};
