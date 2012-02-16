// FEBioTM.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "FEBioTM.h"
#include "MainApp.h"

int APIENTRY _tWinMain(HINSTANCE hInstance,
                     HINSTANCE hPrevInstance,
                     LPTSTR    lpCmdLine,
                     int       nCmdShow)
{
	// get the main application
	CMainApp* papp = FLX_GetMainApp();

	// initialize application
	papp->Init();

	// run the app
	return papp->Run();
}
