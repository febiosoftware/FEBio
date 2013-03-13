// FEBioTM.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "FEBioTM.h"
#include "MainApp.h"
#include <flx_message.h>

#ifdef WIN32
int APIENTRY WinMain(HINSTANCE hInstance,
                     HINSTANCE hPrevInstance,
                     LPSTR     lpCmdLine,
                     int       nCmdShow)
#endif
#ifdef LINUX
int main(int argc, char* argv[])
#endif
#ifdef __APPLE__
int main(int argc, char* argv[])
#endif
{
	// the intial filename (if there is one)
	char szfilename[256] = {0};
	
#ifdef LINUX
	if (argc > 1) strcpy(szfilename, argv[1]);
#endif
#ifdef __APPLE__
	szfilename[0] = 0;
#endif
#ifdef WIN32
	strcpy(szfilename, lpCmdLine);
#endif

	// get the main application
	CMainApp* papp = FLXGetMainApp();

	// initialize the application
	if (papp->Init() == false)
	{
		flx_error("A fatal error has occured. The application will be terminated\n.");
		return -1;
	}

	// load a file if provided on the command line
	if (strlen(szfilename)>0) papp->Load(szfilename);

	// run the app
	return papp->Run();
}
