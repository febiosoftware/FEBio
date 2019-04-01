#include "stdafx.h"
#include "FEBioApp.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

//-----------------------------------------------------------------------------
// The starting point of the application
//
int main(int argc, char* argv[])
{
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#endif

	// create the FEBio app
	FEBioApp febio;

	// initialize the app
	if (febio.Init(argc, argv) == false) return 1;

	// start the FEBio app
	int nret = febio.Run();

	// Don't forget to cleanup
	febio.Finish();

#ifdef USE_MPI
	MPI_Finalize();
#endif

	return nret;
}
