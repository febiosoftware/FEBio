///////////////////////////////////////////////////////////////////////////////
//         ________    _________   _________     __     _________            //
//        |        |\ |        |\ |        |\   |  |\  /         \\          //
//        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||         //
//        |   |\___\| |   |\___\| |   |\_| ||    \_\| |   //  \   ||         //
//        |   ||      |   ||      |   || | ||    __   |  ||    |  ||         //
//        |   ||__    |   ||__    |   ||_| ||   |  |\ |  ||    |  ||         //
//        |       |\  |       |\  |         \\  |  || |  ||    |  ||         //
//        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||         //
//        |   |\__\|  |   |\__\|  |   |\__|  || |  || |  ||    |  ||         //
//        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||         //
//        |   ||      |   ||___   |   ||__|  || |  || |   \\__/   ||         //
//        |   ||      |        |\ |          || |  || |           ||         //
//        |___||      |________|| |__________|| |__||  \_________//          //
//                                                                           //
//      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
// FEBio is a finite element solver that is specifically designed for three
// dimensional biomechanical applications. It solves the nonlinear finite
// element equations using a quasi-Newton method called the BFGS-method. It
// also offers several biologically relevant constitutive models.
//
// This software is developed at the Musculoskeletal Research Laboratories
// at the University of Utah. FEBio is a registered trademark. All rights reserved. 
// Copyright (c) 2006 - 2019
//
// The subversion (svn) revision number of this code can be found in the file
// FEBio/svnrev.h
//
// Main developers:
//  - Steve Maas
//  - Gerard Ateshian
//  - Jeff Weiss
//  - Dave Rawlins
//
// Contributors:
//  - Alexander Veress
//
///////////////////////////////////////////////////////////////////////////////

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
