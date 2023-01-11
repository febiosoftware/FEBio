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

#include "sphericalHarmonics.h"
#include "spherePoints.h"
#include <algorithm>
#include <unordered_map>
#include <vector>

#ifdef HAS_MMG
#include <mmg/mmgs/libmmgs.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

using std::vector;
using std::unordered_map;

enum NUMTYPE { REALTYPE, IMAGTYPE, COMPLEXTYPE };

double fact(int val)
{
    double ans = 1;
    
    for(int i = 1; i <= val; i++)
    {
        ans *= i;
    } 
    
    return ans;
}

void getSphereCoords(int numPts, const double* xCoords, const double* yCoords, const double* zCoords, double* theta, double* phi)
{
    // get spherical coordinates
    for(int index = 0; index < numPts; index++)
    {
        double val = zCoords[index];

        if (val<=-1)
        {
            theta[index] = M_PI;
        } 
        else if (val>=1)
        {
             theta[index] = 0;
        } 
        else 
        {
            theta[index] = acos(zCoords[index]);
        }
    }

    for(int index = 0; index < numPts; index++)
    {
        double x = xCoords[index];
        double y = yCoords[index];

        if(x > 0)
        {
            phi[index] = atan(y/x);
        }
        else if(x < 0 && y >= 0)
        {
            phi[index] = atan(y/x) + M_PI;
        }
        else if(x < 0 && y < 0)
        {
            phi[index] = atan(y/x) - M_PI;
        }
        else if(x == 0 && y > 0)
        {
            phi[index] = M_PI/2;
        }
        else if(x == 0 && y < 0)
        {
            phi[index] = -M_PI/2;
        }
    }

}

std::unique_ptr<matrix> compSH(int order, int numPts, double* theta, double* phi)
{
    int numCols = (order+1)*(order+2)/2;
    std::unique_ptr<matrix> out = std::make_unique<matrix>(numPts, numCols);
    out->fill(0,0,numPts, numCols, 0.0);

    for(int k = 0; k <= order; k+=2)
    {
        for(int m = -k; m <= k; m++)
        {
            int j = (k*k + k + 2)/2 + m - 1;

            int numType = COMPLEXTYPE;
            double factor = 1;
            if(m < 0)
            {
                numType = REALTYPE;
                factor = sqrt(2);
            }
            else if(m > 0)
            {
                numType = IMAGTYPE;
                factor = sqrt(2);
            }

            for(int index = 0; index < numPts; index++)
            {
                (*out)[index][j] = factor*harmonicY(k, m, theta[index], phi[index], numType);
            }
        }
    }

    return std::move(out);
}

double harmonicY(int degree, int order, double theta, double phi, int numType)
{
    // Will be true if order is both positive and odd
    // In that case, we need to negate the answer to match the output
    // in Adam's MATLAB code.
    int negate = 1;
    if(order % 2 == 1) negate = -1;
    
    order = abs(order);

    double a = (2*degree+1)/(4*M_PI);
    double b = fact(degree-order)/fact(degree+order);

    double normalization = sqrt(a*b);

    double e;
    switch (numType)
    {
    case REALTYPE:
        e = cos(order*phi);
        break;
    case IMAGTYPE:
        e = sin(order*phi);
        break;
    case COMPLEXTYPE:
        e = 1;
        break;
    
    default:
        break;
    }

    return normalization*std::assoc_legendre(degree, order, cos(theta))*pow(-1, degree)*e*negate;
}

void reconstructODF(std::vector<double>& sphHarm, std::vector<double>& ODF, int numPts, double* theta, double* phi)
{
    int order = (sqrt(8*sphHarm.size() + 1) - 3)/2;

    ODF.resize(numPts);  
    auto T = compSH(order, numPts, theta, phi);
    (*T).mult(sphHarm, ODF);

    // Normalize ODF
    double sum = 0;
    for(int index = 0; index < ODF.size(); index++)
    {
        if((ODF)[index] < 0)
        {
            ODF[index] = 0;
        }
        
        sum += ODF[index];
    }

    for(int index = 0; index < NPTS; index++)
    {
        ODF[index] /= sum;
    }
}

void altGradient(int order, std::vector<double>& ODF, std::vector<double>& gradient)
{
    gradient.resize(NPTS);
    for(int index = 0; index < NPTS; index++)
    {
        gradient[index] = 0;
    }

    std::vector<int> count(NPTS, 0);

    for(int index = 0; index < NCON; index++)
    {
        int n0 = CONN1[index]-1;
        int n1 = CONN2[index]-1;
        int n2 = CONN3[index]-1;

        double val0 = ODF[n0];
        double val1 = ODF[n1];
        double val2 = ODF[n2];

        double diff0 = abs(val0 - val1);
        double diff1 = abs(val0 - val2);
        double diff2 = abs(val1 - val2);

        gradient[n0] += diff0 + diff1;
        gradient[n1] += diff0 + diff2;
        gradient[n2] += diff1 + diff2;

        count[n0]++;
        count[n1]++;
        count[n2]++;
    }

    for(int index = 0; index < NPTS; index++)
    {
        gradient[index] /= count[index];
    }
}

void remesh(std::vector<double>& gradient, double lengthScale, double hausd, double grad, std::vector<vec3d>& nodePos, std::vector<vec3i>& elems)
{
#ifdef HAS_MMG
	int NN = NPTS;
	int NF = NCON;
    int NC;

    // we only want to remesh half of the sphere, so here we discard 
    // any nodes that have a z coordinate < 0, and any elements defined
    // with those nodes.
    unordered_map<int, int> newNodeIDs;
    int newNodeID = 1;
    for(int index = 0; index < NN; index++)
    {
        if(ZCOORDS[index] >= 0)
        {
            newNodeIDs[index] = newNodeID;
            newNodeID++;
        }
    }

    unordered_map<int, int> newElemIDs;
    int newElemID = 1;
    for(int index = 0; index < NF; index++)
    {
        if(newNodeIDs.count(CONN1[index]-1) == 0 || newNodeIDs.count(CONN2[index]-1) == 0 || newNodeIDs.count(CONN3[index]-1) == 0)
        {
            continue;
        }

        newElemIDs[index] = newElemID;
        newElemID++;
    }

	// build the MMG mesh
	MMG5_pMesh mmgMesh;
	MMG5_pSol  mmgSol;
	mmgMesh = NULL;
	mmgSol = NULL;
	MMGS_Init_mesh(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh,
		MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	// allocate mesh size
	if (MMGS_Set_meshSize(mmgMesh, newNodeID-1, newElemID-1, 0) != 1)
	{
		assert(false);
		// SetError("Error in MMGS_Set_meshSize");
		// return nullptr;
	}

	// build the MMG mesh
	for (int i = 0; i < NN; ++i)
	{
        if(newNodeIDs.count(i))
        {
            MMGS_Set_vertex(mmgMesh, XCOORDS[i], YCOORDS[i], ZCOORDS[i], 0, newNodeIDs[i]);
        }
	}

	for (int i = 0; i < NF; ++i)
	{
        if(newElemIDs.count(i))
        {
            MMGS_Set_triangle(mmgMesh, newNodeIDs[CONN1[i]-1], newNodeIDs[CONN2[i]-1], newNodeIDs[CONN3[i]-1], 0, newElemIDs[i]);
        }
	}
	
    // Now, we build the "solution", i.e. the target element size.
	// If no elements are selected, we set a homogenous remeshing using the element size parameter.
	// set the "solution", i.e. desired element size
	if (MMGS_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, newNodeID-1, MMG5_Scalar) != 1)
	{
		assert(false);
		// SetError("Error in MMG3D_Set_solSize");
		// return nullptr;
	}

    int n0 = CONN1[0]-1;
    int n1 = CONN2[0]-1;

    vec3d pos0(XCOORDS[n0], YCOORDS[n0], ZCOORDS[n0]);
    vec3d pos1(XCOORDS[n1], YCOORDS[n1], ZCOORDS[n1]);

    double minLength = (pos0 - pos1).Length();
    double maxLength = minLength*lengthScale;

    double min = *std::min_element(gradient.begin(), gradient.end());
    double max = *std::max_element(gradient.begin(), gradient.end());
    double range = max-min;


    for (int k = 0; k < NN; k++) {
        if(newNodeIDs.count(k))
        {
            double val = (maxLength - minLength)*(1-(gradient[k] - min)/range) + minLength;
            MMGS_Set_scalarSol(mmgSol, val, newNodeIDs[k]);
        }
    }

	// set the control parameters
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hmin, minLength);
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hausd, hausd);
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hgrad, grad);

    // prevent MMG from outputing information to stdout
    MMGS_Set_iparameter(mmgMesh, mmgSol, MMGS_IPARAM_verbose, -1);

	// run the mesher
	int ier = MMGS_mmgslib(mmgMesh, mmgSol);

	if (ier == MMG5_STRONGFAILURE) {
		// if (min == 0.0) SetError("Element size cannot be zero.");
		// else SetError("MMG was not able to remesh the mesh.");
		// return nullptr;
	}
	else if (ier == MMG5_LOWFAILURE)
	{
		// SetError("MMG return low failure error");
	}

	// convert back to prv mesh
	// FSMesh* newMesh = new FSMesh();

	// get the new mesh sizes
	MMGS_Get_meshSize(mmgMesh, &NN, &NF, &NC);
	// newMesh->Create(NN, NF);

    nodePos.resize(NN);

	// get the vertex coordinates
	for (int i = 0; i < NN; ++i)
	{
		// FSNode& vi = newMesh->Node(i);
		// vec3d& ri = vi.r;

        double x,y,z;
        int g;
		int isCorner = 0;
		MMGS_Get_vertex(mmgMesh, &x, &y, &z, &g, &isCorner, NULL);

        nodePos[i] = vec3d(x,y,z);
	}

    elems.resize(NF);

    // create elements
	for (int i=0; i<NF; ++i)
	{
        int n0, n1, n2, id;

        MMGS_Get_triangle(mmgMesh, &n0, &n1, &n2, &id, NULL);
		n0--;
		n1--;
		n2--;

        elems[i] = vec3i(n0, n1,n2);
	}

	// Clean up
	MMGS_Free_all(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

    // newMesh->RebuildMesh();

	// return newMesh;

#else
	SetError("This version does not have MMG support");
	return nullptr;
#endif
}

void remeshFull(std::vector<double>& gradient, double lengthScale, double hausd, double grad, std::vector<vec3d>& nodePos, std::vector<vec3i>& elems)
{
#ifdef HAS_MMG
	int NN = NPTS;
	int NF = NCON;
    int NC;

	// build the MMG mesh
	MMG5_pMesh mmgMesh;
	MMG5_pSol  mmgSol;
	mmgMesh = NULL;
	mmgSol = NULL;
	MMGS_Init_mesh(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh,
		MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	// allocate mesh size
	if (MMGS_Set_meshSize(mmgMesh, NN, NF, 0) != 1)
	{
		assert(false);
		// SetError("Error in MMGS_Set_meshSize");
		// return nullptr;
	}

	// build the MMG mesh
	for (int i = 0; i < NN; ++i)
	{
        MMGS_Set_vertex(mmgMesh, XCOORDS[i], YCOORDS[i], ZCOORDS[i], 0, i+1);
	}

	for (int i = 0; i < NF; ++i)
	{
        MMGS_Set_triangle(mmgMesh, CONN1[i], CONN2[i], CONN3[i], 0, i+1);
	}
	
    // Now, we build the "solution", i.e. the target element size.
	// If no elements are selected, we set a homogenous remeshing using the element size parameter.
	// set the "solution", i.e. desired element size
	if (MMGS_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, NN, MMG5_Scalar) != 1)
	{
		assert(false);
		// SetError("Error in MMG3D_Set_solSize");
		// return nullptr;
	}

    int n0 = CONN1[0]-1;
    int n1 = CONN2[0]-1;

    vec3d pos0(XCOORDS[n0], YCOORDS[n0], ZCOORDS[n0]);
    vec3d pos1(XCOORDS[n1], YCOORDS[n1], ZCOORDS[n1]);

    double minLength = (pos0 - pos1).Length();
    double maxLength = minLength*lengthScale;

    double min = *std::min_element(gradient.begin(), gradient.end());
    double max = *std::max_element(gradient.begin(), gradient.end());
    double range = max-min;


    for (int k = 0; k < NN; k++)
    {
        double val = (maxLength - minLength)*(1-(gradient[k] - min)/range) + minLength;
        MMGS_Set_scalarSol(mmgSol, val, k+1);
    }

	// set the control parameters
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hmin, minLength);
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hausd, hausd);
	MMGS_Set_dparameter(mmgMesh, mmgSol, MMGS_DPARAM_hgrad, grad);

    // prevent MMG from outputing information to stdout
    MMGS_Set_iparameter(mmgMesh, mmgSol, MMGS_IPARAM_verbose, -1);

	// run the mesher
	int ier = MMGS_mmgslib(mmgMesh, mmgSol);

	if (ier == MMG5_STRONGFAILURE) {
		// if (min == 0.0) SetError("Element size cannot be zero.");
		// else SetError("MMG was not able to remesh the mesh.");
		// return nullptr;
	}
	else if (ier == MMG5_LOWFAILURE)
	{
		// SetError("MMG return low failure error");
	}

	// convert back to prv mesh
	// FSMesh* newMesh = new FSMesh();

	// get the new mesh sizes
	MMGS_Get_meshSize(mmgMesh, &NN, &NF, &NC);
	// newMesh->Create(NN, NF);

    nodePos.resize(NN);

	// get the vertex coordinates
	for (int i = 0; i < NN; ++i)
	{
		// FSNode& vi = newMesh->Node(i);
		// vec3d& ri = vi.r;

        double x,y,z;
        int g;
		int isCorner = 0;
		MMGS_Get_vertex(mmgMesh, &x, &y, &z, &g, &isCorner, NULL);

        nodePos[i] = vec3d(x,y,z);
	}

    elems.resize(NF);

    // create elements
	for (int i=0; i<NF; ++i)
	{
        int n0, n1, n2, id;

        MMGS_Get_triangle(mmgMesh, &n0, &n1, &n2, &id, NULL);
		n0--;
		n1--;
		n2--;

        elems[i] = vec3i(n0, n1,n2);
	}

	// Clean up
	MMGS_Free_all(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

    // newMesh->RebuildMesh();

	// return newMesh;

#else
	SetError("This version does not have MMG support");
	return nullptr;
#endif
}