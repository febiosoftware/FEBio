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
        theta[index] = acos(zCoords[index]);
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