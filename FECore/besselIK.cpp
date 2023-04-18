/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "stdafx.h"
#include "gamma.h"
#include <limits>
#include <math.h>


#include "besselIK.h"

// modified Bessel function of the first kind I0(x) (real x)
double i0(double x)
{
    double ax,ans;
    double y;
    if ((ax=fabs(x)) < 3.75) {
        y=x/3.75;
        y*=y;
        ans=1.0+y*(3.5156229+y*(3.0899424+y*
                                (1.2067492+y*
                                 (0.2659732+y*(0.360768e-1+y*0.45813e-2))
                                 )
                                )
                   );
    }
    else {
        y=3.75/ax;
        ans=(exp(ax)/sqrt(ax))*(0.39894228+y*
                                (0.1328592e-1+y*
                                 (0.225319e-2+y*
                                  (-0.157565e-2+y*(0.916281e-2 +y*
                                                   (-0.2057706e-1+y*
                                                    (0.2635537e-1+y*
                                                     (-0.1647633e-1 +y*0.392377e-2)
                                                     )
                                                    )
                                                   )
                                   )
                                  )
                                 )
                                );
    }
    return ans;
}

// modified Bessel function of the first kind I1(x) (real x)
double i1(double x)
{
    double ax,ans;
    double y;
    if ((ax=fabs(x)) < 3.75) {
        y=x/3.75;
        y*=y;
        ans=ax*(0.5+y*(0.87890594+y*
                       (0.51498869+y*(0.15084934+y*
                                      (0.2658733e-1+y*
                                       (0.301532e-2+y*0.32411e-3)
                                       )
                                      )
                        )
                       )
                );
    }
    else {
        y=3.75/ax;
        ans=0.2282967e-1+y*(-0.2895312e-1+y*
                            (0.1787654e-1-y*0.420059e-2));
        ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*
                                           (0.163801e-2+y*(-0.1031555e-1+y*ans)
                                            )
                                           )
                          );
        ans *= (exp(ax)/sqrt(ax));
    }
    return (x < 0.0) ? -ans : ans;
}

// modified Bessel function of the second kind K0(x) (real x)
double k0(double x)
{
    double y,ans;
    if (x <= 2.0) {
        y=x*x/4.0;
        ans=(-log(x/2.0)*i0(x))+(-0.57721566+y*
                                 (0.42278420 +y*
                                  (0.23069756+y*
                                   (0.3488590e-1+y*
                                    (0.262698e-2 +y*
                                     (0.10750e-3+y*0.74e-5)
                                     )
                                    )
                                   )
                                  )
                                 );
    }
    else {
        y=2.0/x;
        ans=(exp(-x)/sqrt(x))*(1.25331414+y*
                               (-0.7832358e-1 +y*
                                (0.2189568e-1+y*(-0.1062446e-1+y*
                                                 (0.587872e-2 +y*
                                                  (-0.251540e-2+y*0.53208e-3)
                                                  )
                                                 )
                                 )
                                )
                               );
    }
    return ans;
}

// modified Bessel function of the second kind K1(x) (real x)
double k1(double x)
{
    double y,ans;
    if (x <= 2.0) {
        y=x*x/4.0;
        ans=(log(x/2.0)*i1(x))+(1.0/x)*(1.0+y*
                                        (0.15443144 +y*
                                         (-0.67278579+y*
                                          (-0.18156897+y*
                                           (-0.1919402e-1 +y*
                                            (-0.110404e-2+y*
                                             (-0.4686e-4)
                                             )
                                            )
                                           )
                                          )
                                         )
                                        );
    }
    else {
        y=2.0/x;
        ans=(exp(-x)/sqrt(x))*(1.25331414+y*
                               (0.23498619 +y*
                                (-0.3655620e-1+y*
                                 (0.1504268e-1+y*
                                  (-0.780353e-2 +y*
                                   (0.325614e-2+y*
                                    (-0.68245e-3)
                                    )
                                   )
                                  )
                                 )
                                )
                               );
    }
    return ans;
}

