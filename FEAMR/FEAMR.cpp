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
#include "stdafx.h"
#include "FEAMR.h"
#include <FECore/FECoreKernel.h>
#include "FEErosionAdaptor.h"
#include "FEHexRefine.h"
#include "FEHexRefine2D.h"
#include "FETetRefine.h"
#include "FEMMGRemesh.h"
#include "FETestRefine.h"
#include "FEVariableCriterion.h"
#include "FEElementSelectionCriterion.h"
#include "FEScaleAdaptorCriterion.h"
#include "FEFilterAdaptorCriterion.h"
#include "FEDomainErrorCriterion.h"
#include "FEElementDataCriterion.h"

//-----------------------------------------------------------------------------
void FEAMR::InitModule()
{
// mesh adaptors
REGISTER_FECORE_CLASS(FEErosionAdaptor, "erosion");
REGISTER_FECORE_CLASS(FEHexRefine     , "hex_refine");
REGISTER_FECORE_CLASS(FEHexRefine2D   , "hex_refine2d");
REGISTER_FECORE_CLASS(FETetRefine     , "tet_refine");
REGISTER_FECORE_CLASS(FEMMGRemesh     , "mmg_remesh");
REGISTER_FECORE_CLASS(FETestRefine    , "test_refine");

// adaptor criteria
REGISTER_FECORE_CLASS(FEVariableCriterion        , "max_variable");
REGISTER_FECORE_CLASS(FEElementSelectionCriterion, "element_selection");
REGISTER_FECORE_CLASS(FEScaleAdaptorCriterion    , "math");
REGISTER_FECORE_CLASS(FEMinMaxFilterAdaptorCriterion, "min-max filter");
REGISTER_FECORE_CLASS(FEDomainErrorCriterion, "relative error");
REGISTER_FECORE_CLASS(FEElementDataCriterion, "element data");
}
