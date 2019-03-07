#pragma once
#include <FECore/FEModel.h>

#define feLog(...) GetFEModel()->Logf(0, __VA_ARGS__)
#define feLogWarning(...) GetFEModel()->Logf(1, __VA_ARGS__)
#define feLogError(...) GetFEModel()->Logf(2, __VA_ARGS__)

#define feLogEx(fem, ...) fem->Logf(0, __VA_ARGS__)
#define feLogWarningEx(fem, ...) fem->Logf(1, __VA_ARGS__)
#define feLogErrorEx(fem, ...) fem->Logf(2, __VA_ARGS__)
