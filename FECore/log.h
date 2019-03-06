#pragma once
#include <FECore/FEModel.h>

#define feLog(msg, ...) GetFEModel()->Logf(0, msg, __VA_ARGS__)
#define feLogWarning(msg, ...) GetFEModel()->Logf(1, msg, __VA_ARGS__)
#define feLogError(msg, ...) GetFEModel()->Logf(2, msg, __VA_ARGS__)

#define feLogEx(fem, msg, ...) fem->Logf(0, msg, __VA_ARGS__)
#define feLogWarningEx(fem, msg, ...) fem->Logf(1, msg, __VA_ARGS__)
#define feLogErrorEx(fem, msg, ...) fem->Logf(2, msg, __VA_ARGS__)
