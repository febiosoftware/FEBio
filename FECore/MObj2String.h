#pragma once
#include "MathObject.h"
#include <string>
using namespace std;

// this class converts a MathObject to a string
class MObj2String
{
public:
	string Convert(const MathObject& o);

protected:
	string Convert(const MItem* pi);

	string Constant (const MConstant* pc);
	string Fraction (const MFraction* pc);
	string NamedCt  (const MNamedCt*  pc);
	string Variable (const MVarRef*   pv);
	string OpNeg    (const MNeg*      po);
	string OpAdd    (const MAdd*      po);
	string OpSub    (const MSub*      po);
	string OpMul    (const MMul*      po);
	string OpDiv    (const MDiv*      po);
	string OpPow    (const MPow*      po);
	string OpEqual  (const MEquation* po);
	string OpFnc1D  (const MFunc1D*   po);
	string OpFnc2D  (const MFunc2D*   po);
	string OpFncND  (const MFuncND*   po);
	string OpSFnc   (const MSFuncND*  po);
};
