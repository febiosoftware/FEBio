#include "vm.h"
#include "../febcode/module_vec3.h"
#include "../febcode/module_vec2.h"
#include "../febcode/module_math.h"
#include <iostream>
#include <iomanip>

using namespace febcode;

// helper function for integer exponentiation (used in EXP_INT)
static int ipow(int base, int exp)
{
	int result = 1;
	while (exp > 0)
	{
		if (exp & 1)
			result *= base;
		base *= base;
		exp >>= 1;
	}
	return result;
}

const char* IPToString(uint8_t ip)
{
	switch ((OpCode)ip)
	{
	case OpCode::PUSH_VOID     : 
	case OpCode::PUSH_BOOL     :
	case OpCode::PUSH_INT      :
	case OpCode::PUSH_DOUBLE   : return "MOV ";

	case OpCode::GET_GLOBAL_BOOL  :
	case OpCode::GET_GLOBAL_INT   :
	case OpCode::GET_GLOBAL_DOUBLE:
	case OpCode::GET_GLOBAL_VEC2  :
	case OpCode::GET_GLOBAL_VEC3  : 
	case OpCode::GET_GLOBAL_ARRAY : 
	case OpCode::GET_GLOBAL_STRUCT: return "GTG ";

	case OpCode::SET_GLOBAL_BOOL  :
	case OpCode::SET_GLOBAL_INT   :
	case OpCode::SET_GLOBAL_DOUBLE:
	case OpCode::SET_GLOBAL_VEC2  :
	case OpCode::SET_GLOBAL_VEC3  :
	case OpCode::SET_GLOBAL_ARRAY :
	case OpCode::SET_GLOBAL_STRUCT: return "STG ";

	case OpCode::STORE_BOOL  :
	case OpCode::STORE_INT   :
	case OpCode::STORE_DOUBLE:
	case OpCode::STORE_VEC2  :
	case OpCode::STORE_VEC3  :
	case OpCode::STORE_ARRAY :
	case OpCode::STORE_STRUCT:
		return "STRE";

	case OpCode::GET_GLOBAL_REF: return "GREF";
	case OpCode::GET_LOCAL_REF : return "LREF";
	case OpCode::GET_INDEX_REF : return "IREF";
	case OpCode::GET_MEMBER_REF: return "MREF";

	case OpCode::GET_LOCAL_BOOL  :
	case OpCode::GET_LOCAL_INT   :
	case OpCode::GET_LOCAL_DOUBLE:
	case OpCode::GET_LOCAL_VEC2  :
	case OpCode::GET_LOCAL_VEC3  :
	case OpCode::GET_LOCAL_ARRAY :
	case OpCode::GET_LOCAL_STRUCT: return "GETL";

	case OpCode::ADD_INT       : return "ADDI";
	case OpCode::ADD_DOUBLE    : return "ADDF";
	case OpCode::SUB_INT       : return "SUBI";
	case OpCode::SUB_DOUBLE    : return "SUBF";
	case OpCode::MUL_INT       : return "MULI";
	case OpCode::MUL_DOUBLE    : return "MULF";
	case OpCode::DIV_INT       : return "DIVI";
	case OpCode::DIV_DOUBLE    : return "DIVF";
	case OpCode::EXP_INT       : return "EXPI";
	case OpCode::EXP_DOUBLE    : return "EXPF";

	case OpCode::EQUAL_BOOL  :
	case OpCode::EQUAL_INT   :
	case OpCode::EQUAL_DOUBLE: return "EQ  ";

	case OpCode::NEQ_BOOL  :
	case OpCode::NEQ_INT   :
	case OpCode::NEQ_DOUBLE: return "NEQ ";

	case OpCode::GT_INT        : return "GTI ";
	case OpCode::GT_DOUBLE     : return "GTF ";
	case OpCode::LT_INT        : return "LTI ";
	case OpCode::LT_DOUBLE     : return "LTF ";
	case OpCode::GE_INT        : return "GEI ";
	case OpCode::GE_DOUBLE     : return "GEF ";
	case OpCode::LE_INT        : return "LEI ";
	case OpCode::LE_DOUBLE     : return "LEF ";
	case OpCode::NEG_INT       : return "NEGI";
	case OpCode::NEG_DOUBLE    : return "NEGF";
	case OpCode::CREATE_VEC2   : return "VEC2";
	case OpCode::GET_VEC2_X    : return "GV2X";
	case OpCode::GET_VEC2_Y    : return "GV2Y";
	case OpCode::GET_VEC2_SWIZZLE: return "G2SW";
	case OpCode::ADD_VEC2      : return "ADD2";
	case OpCode::SUB_VEC2      : return "SUB2";
	case OpCode::DOT_VEC2      : return "DOT2";
	case OpCode::MUL_VEC2_DOUBLE: return "ML2F";
	case OpCode::MUL_DOUBLE_VEC2: return "MLF2";
	case OpCode::NEG_VEC2      : return "NEG2";
	case OpCode::CREATE_VEC3   : return "VEC3";
	case OpCode::GET_VEC3_X    : return "GV3X";
	case OpCode::GET_VEC3_Y    : return "GV3Y";
	case OpCode::GET_VEC3_Z    : return "GV3Z";
	case OpCode::GET_VEC3_SWIZZLE: return "G3SW";
	case OpCode::ADD_VEC3      : return "ADD3";
	case OpCode::SUB_VEC3      : return "SUB3";
	case OpCode::DOT_VEC3      : return "DOT3";
	case OpCode::MUL_VEC3_DOUBLE: return "ML3F";
	case OpCode::MUL_DOUBLE_VEC3: return "MLF3";
	case OpCode::NEG_VEC3      : return "NEG3";
	case OpCode::NOT           : return "NOT ";
	case OpCode::CREATE_STRUCT : return "STRC";
	case OpCode::COPY_STRUCT   : return "CPYS";

	case OpCode::GET_PROPERTY_BOOL:
	case OpCode::GET_PROPERTY_INT:
	case OpCode::GET_PROPERTY_DOUBLE:
	case OpCode::GET_PROPERTY_VEC2:
	case OpCode::GET_PROPERTY_VEC3:
	case OpCode::GET_PROPERTY_ARRAY:
	case OpCode::GET_PROPERTY_STRUCT:
		return "GETP";

	case OpCode::CREATE_ARRAY  : return "ARR ";
	case OpCode::COPY_ARRAY    : return "CPYA";

	case OpCode::GET_INDEX_BOOL  :
	case OpCode::GET_INDEX_INT   :
	case OpCode::GET_INDEX_DOUBLE:
	case OpCode::GET_INDEX_VEC2  :
	case OpCode::GET_INDEX_VEC3  :
	case OpCode::GET_INDEX_ARRAY :
	case OpCode::GET_INDEX_STRUCT:
		return "GETI";

	case OpCode::JUMP          : return "JMP ";
	case OpCode::JUMP_IF_FALSE : return "JMPF";
	case OpCode::JUMP_IF_TRUE  : return "JMPT";
	case OpCode::LOOP          : return "LOOP";
	case OpCode::CALL          : return "CALL";
	case OpCode::CALL_BINARY   : return "BINO";
	case OpCode::POP           : return "POP ";
	case OpCode::GET_VEC2_X_REF: return "RV2X";
	case OpCode::GET_VEC2_Y_REF: return "RV2Y";
	case OpCode::GET_VEC3_X_REF: return "RV3X";
	case OpCode::GET_VEC3_Y_REF: return "RV3Y";
	case OpCode::GET_VEC3_Z_REF: return "RV3Z";

	case OpCode::RETURN_VOID: 
	case OpCode::RETURN_BOOL: 
	case OpCode::RETURN_INT: 
	case OpCode::RETURN_DOUBLE: 
	case OpCode::RETURN_VEC2: 
	case OpCode::RETURN_VEC3: 
	case OpCode::RETURN_ARRAY: 
	case OpCode::RETURN_STRUCT: 
		return "RET ";
	}
	return "(UNKNOWN)";
}

void printStack(const std::vector<Value>& stack, int numGlobals, int stackSize)
{
	std::cout << "Stack: [";
	for (size_t i = 0; i < numGlobals; i++)
	{
		std::cout << ValueToString(stack[i]);
		if (i < numGlobals - 1)
			std::cout << ",";
	}
	std::cout << "|";
	for (size_t i = numGlobals; i < stackSize; i++)
	{
		std::cout << ValueToString(stack[i]);
		if (i < stack.size() - 1)
			std::cout << ",";
	}
	for (size_t i= stackSize; i < stack.size(); i++)
	{
		std::cout << "_";
		if (i < stack.size() - 1)
			std::cout << ",";
	}
	std::cout << "]" << std::endl;
}

Value VM::execute()
{
	if (frameCount == 0)
		throw std::runtime_error("No function to execute");

	const size_t instructions = m_program->code.size();

	while (currentFrame().ip < instructions)
	{
		uint8_t instruction = readByte();

		if (m_debug)
		{
			std::cout << "IP: " << std::setw(4) << currentFrame().ip - 1;
			std::cout << " | Executing: " << IPToString(instruction) << " | ";
		}

		switch ((OpCode)instruction)
		{
		case OpCode::PUSH_VOID:
		{
			pushVoid();
			break;
		}

		case OpCode::PUSH_BOOL:
		{
			uint8_t idx = readByte();
			pushBool(m_program->constants[idx].b);
			break;
		}

		case OpCode::PUSH_INT:
		{
			uint8_t idx = readByte();
			pushInt(m_program->constants[idx].i);
			break;
		}

		case OpCode::PUSH_DOUBLE:
		{
			uint8_t idx = readByte();
			pushDouble(m_program->constants[idx].d);
			break;
		}

		case OpCode::GET_GLOBAL_BOOL:
		{
			uint8_t slot = readByte();
			pushBool(getBoolAt(slot));
			break;
		}

		case OpCode::GET_GLOBAL_INT:
		{
			uint8_t slot = readByte();
			pushInt(getIntAt(slot));
			break;
		}

		case OpCode::GET_GLOBAL_DOUBLE:
		{
			uint8_t slot = readByte();
			pushDouble(getDoubleAt(slot));
			break;
		}

		case OpCode::GET_GLOBAL_VEC2:
		{
			uint8_t slot = readByte();
			pushVec2(getVec2At(slot));
			break;
		}

		case OpCode::GET_GLOBAL_VEC3:
		{
			uint8_t slot = readByte();
			pushVec3(getVec3At(slot));
			break;
		}

		case OpCode::GET_GLOBAL_ARRAY:
		{
			uint8_t slot = readByte();
			const ArrayValuePtr& arr = m_stack[slot].arrayValue;
			pushArray(arr);
			break;
		}

		case OpCode::GET_GLOBAL_STRUCT:
		{
			uint8_t slot = readByte();
			const StructValuePtr& obj = m_stack[slot].structValue;
			pushStruct(obj);
			break;
		}

		case OpCode::GET_GLOBAL_REF:
		{
			uint8_t slot = readByte();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &m_stack[slot];
			ref.ref.type = RefType::Value;

			push(ref);
			break;
		}

		case OpCode::GET_LOCAL_REF:
		{
			uint8_t slot = readByte();

			CallFrame& frame = currentFrame();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &m_stack[frame.base + slot];
			ref.ref.type = RefType::Value;

			push(ref);
			break;
		}

		case OpCode::GET_MEMBER_REF:
		{
			uint8_t memberIndex = readByte();

			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			StructValue* s = slot->structValue.get();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &s->fields[memberIndex];
			ref.ref.type = RefType::Value;

			push(ref);
			break;
		}

		case OpCode::GET_INDEX_REF:
		{
			const Value& indexVal = pop();
			int index = getInt(indexVal);

			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			ArrayValue* s = slot->arrayValue.get();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &s->elements[index];
			ref.ref.type = RefType::Value;

			push(ref);
			break;
		}

		case OpCode::GET_VEC2_X_REF:
		{
			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			vec2* v = &slot->vec2Value;

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &v->x;
			ref.ref.type = RefType::Double;

			push(ref);
			break;
		}
		case OpCode::GET_VEC2_Y_REF:
		{
			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			vec2* v = &slot->vec2Value;

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &v->y;
			ref.ref.type = RefType::Double;

			push(ref);
			break;
		}

		case OpCode::GET_VEC3_X_REF:
		{
			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			vec3* v = &slot->vec3Value;

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &v->x;
			ref.ref.type = RefType::Double;

			push(ref);
			break;
		}

		case OpCode::GET_VEC3_Y_REF:
		{
			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			vec3* v = &slot->vec3Value;

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &v->y;
			ref.ref.type = RefType::Double;

			push(ref);
			break;
		}

		case OpCode::GET_VEC3_Z_REF:
		{
			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			vec3* v = &slot->vec3Value;

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &v->z;
			ref.ref.type = RefType::Double;

			push(ref);
			break;
		}

		case OpCode::STORE_BOOL:
		{
			bool b = popBool();
			const Ref& r = popRef();
			*(Value*)r.ptr = b;
			pushBool(b);
			break;
		}

		case OpCode::STORE_INT :
		{
			int n = popInt();
			const Ref& r = popRef();
			*(Value*)r.ptr = n;
			pushInt(n);
			break;
		}

		case OpCode::STORE_DOUBLE:
		{
			double d = popDouble();
			const Ref& r = popRef();
			if (r.type == RefType::Double)
				*(double*)r.ptr = d;
			else
				*(Value*)r.ptr = d;
			pushDouble(d);
			break;
		}

		case OpCode::STORE_VEC2:
		{
			vec2 v = popVec2();
			const Ref& r = popRef();
			*(Value*)r.ptr = v;
			pushVec2(v);
			break;
		}

		case OpCode::STORE_VEC3:
		{
			vec3 v = popVec3();
			const Ref& r = popRef();
			*(Value*)r.ptr = v;
			pushVec3(v);
			break;
		}

		case OpCode::STORE_ARRAY:
		case OpCode::STORE_STRUCT:
		{
			const Value& value = pop();
			const Value& ref = pop();
			const Ref& r = ref.ref;
			*(Value*)r.ptr = value;
			push(value);
			break;
		}

		case OpCode::SET_GLOBAL_BOOL:
		{
			uint8_t slot = readByte();
			setBoolAt(slot, peekBool());
			break;
		}

		case OpCode::SET_GLOBAL_INT:
		{
			uint8_t slot = readByte();
			setIntAt(slot, peekInt());
			break;
		}

		case OpCode::SET_GLOBAL_DOUBLE:
		{
			uint8_t slot = readByte();
			setDoubleAt(slot, peekDouble());
			break;
		}

		case OpCode::SET_GLOBAL_VEC2:
		{
			uint8_t slot = readByte();
			setVec2At(slot, peekVec2());
			break;
		}

		case OpCode::SET_GLOBAL_VEC3:
		{
			uint8_t slot = readByte();
			setVec3At(slot, peekVec3());
			break;
		}

		case OpCode::SET_GLOBAL_ARRAY:
		{
			uint8_t slot = readByte();
			m_stack[slot] = peek();
			break;
		}

		case OpCode::SET_GLOBAL_STRUCT:
		{
			uint8_t slot = readByte();
			m_stack[slot] = peek();
			break;
		}

		case OpCode::GET_LOCAL_BOOL:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushBool(m_stack[frame.base + slot].b);
			break;
		}

		case OpCode::GET_LOCAL_INT:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushInt(m_stack[frame.base + slot].i);
			break;
		}

		case OpCode::GET_LOCAL_DOUBLE:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushDouble(m_stack[frame.base + slot].d);
			break;
		}

		case OpCode::GET_LOCAL_VEC2:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushVec2(m_stack[frame.base + slot].vec2Value);
			break;
		}

		case OpCode::GET_LOCAL_VEC3:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushVec3(m_stack[frame.base + slot].vec3Value);
			break;
		}

		case OpCode::GET_LOCAL_ARRAY:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			push(m_stack[frame.base + slot]);
			break;
		}

		case OpCode::GET_LOCAL_STRUCT:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			push(m_stack[frame.base + slot]);
			break;
		}

		// Integer operators
		case OpCode::NEG_INT:
		{
			int a = popInt();
			pushInt(-a);
			break;
		}
		case OpCode::ADD_INT:
		{
			int b = popInt();
			int a = popInt();
			pushInt(a + b);
			break;
		}
		case OpCode::SUB_INT:
		{
			int b = popInt();
			int a = popInt();
			pushInt(a - b);
			break;
		}
		case OpCode::MUL_INT:
		{
			int b = popInt();
			int a = popInt();
			pushInt(a * b);
			break;
		}
		case OpCode::DIV_INT:
		{
			int b = popInt();
			int a = popInt();
			if (b == 0)
				throw std::runtime_error("division by zero.");
			pushInt(a / b);
			break;
		}
		case OpCode::EXP_INT:
		{
			int b = popInt();
			int a = popInt();

			if (b < 0)
				throw std::runtime_error("Negative exponent not supported for integers.");

			pushInt(ipow(a, b));
			break;
		}
		case OpCode::GT_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a > b);
			break;
		}
		case OpCode::LT_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a < b);
			break;
		}
		case OpCode::GE_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a >= b);
			break;
		}
		case OpCode::LE_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a <= b);
			break;
		}

		// Double operators
		case OpCode::NEG_DOUBLE:
		{
			double a = popDouble();
			pushDouble(-a);
			break;
		}

		case OpCode::ADD_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushDouble(a + b);
			break;
		}
		case OpCode::SUB_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushDouble(a - b);
			break;
		}
		case OpCode::MUL_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushDouble(a * b);
			break;
		}
		case OpCode::DIV_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			if (b == 0.0)
				throw std::runtime_error("division by zero.");
			pushDouble(a / b);
			break;
		}
		case OpCode::EXP_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushDouble(std::pow(a, b));
			break;
		}
		case OpCode::GT_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a > b);
			break;
		}
		case OpCode::LT_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a < b);
			break;
		}
		case OpCode::GE_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a >= b);
			break;
		}
		case OpCode::LE_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a <= b);
			break;
		}

		case OpCode::CREATE_VEC2:
		{
			double y = popDouble();
			double x = popDouble();
			pushVec2(vec2(x, y));
			break;
		}
		case OpCode::ADD_VEC2:
		{
			vec2 b = popVec2();
			vec2 a = popVec2();
			pushVec2(a + b);
			break;
		}
		case OpCode::SUB_VEC2:
		{
			vec2 b = popVec2();
			vec2 a = popVec2();
			pushVec2(a - b);
			break;
		}
		case OpCode::DOT_VEC2:
		{
			vec2 b = popVec2();
			vec2 a = popVec2();
			pushDouble(a * b);
			break;
		}
		case OpCode::MUL_VEC2_DOUBLE:
		{
			double scalar = popDouble();
			vec2 vec = popVec2();
			pushVec2(vec * scalar);
			break;
		}
		case OpCode::MUL_DOUBLE_VEC2:
		{
			vec2 vec = popVec2();
			double scalar = popDouble();
			pushVec2(vec * scalar);
			break;
		}
		case OpCode::GET_VEC2_X:
		{
			vec2 vec = popVec2();
			pushDouble(vec.x);
			break;
		}
		case OpCode::GET_VEC2_Y:
		{
			vec2 vec = popVec2();
			pushDouble(vec.y);
			break;
		}
		case OpCode::NEG_VEC2:
		{
			vec2 v = popVec2();
			pushVec2(vec2(-v.x, -v.y));
			break;
		}
		case OpCode::GET_VEC2_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			vec2 v = popVec2();
			double c[4] = { v.x, v.y, 0.0, 0.0 };

			if (size == 2)
			{
				vec2 result;
				result.y = c[mask & 0b11];
				result.x = c[(mask >> 2) & 0b11];
				pushVec2(result);
			}
			else if (size == 3)
			{
				vec3 result;
				result.z = c[mask & 0b11];
				result.y = c[(mask >> 2) & 0b11];
				result.x = c[(mask >> 4) & 0b11];
				pushVec3(result);
			}
			break;
		}
		case OpCode::CREATE_VEC3:
		{
			double z = popDouble();
			double y = popDouble();
			double x = popDouble();
			pushVec3(vec3(x, y, z));
			break;
		}
		case OpCode::ADD_VEC3:
		{
			vec3 b = popVec3();
			vec3 a = popVec3();
			pushVec3(a + b);
			break;
		}
		case OpCode::SUB_VEC3:
		{
			vec3 b = popVec3();
			vec3 a = popVec3();
			pushVec3(a - b);
			break;
		}
		case OpCode::DOT_VEC3:
		{
			vec3 b = popVec3();
			vec3 a = popVec3();
			pushDouble(a * b);
			break;
		}
		case OpCode::MUL_VEC3_DOUBLE:
		{
			double scalar = popDouble();
			vec3 vec = popVec3();
			pushVec3(vec * scalar);
			break;
		}
		case OpCode::MUL_DOUBLE_VEC3:
		{
			vec3 vec = popVec3();
			double scalar = popDouble();
			pushVec3(vec * scalar);
			break;
		}
		case OpCode::GET_VEC3_X:
		{
			vec3 vec = popVec3();
			pushDouble(vec.x);
			break;
		}
		case OpCode::GET_VEC3_Y:
		{
			vec3 vec = popVec3();
			pushDouble(vec.y);
			break;
		}
		case OpCode::GET_VEC3_Z:
		{
			vec3 vec = popVec3();
			pushDouble(vec.z);
			break;
		}
		case OpCode::GET_VEC3_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			vec3 v = popVec3();
			double c[4] = { v.x, v.y, v.z, 0.0 };

			if (size == 2)
			{
				vec2 result;
				result.y = c[mask & 0b11];
				result.x = c[(mask >> 2) & 0b11];
				pushVec2(result);
			}
			else if (size == 3)
			{
				vec3 result;
				result.z = c[mask & 0b11];
				result.y = c[(mask >> 2) & 0b11];
				result.x = c[(mask >> 4) & 0b11];
				pushVec3(result);
			}
			break;
		}

		case OpCode::NEG_VEC3:
		{
			vec3 v = popVec3();
			pushVec3(vec3(-v.x, -v.y, -v.z));
			break;
		}

		// Logical operators
		case OpCode::NOT:
		{
			bool b = popBool();
			pushBool(!b);
			break;
		}

		case OpCode::EQUAL_BOOL:
		{
			bool b = popBool();
			bool a = popBool();
			pushBool(a == b);
			break;
		}
		case OpCode::EQUAL_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a == b);
			break;
		}
		case OpCode::EQUAL_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a == b);
			break;
		}

		case OpCode::NEQ_BOOL:
		{
			bool b = popBool();
			bool a = popBool();
			pushBool(a != b);
			break;
		}

		case OpCode::NEQ_INT:
		{
			int b = popInt();
			int a = popInt();
			pushBool(a != b);
			break;
		}

		case OpCode::NEQ_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();
			pushBool(a != b);
			break;
		}

		case OpCode::CREATE_STRUCT:
		{
			uint8_t typeIndex = readByte();
			Type type = m_program->types.getStructType(typeIndex);
			int fieldCount = (int)type->fields.size();

			auto obj = createStruct(type);
			for (int i = fieldCount - 1; i >= 0; --i)
			{
				obj->fields[i] = pop();
			}

			push(obj);
			break;
		}
		case OpCode::COPY_STRUCT:
		{
			auto src = pop();

			if (!isStruct(src))
				throw std::runtime_error("COPY_STRUCT operand must be a struct");

			const StructValue& srcObj = getStruct(src);
			auto obj = copyStruct(srcObj);

			push(obj);
			break;
		}
		case OpCode::GET_PROPERTY_BOOL:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			pushBool(obj.fields[slot].b);
			break;
		}

		case OpCode::GET_PROPERTY_INT:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			pushInt(obj.fields[slot].i);
			break;
		}

		case OpCode::GET_PROPERTY_DOUBLE:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			pushDouble(obj.fields[slot].d);
			break;
		}

		case OpCode::GET_PROPERTY_VEC2:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			pushVec2(obj.fields[slot].vec2Value);
			break;
		}

		case OpCode::GET_PROPERTY_VEC3:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			pushVec3(obj.fields[slot].vec3Value);
			break;
		}

		case OpCode::GET_PROPERTY_ARRAY:
		case OpCode::GET_PROPERTY_STRUCT:
		{
			const StructValue& obj = popStruct();
			uint8_t slot = readByte();
			push(obj.fields[slot]);
			break;
		}

		case OpCode::CREATE_ARRAY: {
			uint8_t count = readByte();
			auto arr = std::make_shared<ArrayValue>();

			arr->elements.resize(count);

			for (int i = count - 1; i >= 0; --i) {
				arr->elements[i] = pop();
			}

			Type elemType = m_program->types.getBuiltinType(arr->elements[0]);
			arr->type = m_program->types.getArrayType(elemType, count);

			push(arr);
			break;
		}
		case OpCode::COPY_ARRAY: {
			const ArrayValue& arr = popArray();
			auto obj = copyArray(arr);
			push(obj);
			break;
		}
		case OpCode::GET_INDEX_BOOL:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushBool(arr.elements[index].b);
			break;
		}

		case OpCode::GET_INDEX_INT:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushInt(arr.elements[index].i);
			break;
		}

		case OpCode::GET_INDEX_DOUBLE:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushDouble(arr.elements[index].d);
			break;
		}

		case OpCode::GET_INDEX_VEC2:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushVec2(arr.elements[index].vec2Value);
			break;
		}

		case OpCode::GET_INDEX_VEC3:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushVec3(arr.elements[index].vec3Value);
			break;
		}

		case OpCode::GET_INDEX_ARRAY:
		case OpCode::GET_INDEX_STRUCT:
		{
			int index = popInt();
			const ArrayValue& arr = popArray();

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			push(arr.elements[index]);
			break;
		}

		case OpCode::JUMP:
		{
			uint16_t offset = readUint16();
			currentFrame().ip += offset;
			break;
		}

		case OpCode::JUMP_IF_FALSE:
		{
			uint16_t offset = readUint16();
			if (!isTruthy(peek()))
				currentFrame().ip += offset;
			break;
		}

		case OpCode::JUMP_IF_TRUE:
		{
			uint16_t offset = readUint16();
			if (isTruthy(peek()))
				currentFrame().ip += offset;
			break;
		}

		case OpCode::LOOP:
		{
			uint16_t offset = readUint16();
			currentFrame().ip -= offset;
			break;
		}

		case OpCode::POP:
		{
			pop();
			break;
		}

		case OpCode::CALL:
		{
			uint8_t fnIndex = readByte();
			uint8_t args = readByte();
			callFunction(fnIndex, args);
			break;
		}

		case OpCode::CALL_BINARY:
		{
			uint8_t opIndex = readByte();
			callBinaryOperator(opIndex);
			break;
		}
		
		case OpCode::RETURN_VOID:
		case OpCode::RETURN_BOOL:
		case OpCode::RETURN_INT:
		case OpCode::RETURN_DOUBLE:
		case OpCode::RETURN_VEC2:
		case OpCode::RETURN_VEC3:
		case OpCode::RETURN_ARRAY:
		case OpCode::RETURN_STRUCT:
		{
			TypeKind returnType = (TypeKind)(instruction - (uint8_t)OpCode::RETURN_VOID);

			Value result;
			if (stackTop > globalCount)
			{
				switch (returnType)
				{
				case TypeKind::Void  : pop(); break;
				case TypeKind::Bool  : result = popBool  (); break;
				case TypeKind::Int   : result = popInt   (); break;
				case TypeKind::Double: result = popDouble(); break;
				case TypeKind::Vec2  : result = popVec2  (); break;
				case TypeKind::Vec3  : result = popVec3  (); break;
				default:
					result = pop();
				}
			}

			CallFrame frame = m_frames[--frameCount];
			stackTop = frame.base;

			switch (returnType)
			{
			case TypeKind::Void  : push      (result); break;
			case TypeKind::Bool  : pushBool  (result.b); break;
			case TypeKind::Int   : pushInt   (result.i); break;
			case TypeKind::Double: pushDouble(result.d); break;
			case TypeKind::Vec2  : pushVec2  (result.vec2Value); break;
			case TypeKind::Vec3  : pushVec3  (result.vec3Value); break;
			default:
				push(result);
			}

			if (m_debug)
			{
				printStack(m_stack, globalCount, (int)stackTop);
			}

			if (frameCount == 0)
				return result;

			break;
		}

		default:
			throw std::runtime_error("Unknown opcode");
		}

		if (m_debug && (instruction < (uint8_t)OpCode::RETURN_VOID))
		{
			printStack(m_stack, globalCount, (int)stackTop);
		}
	}

	throw std::runtime_error("Unexpected end of code.");
}

void VM::callFunction(int fnIndex, int argCount)
{
	const FunctionInfo& fn = m_program->functions[fnIndex];

	if (fn.isNative)
	{
		Value result = fn.fnc(&m_stack[stackTop - argCount], argCount);
		stackTop -= argCount;
		push(result);
		return;
	}

	if (argCount != fn.args.size())
		throw std::runtime_error("Arity mismatch in call to " + fn.name);

	if (frameCount >= MAX_CALL_DEPTH)
		throw std::runtime_error("Stack overflow: too many nested function calls.");

	CallFrame frame;
	frame.functionIndex = fnIndex;
	frame.ip = fn.entry;
	frame.base = stackTop - argCount;

	m_frames[frameCount++] = frame;
}

void VM::callBinaryOperator(int opIndex)
{
	const BinaryOperatorInfo& fn = m_program->operators[opIndex];

	Value lv = m_stack[stackTop - 2];
	Value rv = m_stack[stackTop - 1];
	stackTop -= 2; // pop the operands

	Value result = fn.fnc(lv, rv);
	push(result);
}

Value febcode::runScript(const std::string& script, bool initModules)
{
	Program prg;

	if (initModules)
	{
		MathModule mathModule;
		mathModule.Register(prg);

		Vec2Module vec2Module;
		vec2Module.Register(prg);

		Vec3Module vec3Module;
		vec3Module.Register(prg);
	}

	CompileSource(prg, script);
	VM vm(prg);
	return vm.run();
}
