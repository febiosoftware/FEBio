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
		return "STOR";

	case OpCode::GET_GLOBAL_REF: 
		return "GREF";

	case OpCode::GET_LOCAL_REF       :
		return "LREF";

	case OpCode::GET_INDEX_REF       : 
	case OpCode::GET_INDEX_REF_BOOL  : 
	case OpCode::GET_INDEX_REF_INT   : 
	case OpCode::GET_INDEX_REF_DOUBLE: 
	case OpCode::GET_INDEX_REF_VEC2  : 
	case OpCode::GET_INDEX_REF_VEC3  : 
		return "IREF";

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
	case OpCode::SQR_DOUBLE    : return "SQR ";
	case OpCode::SQRT_DOUBLE   : return "SQRT";

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
	case OpCode::GET_VEC2_X    : return "GV2X";
	case OpCode::GET_VEC2_Y    : return "GV2Y";
	case OpCode::GET_VEC2_SWIZZLE: return "G2SW";
	case OpCode::ADD_VEC2      : return "ADD2";
	case OpCode::SUB_VEC2      : return "SUB2";
	case OpCode::DOT_VEC2      : return "DOT2";
	case OpCode::MUL_VEC2_DOUBLE: return "ML2F";
	case OpCode::MUL_DOUBLE_VEC2: return "MLF2";
	case OpCode::NEG_VEC2      : return "NEG2";
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

	case OpCode::POP_VOID  : return "POPV";
	case OpCode::POP_BOOL  : return "POPB";
	case OpCode::POP_INT   : return "POPI";
	case OpCode::POP_DOUBLE: return "POPD";
	case OpCode::POP_VEC2  : return "POP2";
	case OpCode::POP_VEC3  : return "POP3";
	case OpCode::POP_ARRAY : return "POPA";
	case OpCode::POP_STRUCT: return "POPS";	

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

void printStack(const std::vector<double>& stack, int numGlobals, int stackSize, double* ref)
{
	std::cout << "Stack: [";
	for (size_t i = 0; i < numGlobals; i++)
	{
		if (ref == &stack[i])
			std::cout << "*";

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
	ref.ptr = nullptr;
	const uint8_t* lastIP = m_program->code.data() + instructions;

	while (ip < lastIP)
	{
		OpCode instruction = (OpCode)readByte();

#ifndef NDEBUG
		if (m_debug)
		{
			std::cout << "IP: " << std::setw(4) << ip - 1;
			std::cout << " | Executing: " << IPToString((uint8_t)instruction) << " | ";
		}
#endif

		switch (instruction)
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
			pushDouble(m_stack[slot]);
			pushDouble(m_stack[slot + 1]);
			break;
		}

		case OpCode::GET_GLOBAL_VEC3:
		{
			uint8_t slot = readByte();
			pushDouble(m_stack[slot  ]);
			pushDouble(m_stack[slot+1]);
			pushDouble(m_stack[slot+2]);
			break;
		}

		case OpCode::GET_GLOBAL_ARRAY:
		{
			uint8_t size = readByte();
			uint8_t slot = readByte();
			copy(stackTop, slot, size);
			break;
		}

		case OpCode::GET_GLOBAL_STRUCT:
		{
			uint8_t size = readByte();
			uint8_t slot = readByte();
			copy(stackTop, slot, size);
			break;
		}

		case OpCode::GET_GLOBAL_REF:
		{
			assert(ref.ptr == nullptr); // make sure we don't overwrite an existing reference
			uint8_t slot = readByte();
			ref.ptr = &m_stack[slot];
			break;
		}

		case OpCode::GET_LOCAL_REF:
		{
			assert(ref.ptr == nullptr); // make sure we don't overwrite an existing reference
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			ref.ptr = &m_stack[frame.base + slot];
			break;
		}

		case OpCode::GET_MEMBER_REF:
		{
			uint8_t offset = readByte();
			double* slot = ref.ptr;
			ref.ptr = slot + offset;
			break;
		}

		case OpCode::GET_INDEX_REF:
		{
			uint8_t elemSize = readByte();
			int index = popInt();

			double* slot = ref.ptr;
			ref.ptr = slot + elemSize*index;
			break;
		}

		case OpCode::GET_INDEX_REF_BOOL:
		case OpCode::GET_INDEX_REF_INT:
		case OpCode::GET_INDEX_REF_DOUBLE:
		{
			int index = popInt();

			double* slot = ref.ptr;
			ref.ptr = slot + index;
			break;
		}

		case OpCode::GET_INDEX_REF_VEC2:
		{
			int index = popInt();

			double* slot = ref.ptr;
			ref.ptr = slot + index * 2;
			break;
		}

		case OpCode::GET_INDEX_REF_VEC3:
		{
			int index = popInt();

			double* slot = ref.ptr;
			ref.ptr = slot + index * 3;
			break;
		}

		case OpCode::GET_VEC2_X_REF:
		{
			double* v = ref.ptr;
			ref.ptr = v;
			break;
		}
		case OpCode::GET_VEC2_Y_REF:
		{
			double* v = ref.ptr;
			ref.ptr = v + 1;
			break;
		}

		case OpCode::GET_VEC3_X_REF:
		{
			double* v = ref.ptr;
			ref.ptr = v;
			break;
		}

		case OpCode::GET_VEC3_Y_REF:
		{
			double* v = ref.ptr;
			ref.ptr = v + 1;
			break;
		}

		case OpCode::GET_VEC3_Z_REF:
		{
			double* v = ref.ptr;
			ref.ptr = v + 2;
			break;
		}

		case OpCode::STORE_BOOL:
		{
			bool b = peekBool();
			*ref.ptr = b;
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_INT :
		{
			int n = peekInt();
			*ref.ptr = n;
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_DOUBLE:
		{
			*ref.ptr = peek();
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_VEC2:
		{
			vec2 v = peekVec2();
			double* xPtr = ref.ptr;
			xPtr[0] = v.x;
			xPtr[1] = v.y;
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_VEC3:
		{
			vec3 v = peekVec3();
			double* xPtr = ref.ptr;
			xPtr[0] = v.x;
			xPtr[1] = v.y;
			xPtr[2] = v.z;
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_ARRAY:
		{
			uint8_t size = readByte();
			memcpy(ref.ptr, m_stack.data() + (stackTop - size), size*sizeof(double));
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_STRUCT:
		{
			uint8_t size = readByte();
			memcpy(ref.ptr, m_stack.data() + (stackTop - size), size * sizeof(double));
			ref.ptr = nullptr;
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
			m_stack[slot] = peek();
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
			uint8_t size = readByte();
			uint8_t slot = readByte();
			copy(slot, (int)(stackTop - size), (int)size);
			break;
		}

		case OpCode::SET_GLOBAL_STRUCT:
		{
			uint8_t size = readByte();
			uint8_t slot = readByte();

			copy(slot, (int)(stackTop - size), (int)size);
			break;
		}

		case OpCode::GET_LOCAL_BOOL:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushBool((bool)m_stack[frame.base + slot]);
			break;
		}

		case OpCode::GET_LOCAL_INT:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushInt((int)m_stack[frame.base + slot]);
			break;
		}

		case OpCode::GET_LOCAL_DOUBLE:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushDouble(m_stack[frame.base + slot]);
			break;
		}

		case OpCode::GET_LOCAL_VEC2:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushVec2(getVec2At((int)frame.base + (int)slot));
			break;
		}

		case OpCode::GET_LOCAL_VEC3:
		{
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			pushVec3(getVec3At((int)frame.base + (int)slot));
			break;
		}

		case OpCode::GET_LOCAL_ARRAY:
		{
			uint8_t size = readByte();
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			copy(stackTop, frame.base + slot, size);
			break;
		}

		case OpCode::GET_LOCAL_STRUCT:
		{
			uint8_t size = readByte();
			uint8_t slot = readByte();
			CallFrame& frame = currentFrame();
			copy(stackTop, frame.base + slot, size);
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
		case OpCode::SUB_INT:
		case OpCode::MUL_INT:
		case OpCode::DIV_INT:
		case OpCode::EXP_INT:
		case OpCode::GT_INT:
		case OpCode::LT_INT:
		case OpCode::GE_INT:
		case OpCode::LE_INT:
		{
			int b = popInt();
			int a = popInt();

			switch (instruction)
			{
			case OpCode::ADD_INT: pushInt(a + b); break;
			case OpCode::SUB_INT: pushInt(a - b); break;
			case OpCode::MUL_INT: pushInt(a * b); break;
			case OpCode::DIV_INT:
				if (b == 0)
					throw std::runtime_error("division by zero.");
				pushInt(a / b);
				break;
			case OpCode::EXP_INT:
				if (b < 0)
					throw std::runtime_error("Negative exponent not supported for integers.");

				pushInt(ipow(a, b));
				break;
			case OpCode::GT_INT: pushBool(a > b); break;
			case OpCode::LT_INT: pushBool(a < b); break;
			case OpCode::GE_INT: pushBool(a >= b); break;
			case OpCode::LE_INT: pushBool(a <= b); break;
			}

			break;
		}

		// Double unary operators
		case OpCode::NEG_DOUBLE:
		case OpCode::SQR_DOUBLE:
		case OpCode::SQRT_DOUBLE:
		{
			double a = popDouble();

			switch (instruction)
			{
			case OpCode::NEG_DOUBLE: pushDouble(-a); break;
			case OpCode::SQR_DOUBLE: pushDouble(a * a); break;
			case OpCode::SQRT_DOUBLE: pushDouble(std::sqrt(a)); break;
			}
			break;
		}

		// double binary operators
		case OpCode::ADD_DOUBLE:
		case OpCode::SUB_DOUBLE:
		case OpCode::MUL_DOUBLE:
		case OpCode::DIV_DOUBLE:
		case OpCode::EXP_DOUBLE:
		case OpCode::GT_DOUBLE:
		case OpCode::LT_DOUBLE:
		case OpCode::GE_DOUBLE:
		case OpCode::LE_DOUBLE:
		{
			double b = popDouble();
			double a = popDouble();

			switch (instruction)
			{
			case OpCode::ADD_DOUBLE: pushDouble(a + b); break;
			case OpCode::SUB_DOUBLE: pushDouble(a - b); break;
			case OpCode::MUL_DOUBLE: pushDouble(a * b); break;
			case OpCode::EXP_DOUBLE: pushDouble(std::pow(a, b)); break;
			case OpCode::DIV_DOUBLE: 
				if (b == 0.0)
					throw std::runtime_error("division by zero.");
				pushDouble(a / b); 
				break;
			case OpCode::GT_DOUBLE: pushBool(a > b); break;
			case OpCode::LT_DOUBLE: pushBool(a < b); break;
			case OpCode::GE_DOUBLE: pushBool(a >= b); break;
			case OpCode::LE_DOUBLE: pushBool(a <= b); break;
			}
			break;
		}

		case OpCode::ADD_VEC2:
		case OpCode::SUB_VEC2:
		case OpCode::DOT_VEC2:
		{
			popVec2_1();
			popVec2_0();

			switch (instruction)
			{
			case OpCode::ADD_VEC2: pushVec2(vec2_0 + vec2_1); break;
			case OpCode::SUB_VEC2: pushVec2(vec2_0 - vec2_1); break;
			case OpCode::DOT_VEC2: pushDouble(vec2_0 * vec2_1); break;
			}
			break;
		}
		case OpCode::MUL_VEC2_DOUBLE:
		{
			double scalar = popDouble();
			popVec2_0();
			pushVec2(vec2_0 * scalar);
			break;
		}
		case OpCode::MUL_DOUBLE_VEC2:
		{
			popVec2_0();
			double scalar = popDouble();
			pushVec2(vec2_0 * scalar);
			break;
		}

		case OpCode::GET_VEC2_X:
		case OpCode::GET_VEC2_Y:
		{
			popVec2_0();
			switch (instruction)
			{
			case OpCode::GET_VEC2_X: pushDouble(vec2_0.x); break;
			case OpCode::GET_VEC2_Y: pushDouble(vec2_0.y); break;
			}
			break;
		}
		case OpCode::NEG_VEC2:
		{
			popVec2_0();
			pushVec2(vec2(-vec2_0.x, -vec2_0.y));
			break;
		}

		case OpCode::GET_VEC2_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			popVec2_0();
			double c[4] = { vec2_0.x, vec2_0.y, 0.0, 0.0 };

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
		case OpCode::ADD_VEC3:
		case OpCode::SUB_VEC3:
		case OpCode::DOT_VEC3:
		{
			popVec3_1();
			popVec3_0();

			switch (instruction)
			{
			case OpCode::ADD_VEC3: pushVec3(vec3_0 + vec3_1); break;
			case OpCode::SUB_VEC3: pushVec3(vec3_0 - vec3_1); break;
			case OpCode::DOT_VEC3: pushDouble(vec3_0 * vec3_1); break;
			}
			break;
		}
		case OpCode::MUL_VEC3_DOUBLE:
		{
			double scalar = popDouble();
			popVec3_0();
			pushVec3(vec3_0 * scalar);
			break;
		}
		case OpCode::MUL_DOUBLE_VEC3:
		{
			popVec3_0();
			double scalar = popDouble();
			pushVec3(vec3_0 * scalar);
			break;
		}
		case OpCode::GET_VEC3_X:
		{
			popVec3_0();
			pushDouble(vec3_0.x);
			break;
		}
		case OpCode::GET_VEC3_Y:
		{
			popVec3_0();
			pushDouble(vec3_0.y);
			break;
		}
		case OpCode::GET_VEC3_Z:
		{
			popVec3_0();
			pushDouble(vec3_0.z);
			break;
		}
		case OpCode::GET_VEC3_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			popVec3_0();
			double c[4] = { vec3_0.x, vec3_0.y, vec3_0.z, 0.0 };

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
			popVec3_0();
			pushVec3(vec3(-vec3_0.x, -vec3_0.y, -vec3_0.z));
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
			// Nothing to do here since structs are created on the stack.
			break;
		}
		case OpCode::COPY_STRUCT:
		{
			// Nothing to do here. 
			break;
		}
		case OpCode::GET_PROPERTY_BOOL:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushBool(obj.fields[slot].b);
			break;
		}

		case OpCode::GET_PROPERTY_INT:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushInt(obj.fields[slot].i);
			break;
		}

		case OpCode::GET_PROPERTY_DOUBLE:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushDouble(obj.fields[slot].d);
			break;
		}

		case OpCode::GET_PROPERTY_VEC2:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushVec2(obj.fields[slot].vec2Value);
			break;
		}

		case OpCode::GET_PROPERTY_VEC3:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushVec3(obj.fields[slot].vec3Value);
			break;
		}

		case OpCode::GET_PROPERTY_ARRAY:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushArray(*obj.fields[slot].arrayValue);
			break;
		}

		case OpCode::GET_PROPERTY_STRUCT:
		{
			uint8_t typeIndex = readByte();
			StructValue obj = popStruct(m_program->types.getStructType(typeIndex));
			uint8_t slot = readByte();
			pushStruct(*obj.fields[slot].structValue);
			break;
		}

		case OpCode::GET_INDEX_BOOL:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

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
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

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
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

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
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

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
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushVec3(arr.elements[index].vec3Value);
			break;
		}

		case OpCode::GET_INDEX_ARRAY:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			Value& element = arr.elements[index];
			pushArray(*element.arrayValue); break;
			break;
		}

		case OpCode::GET_INDEX_STRUCT:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();
			ArrayValue arr = popArray(m_program->types.getArrayType((int)typeIndex));

#ifndef  NDEBUG
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			Value& element = arr.elements[index];
			pushStruct(*element.structValue); break;
			break;
		}

		case OpCode::JUMP:
		{
			uint16_t offset = readUint16();
			ip += offset;
			break;
		}

		case OpCode::JUMP_IF_FALSE:
		{
			uint16_t offset = readUint16();
			if (peek() == 0.0)
				ip += offset;
			break;
		}

		case OpCode::JUMP_IF_TRUE:
		{
			uint16_t offset = readUint16();
			if (peek() != 0.0)
				ip += offset;
			break;
		}

		case OpCode::LOOP:
		{
			uint16_t offset = readUint16();
			ip -= offset;
			break;
		}

		case OpCode::POP_VOID:
		{
			popVoid();
			break;
		}

		case OpCode::POP_BOOL:
		{
			popBool();
			break;
		}

		case OpCode::POP_INT:
		{
			popInt();
			break;
		}

		case OpCode::POP_DOUBLE:
		{
			popDouble();
			break;
		}

		case OpCode::POP_VEC2:
		{
			popVec2();
			break;
		}

		case OpCode::POP_VEC3:
		{
			popVec3();
			break;
		}

		case OpCode::POP_ARRAY:
		{
			uint8_t typeIndex = readByte();
			Type type = m_program->types.getArrayType(typeIndex);
			popValues(type->size());
			break;
		}

		case OpCode::POP_STRUCT:
		{
			uint8_t typeIndex = readByte();
			Type type = m_program->types.getStructType(typeIndex);
			popValues(type->size());
			break;
		}

		case OpCode::CALL:
		{
			uint8_t fnIndex = readByte();
			uint8_t args = readByte();
			callFunction(fnIndex, args);
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
			TypeKind returnType = (TypeKind)((uint8_t)instruction - (uint8_t)OpCode::RETURN_VOID);

			Value result;
			if (stackTop > globalStackSize)
			{
				switch (returnType)
				{
				case TypeKind::Void  : pop(); break;
				case TypeKind::Bool  : result = popBool  (); break;
				case TypeKind::Int   : result = popInt   (); break;
				case TypeKind::Double: result = popDouble(); break;
				case TypeKind::Vec2  : result = popVec2  (); break;
				case TypeKind::Vec3  : result = popVec3  (); break;
				case TypeKind::Array:
				{
					uint8_t typeIndex = readByte();
					result = std::make_shared<ArrayValue>(popArray(m_program->types.getArrayType(typeIndex)));
					break;
				}
				case TypeKind::Struct:
				{
					uint8_t typeIndex = readByte();
					result = std::make_shared<StructValue>(popStruct(m_program->types.getStructType(typeIndex)));
					break;
				}
				default:
					result = pop();
				}
			}

			CallFrame& frame = m_frames[--frameCount];
			stackTop = frame.base;
			ip = frame.ip;

			if (frameCount == 0)
				return result;

			switch (returnType)
			{
			case TypeKind::Void  : pushVoid  (); break;
			case TypeKind::Bool  : pushBool  (result.b); break;
			case TypeKind::Int   : pushInt   (result.i); break;
			case TypeKind::Double: pushDouble(result.d); break;
			case TypeKind::Vec2  : pushVec2  (result.vec2Value); break;
			case TypeKind::Vec3  : pushVec3  (result.vec3Value); break;
			case TypeKind::Array : pushArray (*result.arrayValue); break;
			case TypeKind::Struct: pushStruct(*result.structValue); break;
			default:
				throw std::runtime_error("Unsupported return type");
			}

#ifndef NDEBUG
			if (m_debug)
			{
				printStack(m_stack, (int)globalStackSize, (int)stackTop, ref.ptr);
			}
#endif

			break;
		}

		default:
			throw std::runtime_error("Unknown opcode");
		}

#ifndef NDEBUG
		if (m_debug && (instruction < OpCode::RETURN_VOID))
		{
			printStack(m_stack, (int)globalStackSize, (int)stackTop, ref.ptr);
		}
#endif
	}

	throw std::runtime_error("Unexpected end of code.");
}

void VM::callFunction(int fnIndex, int argCount)
{
	const FunctionInfo& fn = m_program->functions[fnIndex];

	if (argCount != fn.args.size())
		throw std::runtime_error("Arity mismatch in call to " + fn.name);

	if (fn.isNative)
	{
		FuncArgs args;
		args.count = argCount;
		args.stack = &m_stack[stackTop - fn.argSize];

		Value result = fn.fnc(args);
		stackTop -= fn.argSize;

		switch (fn.returnType->kind)
		{
		case TypeKind::Void  : pushVoid  (); break;
		case TypeKind::Bool  : pushBool  (result.b); break;
		case TypeKind::Int   : pushInt   (result.i); break;
		case TypeKind::Double: pushDouble(result.d); break;
		case TypeKind::Vec2  : pushVec2  (result.vec2Value); break;
		case TypeKind::Vec3  : pushVec3  (result.vec3Value); break;
		default:
			throw std::runtime_error("Unsupported return type from native function: " + std::to_string((int)fn.returnType->kind));
		}
		return;
	}

	if (frameCount >= MAX_CALL_DEPTH)
		throw std::runtime_error("Stack overflow: too many nested function calls.");

	CallFrame& frame = m_frames[frameCount++];
	frame.functionIndex = fnIndex;
	frame.base = stackTop - fn.argSize;
	frame.ip = ip;
	ip = &(m_program->code[fn.entry]);
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
