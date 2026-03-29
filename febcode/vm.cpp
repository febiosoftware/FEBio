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
			std::cout << "IP: " << std::setw(4) << (ip - &m_program->code[0]) - 1;
			std::cout << " | Executing: " << OpCodeToString(instruction) << " | ";
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

		case OpCode::GET_GLOBAL_MAT2:
		{
			uint8_t slot = readByte();
			pushDouble(m_stack[slot    ]);
			pushDouble(m_stack[slot + 1]);
			pushDouble(m_stack[slot + 2]);
			pushDouble(m_stack[slot + 3]);
			break;
		}

		case OpCode::GET_GLOBAL_MAT3:
		{
			uint8_t slot = readByte();
			pushFrom(slot, 9);
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

		case OpCode::STORE_MAT2:
		{
			mat2 v = peekMat2();
			double* xPtr = ref.ptr;
			xPtr[0] = v.m[0][0];
			xPtr[1] = v.m[0][1];
			xPtr[2] = v.m[1][0];
			xPtr[3] = v.m[1][1];
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_MAT3:
		{
			double* xPtr = ref.ptr;
			memcpy(xPtr, &m_stack[stackTop - 9], 9 * sizeof(double)); // copy all 9 elements at once
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_ARRAY:
		{
			uint8_t size = readByte();
			memcpy(ref.ptr, &m_stack[stackTop - size], size * sizeof(double));
			ref.ptr = nullptr;
			break;
		}

		case OpCode::STORE_STRUCT:
		{
			uint8_t size = readByte();
			memcpy(ref.ptr, &m_stack[stackTop - size], size * sizeof(double));
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

		case OpCode::SET_GLOBAL_MAT2:
		{
			uint8_t slot = readByte();
			setMat2At(slot);
			break;
		}

		case OpCode::SET_GLOBAL_MAT3:
		{
			uint8_t slot = readByte();
			setMat3At(slot);
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
			vec2& vec2_1 = popVec2();
			vec2& vec2_0 = popVec2();

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
			double* v = peekPtr(2);
			v[0] *= scalar;
			v[1] *= scalar;
			break;
		}
		case OpCode::MUL_DOUBLE_VEC2:
		{
			double* v = popPtr(2);
			double scalar = popDouble();
			pushDouble(v[0] * scalar);
			pushDouble(v[1] * scalar);
			break;
		}

		case OpCode::GET_VEC2_X:
		{
			double* v = popPtr(2);
			pushDouble(v[0]);
			break;
		}
		case OpCode::GET_VEC2_Y:
		{
			double* v = popPtr(2);
			pushDouble(v[1]);
			break;
		}
		case OpCode::NEG_VEC2:
		{
			double* v = peekPtr(2);
			v[0] = -v[0];
			v[1] = -v[1];
			break;
		}

		case OpCode::GET_VEC2_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			double* v = popPtr(2);
			double c[4] = { v[0], v[1], 0.0, 0.0};

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

		case OpCode::GET_VEC2_INDEX:
		{
			int index = popInt();
			double* v = popPtr(2);
#ifndef NDEBUG
			if (index < 0 || index > 1)
				throw std::runtime_error("vec2 index out of bounds.");
#endif
			pushDouble(v[index]);
			break;
		}

		case OpCode::ADD_VEC3:
		{
			double* b = popPtr(3);
			double* a = peekPtr(3);
			a[0] += b[0];
			a[1] += b[1];
			a[2] += b[2];
			break;
		}
		case OpCode::SUB_VEC3:
		{
			double* b = popPtr(3);
			double* a = peekPtr(3);
			a[0] -= b[0];
			a[1] -= b[1];
			a[2] -= b[2];
			break;
		}
		case OpCode::DOT_VEC3:
		{
			double* b = popPtr(3);
			double* a = popPtr(3);
			pushDouble(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
			break;
		}
		case OpCode::MUL_VEC3_DOUBLE:
		{
			double scalar = popDouble();
			double* a = peekPtr(3);
			a[0] *= scalar;
			a[1] *= scalar;
			a[2] *= scalar;
			break;
		}
		case OpCode::MUL_DOUBLE_VEC3:
		{
			double* a = popPtr(3);
			double scalar = popDouble();
			pushDouble(a[0] * scalar);
			pushDouble(a[1] * scalar);
			pushDouble(a[2] * scalar);
			break;
		}
		case OpCode::GET_VEC3_X:
		{
			double* a = popPtr(3);
			pushDouble(a[0]);
			break;
		}
		case OpCode::GET_VEC3_Y:
		{
			double* a = popPtr(3);
			pushDouble(a[1]);
			break;
		}
		case OpCode::GET_VEC3_Z:
		{
			double* a = popPtr(3);
			pushDouble(a[2]);
			break;
		}
		case OpCode::GET_VEC3_SWIZZLE:
		{
			uint8_t mask = readByte();
			uint8_t size = readByte();
			double* a = popPtr(3);
			double c[4] = { a[0], a[1], a[2], 0.0 };

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
			double* a = popPtr(3);
			a[0] = -a[0];
			a[1] = -a[1];
			a[2] = -a[2];
			break;
		}

		case OpCode::GET_VEC3_INDEX:
		{
			int index = popInt();
			double* a = popPtr(3);
#ifndef NDEBUG
			if (index < 0 || index > 2)
				throw std::runtime_error("vec3 index out of bounds.");
#endif
			pushDouble(a[index]);
			break;
		}

		// ----- Mat2 operators ------
		case OpCode::ADD_MAT2:
		{
			mat2& b = popMat2();
			mat2& a = popMat2();
			pushMat2(a + b);
			break;
		}
		case OpCode::SUB_MAT2:
		{
			mat2& b = popMat2();
			mat2& a = popMat2();
			pushMat2(a - b);
			break;
		}
		case OpCode::MUL_MAT2:
		{
			mat2& b = popMat2();
			mat2& a = popMat2();
			pushMat2(a * b);
			break;
		}

		case OpCode::MUL_DOUBLE_MAT2:
		{
			mat2& a = popMat2();
			double scalar = popDouble();
			pushMat2(a * scalar);
			break;
		}
		case OpCode::MUL_MAT2_DOUBLE:
		{
			double scalar = popDouble();
			mat2& a = popMat2();
			pushMat2(a * scalar);
			break;
		}

		case OpCode::MUL_MAT2_VEC2:
		{
			vec2& v = popVec2();
			mat2& A = popMat2();
			pushVec2(A * v);
			break;
		}

		case OpCode::GET_MAT2_INDEX:
		{
			int index = popInt();
			mat2& A = popMat2();
#ifndef NDEBUG
			if (index < 0 || index > 1)
				throw std::runtime_error("mat2 index out of bounds.");
#endif
			pushVec2(*((vec2*)(&A.m[index][0])));
			break;
		}

		// ----- Mat3 operators ------
		case OpCode::ADD_MAT3:
		{
			mat3& B = popMat3();
			mat3& A = popMat3();
			pushMat3(A + B);
			break;
		}
		case OpCode::SUB_MAT3:
		{
			mat3& B = popMat3();
			mat3& A = popMat3();
			pushMat3(A - B);
			break;
		}
		case OpCode::MUL_MAT3:
		{
			mat3& B = popMat3();
			mat3& A = popMat3();
			pushMat3(A * B);
			break;
		}
		case OpCode::MUL_DOUBLE_MAT3:
		{
			mat3& A = popMat3();
			double scalar = popDouble();
			pushMat3(A * scalar);
			break;
		}
		case OpCode::MUL_MAT3_DOUBLE:
		{
			double scalar = popDouble();
			mat3& A = popMat3();
			pushMat3(A * scalar);
			break;
		}

		case OpCode::MUL_MAT3_VEC3:
		{
			vec3& v = popVec3();
			mat3& A = popMat3();
			pushVec3(A * v);
			break;
		}

		case OpCode::GET_MAT3_INDEX:
		{
			int index = popInt();
			mat3& A = popMat3();
#ifndef NDEBUG
			if (index < 0 || index > 2)
				throw std::runtime_error("mat3 index out of bounds.");
#endif
			pushVec3(*((vec3*)(&A.m[index][0])));
			break;
		}

		case OpCode::ADD_GLOBAL_MAT3:
		{
			uint8_t slotA = readByte();
			uint8_t slotB = readByte();
			mat3& A = getMat3At(slotA);
			mat3& B = getMat3At(slotB);
			pushMat3(A + B);
			break;
		}

		case OpCode::SUB_GLOBAL_MAT3:
		{
			uint8_t slotA = readByte();
			uint8_t slotB = readByte();
			mat3& A = getMat3At(slotA);
			mat3& B = getMat3At(slotB);
			pushMat3(A - B);
			break;
		}

		case OpCode::MUL_GLOBAL_MAT3:
		{
			uint8_t slotA = readByte();
			uint8_t slotB = readByte();
			mat3& A = getMat3At(slotA);
			mat3& B = getMat3At(slotB);
			pushMat3(A * B);
			break;
		}

		// ----- Logical operators -----
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

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Bool);

			double* arr = popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushBool(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_INT:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Int);

			double* arr = popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushInt(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_DOUBLE:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Double);

			double* arr = popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushDouble(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_VEC2:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Vec2);

			vec2* arr = (vec2*) popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushVec2(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_VEC3:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Vec3);

			vec3* arr = (vec3*)popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushVec3(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_MAT2:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Mat2);

			mat2* arr = (mat2*)popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushMat2(arr[index]);
			break;
		}

		case OpCode::GET_INDEX_MAT3:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Mat3);

			mat3* arr = (mat3*)popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->size())
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			pushMat3(arr[index]);
			break;
		}
		case OpCode::GET_INDEX_ARRAY:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Array);

			double* arr = popPtr(type->size());

			size_t elemSize = type->elementType->size();

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			push(arr + index * elemSize, elemSize); break;
			break;
		}

		case OpCode::GET_INDEX_STRUCT:
		{
			int index = popInt();
			uint8_t typeIndex = readByte();

			Type type = m_program->types.getArrayType((int)typeIndex);
			assert(type->kind == TypeKind::Array);
			assert(type->elementType->kind == TypeKind::Struct);

			double* arr = popPtr(type->size());

#ifndef  NDEBUG
			if (index >= type->arraySize)
				throw std::runtime_error("Array index out of bounds.");
#endif // ! NDEBUG

			size_t structSize = type->elementType->size();

			push(arr + index * structSize, structSize);
			break;
		}

		case OpCode::GET_GLOBAL_INDEX_DOUBLE:
		{
			int index = popInt();
			uint8_t slot = readByte();
			double* arr = getPtrAt(slot);
			pushDouble(arr[index]);
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
			pop(2);
			break;
		}

		case OpCode::POP_VEC3:
		{
			pop(3);
			break;
		}

		case OpCode::POP_MAT2:
		{
			pop(4);
			break;
		}

		case OpCode::POP_MAT3:
		{
			pop(9);
			break;
		}

		case OpCode::POP_ARRAY:
		{
			uint8_t size = readByte();
			pop(size);
			break;
		}

		case OpCode::POP_STRUCT:
		{
			uint8_t typeIndex = readByte();
			Type type = m_program->types.getStructType(typeIndex);
			pop(type->size());
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
		case OpCode::RETURN_MAT2:
		case OpCode::RETURN_MAT3:
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
				case TypeKind::Mat2  : result = popMat2  (); break;
				case TypeKind::Mat3  : result = popMat3  (); break;
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
