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
	case OpCode::PUSH_CONST    : return "MOV ";
	case OpCode::GET_GLOBAL    : return "GETG";
	case OpCode::SET_GLOBAL    : return "SETG";
	case OpCode::STORE         : return "STRE";
	case OpCode::GET_GLOBAL_REF: return "GREF";
	case OpCode::GET_LOCAL_REF : return "LREF";
	case OpCode::GET_INDEX_REF : return "IREF";
	case OpCode::GET_MEMBER_REF: return "MREF";
	case OpCode::GET_LOCAL     : return "GETL";
	case OpCode::SET_LOCAL     : return "SETL";
	case OpCode::ADD_INT       : return "ADDI";
	case OpCode::ADD_DOUBLE    : return "ADDF";
	case OpCode::ADD_STRING    : return "ADDS";
	case OpCode::SUB_INT       : return "SUBI";
	case OpCode::SUB_DOUBLE    : return "SUBF";
	case OpCode::MUL_INT       : return "MULI";
	case OpCode::MUL_DOUBLE    : return "MULF";
	case OpCode::DIV_INT       : return "DIVI";
	case OpCode::DIV_DOUBLE    : return "DIVF";
	case OpCode::EXP_INT       : return "EXPI";
	case OpCode::EXP_DOUBLE    : return "EXPF";
	case OpCode::EQUAL         : return "EQ  ";
	case OpCode::NOT_EQUAL     : return "NEQ ";
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
	case OpCode::COPY_VEC2     : return "CPV2";
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
	case OpCode::COPY_VEC3     : return "CPV3";
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
	case OpCode::GET_PROPERTY  : return "GETP";
	case OpCode::CREATE_ARRAY  : return "ARR ";
	case OpCode::COPY_ARRAY    : return "CPYA";
	case OpCode::GET_INDEX     : return "GETI";
	case OpCode::JUMP          : return "JMP ";
	case OpCode::JUMP_IF_FALSE : return "JMPF";
	case OpCode::JUMP_IF_TRUE  : return "JMPT";
	case OpCode::LOOP          : return "LOOP";
	case OpCode::CALL          : return "CALL";
	case OpCode::CALL_BINARY   : return "BINO";
	case OpCode::PRINT         : return "PRNT";
	case OpCode::RETURN        : return "RET ";
	case OpCode::POP           : return "POP ";
	case OpCode::GET_VEC2_X_REF: return "RV2X";
	case OpCode::GET_VEC2_Y_REF: return "RV2Y";
	case OpCode::GET_VEC3_X_REF: return "RV3X";
	case OpCode::GET_VEC3_Y_REF: return "RV3Y";
	case OpCode::GET_VEC3_Z_REF: return "RV3Z";
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
	if (m_frames.empty())
		throw std::runtime_error("No function to execute");

	const size_t instructions = m_program.code.size();

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
		case OpCode::PUSH_CONST:
		{
			uint16_t idx = readUint16();
			push(m_program.constants[idx]);
			break;
		}

		case OpCode::GET_GLOBAL:
		{
			uint16_t slot = readUint16();
			push(m_stack[slot]);
			break;
		}

		case OpCode::GET_GLOBAL_REF:
		{
			uint16_t slot = readUint16();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &m_stack[slot];
			ref.ref.type = TypeKind::Value;

			push(ref);
			break;
		}

		case OpCode::GET_LOCAL_REF:
		{
			uint16_t slot = readUint16();

			CallFrame& frame = currentFrame();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &m_stack[frame.base + slot];
			ref.ref.type = TypeKind::Value;

			push(ref);
			break;
		}

		case OpCode::GET_MEMBER_REF:
		{
			uint16_t memberIndex = readUint16();

			const Value& objRef = pop();
			Value* slot = (Value*)objRef.ref.ptr;

			StructValue* s = slot->structValue.get();

			Value ref;
			ref.index = ValueIndex::REF;
			ref.ref.ptr = &s->fields[memberIndex];
			ref.ref.type = TypeKind::Value;

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
			ref.ref.type = TypeKind::Value;

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
			ref.ref.type = TypeKind::Double;

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
			ref.ref.type = TypeKind::Double;

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
			ref.ref.type = TypeKind::Double;

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
			ref.ref.type = TypeKind::Double;

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
			ref.ref.type = TypeKind::Double;

			push(ref);
			break;
		}

		case OpCode::STORE:
		{
			const Value& value = pop();
			const Value& ref = pop();

			const Ref& r = ref.ref;

			switch (r.type)
			{
			case TypeKind::Value:
				*(Value*)r.ptr = value;
				break;

			case TypeKind::Double:
				*(double*)r.ptr = value.d;
				break;

			default:
				throw std::runtime_error("Invalid store target");
			}

			push(value);
			break;
		}

		case OpCode::SET_GLOBAL:
		{
			uint16_t slot = readUint16();
			m_stack[slot] = peek();
			break;
		}

		case OpCode::GET_LOCAL:
		{
			uint16_t slot = readUint16();
			CallFrame& frame = currentFrame();
			push(m_stack[frame.base + slot]);
			break;
		}

		case OpCode::SET_LOCAL:
		{
			uint16_t slot = readUint16();
			CallFrame& frame = currentFrame();
			m_stack[frame.base + slot] = peek();
			break;
		}

		// Integer operators
		case OpCode::NEG_INT:
		{
			const Value& a = pop();
			push(-getInt(a));
			break;
		}
		case OpCode::ADD_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) + getInt(b));
			break;
		}
		case OpCode::SUB_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) - getInt(b));
			break;
		}
		case OpCode::MUL_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) * getInt(b));
			break;
		}
		case OpCode::DIV_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			int divisor = getInt(b);
			if (divisor == 0)
				throw std::runtime_error("division by zero.");
			push(getInt(a) / divisor);
			break;
		}
		case OpCode::EXP_INT:
		{
			const Value& b = pop();
			const Value& a = pop();

			int e = getInt(b);
			if (e < 0)
				throw std::runtime_error("Negative exponent not supported for integers.");

			push(ipow(getInt(a), e));
			break;
		}
		case OpCode::GT_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) > getInt(b));
			break;
		}
		case OpCode::LT_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) < getInt(b));
			break;
		}
		case OpCode::GE_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) >= getInt(b));
			break;
		}
		case OpCode::LE_INT:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getInt(a) <= getInt(b));
			break;
		}

		// Double operators
		case OpCode::NEG_DOUBLE:
		{
			const Value& a = pop();
			push(-getDouble(a));
			break;
		}

		case OpCode::ADD_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) + getDouble(b));
			break;
		}
		case OpCode::SUB_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) - getDouble(b));
			break;
		}
		case OpCode::MUL_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) * getDouble(b));
			break;
		}
		case OpCode::DIV_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			double divisor = getDouble(b);
			if (divisor == 0.0)
				throw std::runtime_error("division by zero.");
			push(getDouble(a) / divisor);
			break;
		}
		case OpCode::EXP_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(std::pow(getDouble(a), getDouble(b)));
			break;
		}
		case OpCode::GT_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) > getDouble(b));
			break;
		}
		case OpCode::LT_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) < getDouble(b));
			break;
		}
		case OpCode::GE_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) >= getDouble(b));
			break;
		}
		case OpCode::LE_DOUBLE:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getDouble(a) <= getDouble(b));
			break;
		}

		case OpCode::CREATE_VEC2:
		{
			const Value& y = pop();
			const Value& x = pop();
			push(vec2(getDouble(x), getDouble(y)));
			break;
		}
		case OpCode::COPY_VEC2:
		{
			const Value& src = pop();
			const vec2& v = getVec2(src);
			push(vec2(v.x, v.y));
			break;
		}
		case OpCode::ADD_VEC2:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec2(a) + getVec2(b));
			break;
		}
		case OpCode::SUB_VEC2:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec2(a) - getVec2(b));
			break;
		}
		case OpCode::DOT_VEC2:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec2(a) * getVec2(b));
			break;
		}
		case OpCode::MUL_VEC2_DOUBLE:
		{
			const Value& scalar = pop();
			const Value& vec = pop();
			push(getVec2(vec) * getDouble(scalar));
			break;
		}
		case OpCode::MUL_DOUBLE_VEC2:
		{
			const Value& vec = pop();
			const Value& scalar = pop();
			push(getVec2(vec) * getDouble(scalar));
			break;
		}
		case OpCode::GET_VEC2_X:
		{
			const Value& vec = pop();
			push(getVec2(vec).x);
			break;
		}
		case OpCode::GET_VEC2_Y:
		{
			const Value& vec = pop();
			push(getVec2(vec).y);
			break;
		}
		case OpCode::NEG_VEC2:
		{
			const Value& vec = pop();
			const vec2& v = getVec2(vec);
			push(vec2(-v.x, -v.y));
			break;
		}
		case OpCode::GET_VEC2_SWIZZLE:
		{
			uint16_t mask = readUint16();
			uint16_t size = readUint16();
			const Value& vec = pop();
			const vec2& v = getVec2(vec);
			double c[4] = { v.x, v.y, 0.0, 0.0 };

			if (size == 2)
			{
				vec2 result;
				result.y = c[mask & 0b11];
				result.x = c[(mask >> 2) & 0b11];
				push(result);
			}
			else if (size == 3)
			{
				vec3 result;
				result.z = c[mask & 0b11];
				result.y = c[(mask >> 2) & 0b11];
				result.x = c[(mask >> 4) & 0b11];
				push(result);
			}
			break;
		}
		case OpCode::CREATE_VEC3:
		{
			const Value& z = pop();
			const Value& y = pop();
			const Value& x = pop();
			push(vec3(getDouble(x), getDouble(y), getDouble(z)));
			break;
		}
		case OpCode::COPY_VEC3:
		{
			const Value& src = pop();
			const vec3& v = getVec3(src);
			push(vec3(v.x, v.y, v.z));
			break;
		}
		case OpCode::ADD_VEC3:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec3(a) + getVec3(b));
			break;
		}
		case OpCode::SUB_VEC3:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec3(a) - getVec3(b));
			break;
		}
		case OpCode::DOT_VEC3:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(getVec3(a) * getVec3(b));
			break;
		}
		case OpCode::MUL_VEC3_DOUBLE:
		{
			const Value& scalar = pop();
			const Value& vec = pop();
			push(getVec3(vec) * getDouble(scalar));
			break;
		}
		case OpCode::MUL_DOUBLE_VEC3:
		{
			const Value& vec = pop();
			const Value& scalar = pop();
			push(getVec3(vec) * getDouble(scalar));
			break;
		}
		case OpCode::GET_VEC3_X:
		{
			const Value& vec = pop();
			push(getVec3(vec).x);
			break;
		}
		case OpCode::GET_VEC3_Y:
		{
			const Value& vec = pop();
			push(getVec3(vec).y);
			break;
		}
		case OpCode::GET_VEC3_Z:
		{
			const Value& vec = pop();
			push(getVec3(vec).z);
			break;
		}
		case OpCode::GET_VEC3_SWIZZLE:
		{
			uint16_t mask = readUint16();
			uint16_t size = readUint16();
			const Value& vec = pop();
			const vec3& v = getVec3(vec);
			double c[4] = { v.x, v.y, v.z, 0.0 };

			if (size == 2)
			{
				vec2 result;
				result.y = c[mask & 0b11];
				result.x = c[(mask >> 2) & 0b11];
				push(result);
			}
			else if (size == 3)
			{
				vec3 result;
				result.z = c[mask & 0b11];
				result.y = c[(mask >> 2) & 0b11];
				result.x = c[(mask >> 4) & 0b11];
				push(result);
			}
			break;
		}

		case OpCode::NEG_VEC3:
		{
			const Value& vec = pop();
			const vec3& v = getVec3(vec);
			push(vec2(-v.x, -v.y));
			break;
		}

		// string operators
		case OpCode::ADD_STRING:
		{
			const Value& a = pop();
			const Value& b = pop();
			push(getString(a) + getString(b));
			break;
		}

		// Logical operators
		case OpCode::NOT:
		{
			const Value& a = pop();
			bool b = isTruthy(a);
			push(!b);
			break;
		}
		case OpCode::EQUAL:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(a == b);
			break;
		}
		case OpCode::NOT_EQUAL:
		{
			const Value& b = pop();
			const Value& a = pop();
			push(a != b);
			break;
		}

		case OpCode::CREATE_STRUCT:
		{
			uint16_t typeIndex = readUint16();
			Type type = m_program.types.getStructType(typeIndex);
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
		case OpCode::GET_PROPERTY:
		{
			uint16_t slot = readUint16();

			const Value& objVal = pop();

			switch (ValueType(objVal))
			{
			case TypeKind::Struct:
			{
				const auto& structVal = getStruct(objVal);

				if (slot >= structVal.fields.size())
					throw std::runtime_error("Invalid property slot");

				push(structVal.fields[slot]);
				break;
			}
			default:
				throw std::runtime_error("GET_PROPERTY on non-struct");
			}

			break;
		}

		case OpCode::CREATE_ARRAY: {
			uint16_t count = readUint16();
			auto arr = std::make_shared<ArrayValue>();

			arr->elements.resize(count);

			for (int i = count - 1; i >= 0; --i) {
				arr->elements[i] = pop();
			}

			Type elemType = m_program.types.getBuiltinType(arr->elements[0]);
			arr->type = m_program.types.getArrayType(elemType, count);

			push(arr);
			break;
		}
		case OpCode::COPY_ARRAY: {
			auto src = pop();

			if (!isArray(src))
				throw std::runtime_error("COPY_ARRAY operand must be an array");

			const ArrayValue& srcObj = getArray(src);
			auto obj = copyArray(srcObj);

			push(obj);
			break;
		}
		case OpCode::GET_INDEX: {
			const Value& indexVal = pop();
			const Value& arrayVal = pop();

			auto& arr = getArray(arrayVal);

			int idxNum = getInt(indexVal);

			size_t index = static_cast<size_t>(idxNum);
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");

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

		case OpCode::CALL:
		{
			uint16_t fnIndex = readUint16();
			uint16_t args = readUint16();
			callFunction(fnIndex, args);
			break;
		}

		case OpCode::CALL_BINARY:
		{
			uint16_t opIndex = readUint16();
			callBinaryOperator(opIndex);
			break;
		}
		
		case OpCode::PRINT:
		{
			uint16_t argCount = readUint16();
			std::vector<Value> args;
			for (int i = 0; i < argCount; ++i)
			{
				args.push_back(pop());
			}

			auto it = m_program.m_specialFns.find("print");
			if (it != m_program.m_specialFns.end())
			{
				it->second(args.data(), argCount);
			}

			// leave empty value on stack after print (for chaining)
			push(Value());

			break;
		}

		case OpCode::RETURN:
		{
			Value result;

			if (stackTop > globalCount)
				result = pop();

			CallFrame frame = m_frames.back();
			m_frames.pop_back();

			stackTop = frame.base;
			push(result);

			if (m_debug)
			{
				printStack(m_stack, globalCount, (int)stackTop);
			}

			if (m_frames.empty())
				return result;

			break;
		}

		case OpCode::POP:
		{
			pop();
			break;
		}

		default:
			throw std::runtime_error("Unknown opcode");
		}

		if (m_debug && (instruction != (uint8_t)OpCode::RETURN))
		{
			printStack(m_stack, globalCount, (int)stackTop);
		}
	}

	throw std::runtime_error("Unexpected end of code.");
}

void VM::callFunction(int fnIndex, int argCount)
{
	const FunctionInfo& fn = m_program.functions[fnIndex];

	if (fn.isNative)
	{
		Value result = fn.fnc(&m_stack[stackTop - argCount], argCount);
		stackTop -= argCount;
		push(result);
		return;
	}

	if (argCount != fn.args.size())
		throw std::runtime_error("Arity mismatch in call to " + fn.name);

	if (m_frames.size() > MAX_CALL_DEPTH)
		throw std::runtime_error("Stack overflow: too many nested function calls.");

	CallFrame frame;
	frame.functionIndex = fnIndex;
	frame.ip = fn.entry;
	frame.base = stackTop - argCount;

	m_frames.push_back(frame);
}

void VM::callBinaryOperator(int opIndex)
{
	const BinaryOperatorInfo& fn = m_program.operators[opIndex];

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
