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
	case OpCode::SET_VEC2_X    : return "SV2X";
	case OpCode::SET_VEC2_Y    : return "SV2Y";
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
	case OpCode::SET_VEC3_X    : return "SV3X";
	case OpCode::SET_VEC3_Y    : return "SV3Y";
	case OpCode::SET_VEC3_Z    : return "SV3Z";
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
	case OpCode::SET_PROPERTY  : return "SETP";
	case OpCode::CREATE_ARRAY  : return "ARR ";
	case OpCode::COPY_ARRAY    : return "CPYA";
	case OpCode::SET_INDEX     : return "SETI";
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
	}
	return "(UNKNOWN)";
}

void printStack(const std::vector<Value>& stack, int numGlobals)
{
	std::cout << "Stack: [";
	for (size_t i = 0; i < numGlobals; i++)
	{
		std::cout << ValueToString(stack[i]);
		if (i < numGlobals - 1)
			std::cout << ",";
	}
	std::cout << "|";
	for (size_t i = numGlobals; i < stack.size(); i++)
	{
		std::cout << ValueToString(stack[i]);
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
			m_stack.push_back(m_program.constants[idx]);
			break;
		}

		case OpCode::GET_GLOBAL:
		{
			uint16_t slot = readUint16();
			m_stack.push_back(m_stack[slot]);
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
			m_stack.push_back(m_stack[frame.base + slot]);
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
			Value a = pop();
			m_stack.push_back(-std::get<int>(a));
			break;
		}
		case OpCode::ADD_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) + getInt(b));
			break;
		}
		case OpCode::SUB_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) - getInt(b));
			break;
		}
		case OpCode::MUL_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) * getInt(b));
			break;
		}
		case OpCode::DIV_INT:
		{
			Value b = pop();
			Value a = pop();
			int divisor = getInt(b);
			if (divisor == 0)
				throw std::runtime_error("division by zero.");
			m_stack.push_back(getInt(a) / divisor);
			break;
		}
		case OpCode::EXP_INT:
		{
			Value b = pop();
			Value a = pop();

			int e = getInt(b);
			if (e < 0)
				throw std::runtime_error("Negative exponent not supported for integers.");

			m_stack.push_back(ipow(getInt(a), e));
			break;
		}
		case OpCode::GT_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) > getInt(b));
			break;
		}
		case OpCode::LT_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) < getInt(b));
			break;
		}
		case OpCode::GE_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) >= getInt(b));
			break;
		}
		case OpCode::LE_INT:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getInt(a) <= getInt(b));
			break;
		}

		// Double operators
		case OpCode::NEG_DOUBLE:
		{
			Value a = pop();
			m_stack.push_back(-std::get<double>(a));
			break;
		}

		case OpCode::ADD_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) + getDouble(b));
			break;
		}
		case OpCode::SUB_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) - getDouble(b));
			break;
		}
		case OpCode::MUL_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) * getDouble(b));
			break;
		}
		case OpCode::DIV_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			double divisor = getDouble(b);
			if (divisor == 0.0)
				throw std::runtime_error("division by zero.");
			m_stack.push_back(getDouble(a) / divisor);
			break;
		}
		case OpCode::EXP_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(std::pow(getDouble(a), getDouble(b)));
			break;
		}
		case OpCode::GT_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) > getDouble(b));
			break;
		}
		case OpCode::LT_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) < getDouble(b));
			break;
		}
		case OpCode::GE_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) >= getDouble(b));
			break;
		}
		case OpCode::LE_DOUBLE:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getDouble(a) <= getDouble(b));
			break;
		}

		case OpCode::CREATE_VEC2:
		{
			Value y = pop();
			Value x = pop();
			m_stack.push_back(vec2(getDouble(x), getDouble(y)));
			break;
		}
		case OpCode::COPY_VEC2:
		{
			Value src = pop();
			const vec2& v = getVec2(src);
			m_stack.push_back(vec2(v.x, v.y));
			break;
		}
		case OpCode::ADD_VEC2:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec2(a) + getVec2(b));
			break;
		}
		case OpCode::SUB_VEC2:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec2(a) - getVec2(b));
			break;
		}
		case OpCode::DOT_VEC2:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec2(a) * getVec2(b));
			break;
		}
		case OpCode::MUL_VEC2_DOUBLE:
		{
			Value scalar = pop();
			Value vec = pop();
			m_stack.push_back(getVec2(vec) * getDouble(scalar));
			break;
		}
		case OpCode::MUL_DOUBLE_VEC2:
		{
			Value vec = pop();
			Value scalar = pop();
			m_stack.push_back(getVec2(vec) * getDouble(scalar));
			break;
		}
		case OpCode::GET_VEC2_X:
		{
			Value vec = pop();
			m_stack.push_back(getVec2(vec).x);
			break;
		}
		case OpCode::GET_VEC2_Y:
		{
			Value vec = pop();
			m_stack.push_back(getVec2(vec).y);
			break;
		}
		case OpCode::SET_VEC2_X:
		{
			uint16_t slot = readUint16();
			Value x = pop();
			Value vec = pop();
			vec2& v = getVec2(vec);
			v.x = getDouble(x);
			m_stack[slot] = v;
			m_stack.push_back(v.x);
			break;
		}
		case OpCode::SET_VEC2_Y:
		{
			uint16_t slot = readUint16();
			Value y = pop();
			Value vec = pop();
			vec2& v = getVec2(vec);
			v.y = getDouble(y);
			m_stack[slot] = v;
			m_stack.push_back(v.y);
			break;
		}
		case OpCode::NEG_VEC2:
		{
			Value vec = pop();
			vec2& v = getVec2(vec);
			m_stack.push_back(vec2(-v.x, -v.y));
			break;
		}
		case OpCode::CREATE_VEC3:
		{
			Value z = pop();
			Value y = pop();
			Value x = pop();
			m_stack.push_back(vec3(getDouble(x), getDouble(y), getDouble(z)));
			break;
		}
		case OpCode::COPY_VEC3:
		{
			Value src = pop();
			const vec3& v = getVec3(src);
			m_stack.push_back(vec3(v.x, v.y, v.z));
			break;
		}
		case OpCode::ADD_VEC3:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec3(a) + getVec3(b));
			break;
		}
		case OpCode::SUB_VEC3:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec3(a) - getVec3(b));
			break;
		}
		case OpCode::DOT_VEC3:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(getVec3(a) * getVec3(b));
			break;
		}
		case OpCode::MUL_VEC3_DOUBLE:
		{
			Value scalar = pop();
			Value vec = pop();
			m_stack.push_back(getVec3(vec) * getDouble(scalar));
			break;
		}
		case OpCode::MUL_DOUBLE_VEC3:
		{
			Value vec = pop();
			Value scalar = pop();
			m_stack.push_back(getVec3(vec) * getDouble(scalar));
			break;
		}
		case OpCode::GET_VEC3_X:
		{
			Value vec = pop();
			m_stack.push_back(getVec3(vec).x);
			break;
		}
		case OpCode::GET_VEC3_Y:
		{
			Value vec = pop();
			m_stack.push_back(getVec3(vec).y);
			break;
		}
		case OpCode::GET_VEC3_Z:
		{
			Value vec = pop();
			m_stack.push_back(getVec3(vec).z);
			break;
		}
		case OpCode::SET_VEC3_X:
		{
			uint16_t slot = readUint16();
			Value x = pop();
			Value vec = pop();
			vec3& v = getVec3(vec);
			v.x = getDouble(x);
			m_stack[slot] = v;
			m_stack.push_back(v.x);
			break;
		}
		case OpCode::SET_VEC3_Y:
		{
			uint16_t slot = readUint16();
			Value y = pop();
			Value vec = pop();
			vec3& v = getVec3(vec);
			v.y = getDouble(y);
			m_stack[slot] = v;
			m_stack.push_back(v.y);
			break;
		}
		case OpCode::SET_VEC3_Z:
		{
			uint16_t slot = readUint16();
			Value z = pop();
			Value vec = pop();
			vec3& v = getVec3(vec);
			v.z = getDouble(z);
			m_stack[slot] = v;
			m_stack.push_back(v.z);
			break;
		}
		case OpCode::NEG_VEC3:
		{
			Value vec = pop();
			vec3& v = getVec3(vec);
			m_stack.push_back(vec2(-v.x, -v.y));
			break;
		}

		// string operators
		case OpCode::ADD_STRING:
		{
			Value a = pop();
			Value b = pop();
			m_stack.push_back(getString(a) + getString(b));
			break;
		}

		// Logical operators
		case OpCode::NOT:
		{
			Value a = pop();
			bool b = isTruthy(a);
			m_stack.push_back(!b);
			break;
		}
		case OpCode::EQUAL:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(a == b);
			break;
		}
		case OpCode::NOT_EQUAL:
		{
			Value b = pop();
			Value a = pop();
			m_stack.push_back(a != b);
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

			m_stack.push_back(obj);
			break;
		}
		case OpCode::COPY_STRUCT:
		{
			auto src = pop();

			if (!isStruct(src))
				throw std::runtime_error("COPY_STRUCT operand must be a struct");

			const StructValue& srcObj = getStruct(src);
			auto obj = copyStruct(srcObj);

			m_stack.push_back(obj);
			break;
		}
		case OpCode::GET_PROPERTY:
		{
			uint16_t slot = readUint16();

			Value objVal = pop();

			switch (ValueType(objVal))
			{
			case TypeKind::Struct:
			{
				const auto& structVal = getStruct(objVal);

				if (slot >= structVal.fields.size())
					throw std::runtime_error("Invalid property slot");

				m_stack.push_back(structVal.fields[slot]);
				break;
			}
			default:
				throw std::runtime_error("GET_PROPERTY on non-struct");
			}

			break;
		}

		case OpCode::SET_PROPERTY:
		{
			uint16_t slot = readUint16();

			Value value = pop();
			Value objVal = pop();

			switch (ValueType(objVal))
			{
			case TypeKind::Struct:
			{
				auto& structVal = getStruct(objVal);

				if (slot >= structVal.fields.size())
					throw std::runtime_error("Invalid property slot");

				structVal.fields[slot] = value;

				// Push assigned value (JS-style assignment semantics)
				m_stack.push_back(value);
				break;
			}
			default:
				throw std::runtime_error("SET_PROPERTY on non-struct");
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

			m_stack.push_back(arr);
			break;
		}
		case OpCode::COPY_ARRAY: {
			auto src = pop();

			if (!isArray(src))
				throw std::runtime_error("COPY_ARRAY operand must be an array");

			const ArrayValue& srcObj = getArray(src);
			auto obj = copyArray(srcObj);

			m_stack.push_back(obj);
			break;
		}
		case OpCode::GET_INDEX: {
			Value indexVal = pop();
			Value arrayVal = pop();

			auto& arr = getArray(arrayVal);

			int idxNum = std::get<int>(indexVal);

			size_t index = static_cast<size_t>(idxNum);
			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");

			m_stack.push_back(arr.elements[index]);
			break;
		}
		case OpCode::SET_INDEX: {
			Value value = pop();
			Value indexVal = pop();
			Value arrayVal = pop();

			auto& arr = getArray(arrayVal);

			int idxNum = std::get<int>(indexVal);

			size_t index = static_cast<size_t>(idxNum);

			if (index >= arr.size())
				throw std::runtime_error("Array index out of bounds.");

			arr.elements[index] = value;

			m_stack.push_back(value);
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
				it->second(args);
			}

			// leave empty value on stack after print (for chaining)
			m_stack.push_back(std::monostate{});

			break;
		}

		case OpCode::RETURN:
		{
			Value result;

			if (m_stack.size() > globalCount)
				result = pop();

			CallFrame frame = m_frames.back();
			m_frames.pop_back();

			m_stack.resize(frame.base);
			m_stack.push_back(result);

			if (m_debug)
			{
				printStack(m_stack, globalCount);
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
			printStack(m_stack, globalCount);
		}
	}

	throw std::runtime_error("Unexpected end of code.");
}

void VM::callFunction(int fnIndex, int argCount)
{
	const FunctionInfo& fn = m_program.functions[fnIndex];

	if (fn.isNative)
	{
		std::vector<Value> args;
		for (int i = 0; i < argCount; ++i)
			args.push_back(m_stack[m_stack.size() - argCount + i]);

		m_stack.resize(m_stack.size() - argCount);

		Value result = fn.fnc(args);
		m_stack.push_back(result);
		return;
	}

	if (argCount != fn.args.size())
		throw std::runtime_error("Arity mismatch in call to " + fn.name);

	if (m_frames.size() > MAX_CALL_DEPTH)
		throw std::runtime_error("Stack overflow: too many nested function calls.");

	CallFrame frame;
	frame.functionIndex = fnIndex;
	frame.ip = fn.entry;
	frame.base = m_stack.size() - argCount;

	m_frames.push_back(frame);
}

void VM::callBinaryOperator(int opIndex)
{
	const BinaryOperatorInfo& fn = m_program.operators[opIndex];

	Value lv = m_stack[m_stack.size() - 2];
	Value rv = m_stack[m_stack.size() - 1];
	m_stack.resize(m_stack.size() - 2);

	Value result = fn.fnc(lv, rv);
	m_stack.push_back(result);
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
