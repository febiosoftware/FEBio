#include "tools.h"
#include "compiler.h"
#include <iomanip>
#include <bitset>

using namespace std;

void febcode::writeProgramOpCodes(std::ostream& os, const Program& prg)
{
	os << "Stack Size : " << prg.maxStackSize << "\n";
	os << "Constants:\n";
	for (size_t i = 0; i < prg.constants.size(); ++i)
	{
		os << "  [" << i << "] " << prg.constants[i] << "\n";
	}

	os << "\nGlobals:\n";
	for (const auto& [name, index] : prg.globalIndices)
	{
		int slot = prg.globals[index].slot;
		os << "  [" << slot << "] " << name << "\n";
	}

	os << "\nBytecode:\n";
	const std::vector<uint8_t>& code = prg.code;
	for (int i = 0; i < code.size(); ++i)
	{
		os << std::setw(4) << i << ": ";

		febcode::OpCode opcode = (febcode::OpCode)code[i];

		const char* opcodeName = febcode::OpCodeToString(opcode);

		switch (opcode)
		{
		case febcode::OpCode::ADD_INT: os << opcodeName; break;
		case febcode::OpCode::ADD_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::SUB_INT: os << opcodeName; break;
		case febcode::OpCode::SUB_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_INT: os << opcodeName; break;
		case febcode::OpCode::MUL_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::DIV_INT: os << opcodeName; break;
		case febcode::OpCode::DIV_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::EXP_INT: os << opcodeName; break;
		case febcode::OpCode::EXP_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::SQR_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::SQRT_DOUBLE: os << opcodeName; break;

		case febcode::OpCode::EQUAL_BOOL:
		case febcode::OpCode::EQUAL_INT:
		case febcode::OpCode::EQUAL_DOUBLE:
			os << opcodeName; break;

		case febcode::OpCode::NEQ_BOOL:
		case febcode::OpCode::NEQ_INT:
		case febcode::OpCode::NEQ_DOUBLE:
			os << opcodeName; break;

		case febcode::OpCode::GT_INT: os << opcodeName; break;
		case febcode::OpCode::GT_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::LT_INT: os << opcodeName; break;
		case febcode::OpCode::LT_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::GE_INT: os << opcodeName; break;
		case febcode::OpCode::GE_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::LE_INT: os << opcodeName; break;
		case febcode::OpCode::LE_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::NEG_INT: os << opcodeName; break;
		case febcode::OpCode::NEG_DOUBLE: os << opcodeName; break;

		case febcode::OpCode::CREATE_VEC2_1ARG: os << opcodeName; break;
		case febcode::OpCode::GET_VEC2_X: os << opcodeName; break;
		case febcode::OpCode::GET_VEC2_Y: os << opcodeName; break;
		case febcode::OpCode::NEG_VEC2: os << opcodeName; break;
		case febcode::OpCode::ADD_VEC2: os << opcodeName; break;
		case febcode::OpCode::SUB_VEC2: os << opcodeName; break;
		case febcode::OpCode::DOT_VEC2: os << opcodeName; break;
		case febcode::OpCode::MUL_VEC2_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_DOUBLE_VEC2: os << opcodeName; break;
		case febcode::OpCode::DIV_VEC2_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_X: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_Y: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_Z: os << opcodeName; break;
		case febcode::OpCode::GET_VEC2_SWIZZLE:
		{
			i++;
			int mask = code[i++];
			int size = code[i];
			os << opcodeName << "," << size << "," << std::bitset<8>(mask);
			break;
		}
		case febcode::OpCode::GET_VEC3_SWIZZLE:
		{
			i++;
			int mask = code[i++];
			int size = code[i];
			os << opcodeName << "," << size << "," << std::bitset<8>(mask);
			break;
		}
		case febcode::OpCode::GET_VEC2_INDEX:
		{
			i++;
			os << opcodeName << " [" << (int)code[i] << "]";
			break;
		}
		case febcode::OpCode::GET_VEC3_INDEX:
		{
			i++;
			os << opcodeName << " [" << (int)code[i] << "]";
			break;
		}
		case febcode::OpCode::GET_MAT2_INDEX:
		{
			i++;
			os << opcodeName << " [" << (int)code[i] << "]";
			break;
		}

		case febcode::OpCode::CREATE_VEC3_1ARG: os << opcodeName; break;
		case febcode::OpCode::NEG_VEC3: os << opcodeName; break;
		case febcode::OpCode::ADD_VEC3: os << opcodeName; break;
		case febcode::OpCode::SUB_VEC3: os << opcodeName; break;
		case febcode::OpCode::DOT_VEC3: os << opcodeName; break;
		case febcode::OpCode::MUL_VEC3_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_DOUBLE_VEC3: os << opcodeName; break;
		case febcode::OpCode::DIV_VEC3_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::NOT: os << opcodeName; break;
		case febcode::OpCode::CREATE_MAT3_DIAG: os << opcodeName; break;
		case febcode::OpCode::CREATE_MAT3_VEC3: os << opcodeName; break;

		case febcode::OpCode::POP_VOID: os << opcodeName; break;
		case febcode::OpCode::POP_BOOL: os << opcodeName; break;
		case febcode::OpCode::POP_INT: os << opcodeName; break;
		case febcode::OpCode::POP_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::POP_VEC2: os << opcodeName; break;
		case febcode::OpCode::POP_VEC3: os << opcodeName; break;
		case febcode::OpCode::POP_MAT2: os << opcodeName; break;
		case febcode::OpCode::POP_MAT3: os << opcodeName; break;

		case febcode::OpCode::ADD_MAT2: os << opcodeName; break;
		case febcode::OpCode::SUB_MAT2: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT2: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT2_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::DIV_MAT2_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT2_VEC2: os << opcodeName; break;

		case febcode::OpCode::ADD_MAT3: os << opcodeName; break;
		case febcode::OpCode::SUB_MAT3: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT3: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT3_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_DOUBLE_MAT3: os << opcodeName; break;
		case febcode::OpCode::DIV_MAT3_DOUBLE: os << opcodeName; break;
		case febcode::OpCode::MUL_MAT3_VEC3: os << opcodeName; break;

		case febcode::OpCode::ADD_GLOBAL_MAT3:
		case febcode::OpCode::SUB_GLOBAL_MAT3:
		case febcode::OpCode::MUL_GLOBAL_MAT3:
		{
			i++;
			int a = code[i++];
			int b = code[i];
			os << opcodeName << " [" << a << "] [" << b << ']';
			break;
		}

		case febcode::OpCode::POP_ARRAY:
		{
			i++;
			int size = code[i];
			os << opcodeName << '[' << size << ']';
			break;
		}

		case febcode::OpCode::POP_STRUCT:
		{
			i++;
			int size = code[i];
			os << opcodeName << '[' << size << ']';
			break;
		}

		case febcode::OpCode::RETURN_VOID:
		case febcode::OpCode::RETURN_BOOL:
		case febcode::OpCode::RETURN_INT:
		case febcode::OpCode::RETURN_DOUBLE:
		case febcode::OpCode::RETURN_VEC2:
		case febcode::OpCode::RETURN_VEC3:
		case febcode::OpCode::RETURN_MAT2:
		case febcode::OpCode::RETURN_MAT3:
		case febcode::OpCode::RETURN_STRUCT:
			os << opcodeName; break;
		case febcode::OpCode::RETURN_ARRAY:
		{
			i++;
			int size = code[i];
			os << opcodeName << '[' << size << ']';
			break;
		}

		case febcode::OpCode::GET_INDEX_BOOL:
		case febcode::OpCode::GET_INDEX_INT:
		case febcode::OpCode::GET_INDEX_DOUBLE:
		case febcode::OpCode::GET_INDEX_VEC2:
		case febcode::OpCode::GET_INDEX_VEC3:
		case febcode::OpCode::GET_INDEX_ARRAY:
		case febcode::OpCode::GET_INDEX_STRUCT:
			os << opcodeName; break;

		case febcode::OpCode::GET_VEC2_X_REF: os << opcodeName; break;
		case febcode::OpCode::GET_VEC2_Y_REF: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_X_REF: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_Y_REF: os << opcodeName; break;
		case febcode::OpCode::GET_VEC3_Z_REF: os << opcodeName; break;

		case febcode::OpCode::STORE_BOOL:
		case febcode::OpCode::STORE_INT:
		case febcode::OpCode::STORE_DOUBLE:
		case febcode::OpCode::STORE_VEC2:
		case febcode::OpCode::STORE_VEC3:
		case febcode::OpCode::STORE_MAT2:
		case febcode::OpCode::STORE_MAT3:
		case febcode::OpCode::STORE_ARRAY:
		case febcode::OpCode::STORE_STRUCT:
			os << opcodeName; break;

		case febcode::OpCode::PUSH_VOID: os << opcodeName; break;
		case febcode::OpCode::PUSH_BOOL:
		case febcode::OpCode::PUSH_INT:
		case febcode::OpCode::PUSH_DOUBLE:
		case febcode::OpCode::PUSH_VEC2:
		case febcode::OpCode::PUSH_VEC3:
		case febcode::OpCode::PUSH_MAT2:
		case febcode::OpCode::PUSH_MAT3:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			const Value& v = prg.constants[offset];
			os << " [" << offset << "] (=" << v << ")";
		}
		break;
		case febcode::OpCode::PUSH_LOCAL:
		{
			os << opcodeName;
			i++;
			int size = code[i];
			os << " " << size;
			break;
		}
		case febcode::OpCode::GET_GLOBAL_BOOL:
		case febcode::OpCode::GET_GLOBAL_INT:
		case febcode::OpCode::GET_GLOBAL_DOUBLE:
		case febcode::OpCode::GET_GLOBAL_VEC2:
		case febcode::OpCode::GET_GLOBAL_VEC3:
		case febcode::OpCode::GET_GLOBAL_MAT2:
		case febcode::OpCode::GET_GLOBAL_MAT3:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::GET_GLOBAL_ARRAY:
		{
			os << opcodeName;
			i++;
			int size = code[i++];
			int slot = code[i];
			os << " " << slot << '[' << size << ']';
			break;
		}
		case febcode::OpCode::GET_GLOBAL_STRUCT:
		{
			os << opcodeName;
			i++;
			int size = code[i++];
			int slot = code[i];
			os << " " << slot << '[' << size << ']';
			break;
		}
		case febcode::OpCode::SET_GLOBAL_BOOL:
		case febcode::OpCode::SET_GLOBAL_INT:
		case febcode::OpCode::SET_GLOBAL_DOUBLE:
		case febcode::OpCode::SET_GLOBAL_VEC2:
		case febcode::OpCode::SET_GLOBAL_VEC3:
		case febcode::OpCode::SET_GLOBAL_MAT2:
		case febcode::OpCode::SET_GLOBAL_MAT3:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::SET_GLOBAL_ARRAY:
		case febcode::OpCode::SET_GLOBAL_STRUCT:
		{
			os << opcodeName;
			i++;
			int size = code[i++];
			int slot = code[i];
			os << " " << slot << '[' << size << ']';
			break;
		}

		case febcode::OpCode::GET_LOCAL_BOOL:
		case febcode::OpCode::GET_LOCAL_INT:
		case febcode::OpCode::GET_LOCAL_DOUBLE:
		case febcode::OpCode::GET_LOCAL_VEC2:
		case febcode::OpCode::GET_LOCAL_VEC3:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::GET_LOCAL_ARRAY:
		case febcode::OpCode::GET_LOCAL_STRUCT:
		{
			os << opcodeName;
			i++;
			int size = code[i++];
			int slot = code[i];
			os << " " << slot << '[' << size << ']';
			break;

		}
		case febcode::OpCode::GET_PROPERTY_BOOL:
		case febcode::OpCode::GET_PROPERTY_INT:
		case febcode::OpCode::GET_PROPERTY_DOUBLE:
		case febcode::OpCode::GET_PROPERTY_VEC2:
		case febcode::OpCode::GET_PROPERTY_VEC3:
		case febcode::OpCode::GET_PROPERTY_ARRAY:
		case febcode::OpCode::GET_PROPERTY_STRUCT:
		{
			os << opcodeName;
			i++;
			int type = code[i++];
			int slot = code[i];
			os << " " << slot;
			break;
		}

		case febcode::OpCode::CALL:
		{
			os << opcodeName;
			i++;
			int offset = code[i++];
			int args = code[i];
			std::string fnname = prg.functions[offset].name;
			int jmp = (int)prg.functions[offset].entry;
			os << " " << offset << " (-->" << jmp << ":" << fnname << ")";
			break;
		}
		case febcode::OpCode::JUMP:
		{
			os << opcodeName;
			i++;
			int offset = (code[i++] << 8) | code[i];
			os << " " << offset << " (-->" << i + offset + 1 << ")";
			break;
		}
		case febcode::OpCode::JUMP_IF_FALSE:
		{
			os << opcodeName;
			i++;
			int offset = (code[i++] << 8) | code[i];
			os << " " << offset << " (-->" << i + offset + 1 << ")";
			break;
		}
		case febcode::OpCode::JUMP_IF_TRUE:
		{
			os << opcodeName;
			i++;
			int offset = (code[i++] << 8) | code[i];
			os << " " << offset << " (-->" << i + offset + 1 << ")";
			break;
		}
		case febcode::OpCode::LOOP:
		{
			os << opcodeName;
			i++;
			int offset = (code[i++] << 8) | code[i];
			os << "," << offset;
			break;
		}
		case febcode::OpCode::GET_GLOBAL_REF:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::GET_LOCAL_REF:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::GET_MEMBER_REF:
		{
			os << opcodeName;
			i++;
			int offset = code[i];
			os << " " << offset;
			break;
		}
		case febcode::OpCode::GET_INDEX_REF_BOOL:
		case febcode::OpCode::GET_INDEX_REF_INT:
		case febcode::OpCode::GET_INDEX_REF_DOUBLE:
		case febcode::OpCode::GET_INDEX_REF_VEC2:
		case febcode::OpCode::GET_INDEX_REF_VEC3:
			os << opcodeName; break;
		case febcode::OpCode::GET_INDEX_REF:
		{
			os << opcodeName;
			i++;
			int size = code[i];
			os << " " << size;
			break;
		}
		default:
			os << "UNKNOWN(" << (int)code[i] << ")\n";
			return;
			break;
		}
		os << std::endl;
	}
}
