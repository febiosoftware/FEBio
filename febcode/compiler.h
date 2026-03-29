#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "program.h"

namespace febcode
{
	enum OpCode : uint8_t
	{
		PUSH_VOID, // just push a dummy value, used for return without value
		PUSH_BOOL,
		PUSH_INT,
		PUSH_DOUBLE,

		GET_GLOBAL_BOOL,
		GET_GLOBAL_INT,
		GET_GLOBAL_DOUBLE,
		GET_GLOBAL_VEC2,
		GET_GLOBAL_VEC3,
		GET_GLOBAL_MAT2,
		GET_GLOBAL_MAT3,
		GET_GLOBAL_ARRAY,
		GET_GLOBAL_STRUCT,

		SET_GLOBAL_BOOL,
		SET_GLOBAL_INT,
		SET_GLOBAL_DOUBLE,
		SET_GLOBAL_VEC2,
		SET_GLOBAL_VEC3,
		SET_GLOBAL_MAT2,
		SET_GLOBAL_MAT3,
		SET_GLOBAL_ARRAY,
		SET_GLOBAL_STRUCT,

		GET_LOCAL_BOOL,
		GET_LOCAL_INT,
		GET_LOCAL_DOUBLE,
		GET_LOCAL_VEC2,
		GET_LOCAL_VEC3,
		GET_LOCAL_MAT2,
		GET_LOCAL_MAT3,
		GET_LOCAL_ARRAY,
		GET_LOCAL_STRUCT,

		GET_GLOBAL_REF,

		GET_LOCAL_REF,

		// struct codes
		CREATE_STRUCT,
		COPY_STRUCT,

		GET_PROPERTY_BOOL,
		GET_PROPERTY_INT,
		GET_PROPERTY_DOUBLE,
		GET_PROPERTY_VEC2,
		GET_PROPERTY_VEC3,
		GET_PROPERTY_MAT2,
		GET_PROPERTY_MAT3,
		GET_PROPERTY_ARRAY,
		GET_PROPERTY_STRUCT,

		GET_MEMBER_REF,

		GET_INDEX_BOOL,
		GET_INDEX_INT,
		GET_INDEX_DOUBLE,
		GET_INDEX_VEC2,
		GET_INDEX_VEC3,
		GET_INDEX_MAT2,
		GET_INDEX_MAT3,
		GET_INDEX_ARRAY,
		GET_INDEX_STRUCT,

		GET_GLOBAL_INDEX_DOUBLE,

		GET_INDEX_REF,
		GET_INDEX_REF_BOOL,
		GET_INDEX_REF_INT,
		GET_INDEX_REF_DOUBLE,
		GET_INDEX_REF_VEC2,
		GET_INDEX_REF_VEC3,
		GET_INDEX_REF_MAT2,
		GET_INDEX_REF_MAT3,

		// vec2 codes
		GET_VEC2_X,
		GET_VEC2_Y,
		GET_VEC2_X_REF,
		GET_VEC2_Y_REF,
		GET_VEC2_SWIZZLE,
		GET_VEC2_INDEX,

		// vec3 codes
		GET_VEC3_X,
		GET_VEC3_Y,
		GET_VEC3_Z,
		GET_VEC3_X_REF,
		GET_VEC3_Y_REF,
		GET_VEC3_Z_REF,
		GET_VEC3_SWIZZLE,
		GET_VEC3_INDEX,

		// int operators
		NEG_INT,
		ADD_INT,
		SUB_INT,
		MUL_INT,
		DIV_INT,
		EXP_INT,
		GT_INT,
		LT_INT,
		GE_INT,
		LE_INT,

		// double operators
		NEG_DOUBLE,
		ADD_DOUBLE,
		SUB_DOUBLE,
		MUL_DOUBLE,
		DIV_DOUBLE,
		EXP_DOUBLE,
		SQR_DOUBLE,
		SQRT_DOUBLE,
		GT_DOUBLE,
		LT_DOUBLE,
		GE_DOUBLE,
		LE_DOUBLE,

		// vec2 operators
		NEG_VEC2,
		ADD_VEC2,
		SUB_VEC2,
		DOT_VEC2,
		MUL_VEC2_DOUBLE,
		MUL_DOUBLE_VEC2,

		// vec3 operators
		NEG_VEC3,
		ADD_VEC3,
		SUB_VEC3,
		DOT_VEC3,
		MUL_VEC3_DOUBLE,
		MUL_DOUBLE_VEC3,

		// mat2 operators
		NEG_MAT2,
		ADD_MAT2,
		SUB_MAT2,
		MUL_MAT2,
		MUL_MAT2_DOUBLE,
		MUL_DOUBLE_MAT2,
		MUL_MAT2_VEC2,
		GET_MAT2_INDEX,

		// mat3 operators
		NEG_MAT3,
		ADD_MAT3,
		SUB_MAT3,
		MUL_MAT3,
		MUL_MAT3_DOUBLE,
		MUL_DOUBLE_MAT3,
		MUL_MAT3_VEC3,
		GET_MAT3_INDEX,

		ADD_GLOBAL_MAT3,
		SUB_GLOBAL_MAT3,
		MUL_GLOBAL_MAT3,

		// logical operators
		NOT,

		// comparison operators
		EQUAL_BOOL,
		EQUAL_INT,
		EQUAL_DOUBLE,

		NEQ_BOOL,
		NEQ_INT,
		NEQ_DOUBLE,

		// control flow
		JUMP,
		JUMP_IF_FALSE,
		JUMP_IF_TRUE,
		LOOP,

		// store value in variable (local or global)
		STORE_BOOL,
		STORE_INT,
		STORE_DOUBLE,
		STORE_VEC2,
		STORE_VEC3,
		STORE_MAT2,
		STORE_MAT3,
		STORE_ARRAY,
		STORE_STRUCT,

		POP_VOID,
		POP_BOOL,
		POP_INT,
		POP_DOUBLE,
		POP_VEC2,
		POP_VEC3,
		POP_MAT2,
		POP_MAT3,
		POP_ARRAY,
		POP_STRUCT,

		CALL,		// call function

		// returns (make sure these are the last opcodes, since they have special handling in the VM)
		RETURN_VOID,
		RETURN_BOOL,
		RETURN_INT,
		RETURN_DOUBLE,
		RETURN_VEC2,
		RETURN_VEC3,
		RETURN_MAT2,
		RETURN_MAT3,
		RETURN_ARRAY,
		RETURN_STRUCT,

		LAST_OPCODE, // not really an opcode, just a marker for the end of the enum. Also, make sure this is less than 255!
	};

	//
	// ================= COMPILER =================
	//

	class Compiler
	{
	public:
		Compiler(Program& program);

		void compile();

	private:

		// ===== Scope =====

		struct Local
		{
			std::string name;
			Type type;
			int depth;
			int slot;
		};

		void beginScope();
		void endScope();
		int resolveLocal(const std::string& name);
		int resolveGlobal(const std::string& name);
		int resolveFunction(const std::string& name, std::vector<Type> args);
		bool isNativeFunction(const std::string& name);

		// ===== Compile =====

		void compileStatement(Statement* stmt);
		void compileExprStmt (ExpressionStmt* stmt);
		void compileBlock    (BlockStmt* stmt);
		void compileVarDecl  (VarDeclStmt* decl);
		void compileFunction (FunctionStmt* fn);
		void compileIf       (IfStmt* stmt);
		void compileWhile    (WhileStmt* stmt);
		void compileReturn   (ReturnStmt* stmt);
		void compileStruct   (StructStmt* stmt);

		Type compileExpression   (Expression* expr);
		Type compileBinary       (BinaryExpr* expr);
		Type compileUnary        (UnaryExpr* expr);
		Type compileLiteral      (LiteralExpr* expr);
		Type compileVariable     (VariableExpr* expr);
		Type compileAssign       (AssignExpr* expr);
		Type compileCall         (CallExpr* expr);
		Type compileMember       (MemberExpr* expr);
		Type compileIndex        (IndexExpr* expr);

		Type compileLValue(Expression* expr);
		Type compileVariableRef(VariableExpr* expr);
		Type compileMemberRef(MemberExpr* expr);
		Type compileIndexRef(IndexExpr* expr);

		int resolveMember(Type type, const std::string& member);
		int resolveMemberOffset(Type type, const std::string& member);
		Type memberType(Type type, int memberIndex);

		void compileInitializer(Expression* expr, Type expectedType);

		std::vector<Type> compileFncArgs(std::vector<std::unique_ptr<Expression>>& args);

		// ===== Bytecode =====

		void emit(OpCode op, int arg = 0);
		void emitUint8(uint8_t v);
		void emitUint16(uint16_t v);

		int stackEffect(OpCode op, int arg);

		uint8_t addConstant(const Value& v);

		int emitJump(OpCode op);
		void patchJump(int offset);
		void emitLoop(int loopStart);

	private:
		Program& prg;

		bool hasReturn = false; // see if the program or function has an explicit return statement. If not, we will add a default return at the end.
		Type expectedReturnType = nullptr;

		std::vector<Local> m_locals;
		int m_scopeDepth = 0;
		int localStackSize = 0;

		int stackDepth = 0;
		int maxStackDepth = 0;

		int currentFunction = -1;
	};

	void CompileSource(Program& prg, const std::string& source);

	const char* OpCodeToString(febcode::OpCode ip);
} // namespace febcode
