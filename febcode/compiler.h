#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <variant>
#include <cstdint>
#include <stdexcept>
#include <functional>
#include "program.h"

namespace febcode
{
	enum class OpCode : uint8_t
	{
		PUSH_CONST,

		GET_GLOBAL,
		SET_GLOBAL,

		GET_LOCAL,
		SET_LOCAL,

		// struct codes
		CREATE_STRUCT,
		COPY_STRUCT,
		GET_PROPERTY,
		SET_PROPERTY,

		// array codes
		CREATE_ARRAY,
		COPY_ARRAY,
		GET_INDEX,
		SET_INDEX,

		// vec2 codes
		CREATE_VEC2,
		COPY_VEC2,
		GET_VEC2_X,
		GET_VEC2_Y,
		SET_VEC2_X,
		SET_VEC2_Y,

		// vec3 codes
		CREATE_VEC3,
		COPY_VEC3,
		GET_VEC3_X,
		GET_VEC3_Y,
		GET_VEC3_Z,
		GET_VEC3_SWIZZLE,
		SET_VEC3_X,
		SET_VEC3_Y,
		SET_VEC3_Z,

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

		// logical operators
		NOT,

		// string operators
		ADD_STRING,

		// comparison operators
		EQUAL,
		NOT_EQUAL,

		// control flow
		JUMP,
		JUMP_IF_FALSE,
		JUMP_IF_TRUE,
		LOOP,

		CALL,		// call function
		CALL_BINARY, // call binary operator
		RETURN,
		PRINT,

		POP
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
			bool isInitialized = false;
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
		Type compileSetMember    (SetMemberExpr* expr);
		Type compileIndex        (IndexExpr* expr);
		Type compileSetIndex     (SetIndexExpr* expr);

		void compileInitializer(Expression* expr, Type expectedType);

		std::vector<Type> compileFncArgs(std::vector<std::unique_ptr<Expression>>& args, bool copyArgs = true);

		// ===== Bytecode =====

		void emit(OpCode op);
		void emitUint16(uint16_t v);

		uint16_t addConstant(const Value& v);

		int emitJump(OpCode op);
		void patchJump(int offset);
		void emitLoop(int loopStart);

	private:
		Program& prg;

		Type expectedReturnType = nullptr;

		std::vector<Local> m_locals;
		int m_scopeDepth = 0;
	};

	void CompileSource(Program& prg, const std::string& source);
} // namespace febcode
