#pragma once
#include <memory>
#include <unordered_map>
#include <string>
#include "ast.h"
#include "program.h"

namespace febcode {
	class Differentiator {
	public:
		Differentiator(Program& prg) : prg(prg) {}

		// differentiate an AST to produce a new AST representing the derivative
		std::unique_ptr<AST> differentiate(const AST& ast, const std::string& var);

		bool DependencyFound() const { return dependencyFound; }

	private:
		void differentiateStmt(BlockStmt& ast, Statement* stmt, const std::string& var);

		void diffExpressionStmt(BlockStmt& ast, ExpressionStmt* stmt, const std::string& var);
		void diffReturnStmt    (BlockStmt& ast, ReturnStmt*     stmt, const std::string& var);
		void diffStructStmt    (BlockStmt& ast, StructStmt*     stmt, const std::string& var);
		void diffVarDeclStmt   (BlockStmt& ast, VarDeclStmt*    stmt, const std::string& var);
		void diffIfStmt        (BlockStmt& ast, IfStmt*         stmt, const std::string& var);

	private:
		// Differentiate an expression with respect to a variable
		std::unique_ptr<Expression> differentiate(const Expression* expr, const std::string& var);

		std::unique_ptr<Expression> diffLiteral (const LiteralExpr*     literal , const std::string& var);
		std::unique_ptr<Expression> diffVariable(const VariableExpr*    variable, const std::string& var);
		std::unique_ptr<Expression> diffUnary   (const UnaryExpr*       unary   , const std::string& var);
		std::unique_ptr<Expression> diffBinary  (const BinaryExpr*      binary  , const std::string& var);
		std::unique_ptr<Expression> diffCall    (const CallExpr*        call    , const std::string& var);
		std::unique_ptr<Expression> diffInit    (const InitializerExpr* init    , const std::string& var);
		std::unique_ptr<Expression> diffAssign  (const AssignExpr*      assign  , const std::string& var);
		std::unique_ptr<Expression> diffMember  (const MemberExpr*      member  , const std::string& var);

	private:
		bool dependencyFound = false; // flag to indicate if we found a dependency on the variable we're differentiating with respect to
		std::unordered_map<std::string, std::string> deriveVars; // map of derivative variables
		Program& prg;
	};

	bool dependsOn(const Statement* stmt, const std::string& varName);
} // namespace febcode
