#pragma once
#include <memory>
#include <unordered_map>
#include <string>
#include "ast.h"
#include "program.h"

namespace febcode {
	class Differentiator {

		struct deriveVar {
			TypeKind type = TypeKind::Double;
			std::string name;
		};

	public:
		Differentiator(Program& prg) : prg(prg) {}

		// differentiate an AST to produce a new AST representing the derivative
		std::unique_ptr<AST> differentiate(const AST& ast, const deriveVar& var);

		bool DependencyFound() const { return dependencyFound; }

	private:
		void differentiateStmt(BlockStmt& ast, Statement* stmt, const deriveVar& var);

		void diffExpressionStmt(BlockStmt& ast, ExpressionStmt* stmt, const deriveVar& var);
		void diffReturnStmt    (BlockStmt& ast, ReturnStmt*     stmt, const deriveVar& var);
		void diffStructStmt    (BlockStmt& ast, StructStmt*     stmt, const deriveVar& var);
		void diffVarDeclStmt   (BlockStmt& ast, VarDeclStmt*    stmt, const deriveVar& var);
		void diffIfStmt        (BlockStmt& ast, IfStmt*         stmt, const deriveVar& var);
		void diffBlockStmt     (BlockStmt& ast, BlockStmt*      stmt, const deriveVar& var);

	private:
		// Differentiate an expression with respect to a variable
		std::unique_ptr<Expression> differentiate(const Expression* expr, const deriveVar& var);

		std::unique_ptr<Expression> diffLiteral    (const LiteralExpr*     literal , const deriveVar& var);
		std::unique_ptr<Expression> diffVariable   (const VariableExpr*    variable, const deriveVar& var);
		std::unique_ptr<Expression> diffUnary      (const UnaryExpr*       unary   , const deriveVar& var);
		std::unique_ptr<Expression> diffBinary     (const BinaryExpr*      binary  , const deriveVar& var);
		std::unique_ptr<Expression> diffCall       (const CallExpr*        call    , const deriveVar& var);
		std::unique_ptr<Expression> diffInit       (const InitializerExpr* init    , const deriveVar& var);
		std::unique_ptr<Expression> diffConstructor(const ConstructorExpr* ctor    , const deriveVar& var);
		std::unique_ptr<Expression> diffAssign     (const AssignExpr*      assign  , const deriveVar& var);
		std::unique_ptr<Expression> diffIndex      (const IndexExpr*       index   , const deriveVar& var);
		std::unique_ptr<Expression> diffMember     (const MemberExpr*      member  , const deriveVar& var);

		Type getDerivativeType(Type varType, TypeKind derivType);

	private:
		bool dependencyFound = false; // flag to indicate if we found a dependency on the variable we're differentiating with respect to
		std::unordered_map<std::string, std::string> deriveVars; // map of derivative variables
		std::unordered_map<std::string, Type> varTypes; // map of variable types for variables in the original program, used to determine the type of the derivative variables.
		Program& prg;
	};

	bool dependsOn(const Statement* stmt, const std::string& varName);
} // namespace febcode
