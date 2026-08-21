#pragma once
#include <memory>
#include <unordered_map>
#include <string>
#include "ast.h"
#include "modifier.h"
#include "simplifier.h"

namespace febcode {
	class Differentiator : public Modifier {

	public:
		enum DiffMode {
			SCALAR,		// scalar differentiation (for scalar variables or components of non-scalar variables)
			GRADIENT,	// gradient differentiation (for vector variables)
			DIRECTIONAL // directional differentiation (for directional derivatives)
		};

	private:
		struct DerivVar
		{
			DiffMode mode; // differentiation mode
			std::string name; // name of the derivative variable
			Type type; // type of the derivative variable
			int component = -1; // -1 for full derivative, or component index for partial derivative of nonscalar variables
			std::string tangentVar; // name of the tangent variable for directional differentiation
		};

	public:
		Differentiator(Program& prg) : Modifier(prg), simplifier(prg) {}

		void differentiate(DiffMode dm, const std::string& var, int component = -1, const std::string& tangentVar = "");

		bool DependencyFound() const { return dependencyFound; }

		void SetSimplify(bool value) { doSimplify = value; }

	private:

		// differentiate an AST to produce a new AST representing the derivative
		std::unique_ptr<AST> differentiate(const AST& ast, const DerivVar& var);

		void differentiateStmt(BlockStmt& ast, Statement* stmt, const DerivVar& var);

		void diffExpressionStmt(BlockStmt& ast, ExpressionStmt* stmt, const DerivVar& var);
		void diffReturnStmt    (BlockStmt& ast, ReturnStmt*     stmt, const DerivVar& var);
		void diffStructStmt    (BlockStmt& ast, StructStmt*     stmt, const DerivVar& var);
		void diffVarDeclStmt   (BlockStmt& ast, VarDeclStmt*    stmt, const DerivVar& var);
		void diffIfStmt        (BlockStmt& ast, IfStmt*         stmt, const DerivVar& var);
		void diffBlockStmt     (BlockStmt& ast, BlockStmt*      stmt, const DerivVar& var);

	private:
		// Differentiate an expression with respect to a variable
		std::unique_ptr<Expression> differentiate(const Expression* expr, const DerivVar& var);

		std::unique_ptr<Expression> diffLiteral    (const LiteralExpr*     literal , const DerivVar& var);
		std::unique_ptr<Expression> diffVariable   (const VariableExpr*    variable, const DerivVar& var);
		std::unique_ptr<Expression> diffUnary      (const UnaryExpr*       unary   , const DerivVar& var);
		std::unique_ptr<Expression> diffBinary     (const BinaryExpr*      binary  , const DerivVar& var);
		std::unique_ptr<Expression> diffCall       (const CallExpr*        call    , const DerivVar& var);
		std::unique_ptr<Expression> diffInit       (const InitExpr*        init    , const DerivVar& var);
		std::unique_ptr<Expression> diffConstructor(const ConstructorExpr* ctor    , const DerivVar& var);
		std::unique_ptr<Expression> diffAssign     (const AssignExpr*      assign  , const DerivVar& var);
		std::unique_ptr<Expression> diffIndex      (const IndexExpr*       index   , const DerivVar& var);
		std::unique_ptr<Expression> diffMember     (const MemberExpr*      member  , const DerivVar& var);

		Type getDerivativeType(Type varType, const DerivVar& dvar);

		ExprPtr simplify(const Expression* expr) 
		{ 
			if (doSimplify)
				return simplifier.simplify(expr); 
			else
				return clone(expr);
		}
		ExprPtr simplify(const ExprPtr& expr) { return simplify(expr.get()); }

	private:
		bool dependencyFound = false; // flag to indicate if we found a dependency on the variable we're differentiating with respect to
		std::unordered_map<std::string, std::string> deriveVars; // map of derivative variables
		std::unordered_map<std::string, Type> varTypes; // map of variable types for variables in the original program, used to determine the type of the derivative variables.

		bool doSimplify = true;
		Simplifier simplifier;
	};

	bool dependsOn(const Statement* stmt, const std::string& varName);
	bool dependsOn(const Expression* expr, const std::string& varName);
} // namespace febcode
