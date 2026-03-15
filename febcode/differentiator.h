#pragma once
#include <memory>
#include <unordered_map>
#include <string>
#include "ast.h"

namespace febcode {
	class Differentiator {
	public:
		// differentiate an AST to produce a new AST representing the derivative
		std::unique_ptr<AST> differentiate(const AST& ast, const std::string& var);

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
		std::unordered_map<std::string, std::string> deriveVars; // map of derivative variables
	};
} // namespace febcode
