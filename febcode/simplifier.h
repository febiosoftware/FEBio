#pragma once
#include "ast.h"
#include "modifier.h"
#include <functional>
#include <vector>

namespace febcode {

	class Simplifier : public Modifier {
	public:
		Simplifier(Program& prg);

		ExprPtr simplify(const Expression* expr);
		inline ExprPtr simplify(const std::unique_ptr<Expression>& expr) { return simplify(expr.get()); }

	private:
		ExprPtr simplifyConstructor(const ConstructorExpr* expr);
		ExprPtr simplifyUnary      (const UnaryExpr* expr);
		ExprPtr simplifyInitializer(const InitExpr* expr);
		ExprPtr simplifyCall       (const CallExpr* expr);
		ExprPtr simplifyBinary     (const BinaryExpr* expr);
		ExprPtr simplifyAssign     (const AssignExpr* expr);

		ExprPtr applyRules(const Expression* expr);

		using RuleFn = std::function<ExprPtr (const Expression* e)>;
		std::vector<RuleFn> rules;
	};
}
