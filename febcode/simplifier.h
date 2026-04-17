#pragma once
#include "ast.h"
#include "modifier.h"

namespace febcode {

	class Simplifier : public Modifier {
	public:
		Simplifier(Program& prg) : Modifier(prg) {}

		ExprPtr simplify(const Expression* expr);
		inline ExprPtr simplify(const std::unique_ptr<Expression>& expr) { return simplify(expr.get()); }

	private:
		ExprPtr simplifyConstructor(const ConstructorExpr* expr);
		ExprPtr simplifyUnary      (const UnaryExpr* expr);
		ExprPtr simplifyInitializer(const InitExpr* expr);
		ExprPtr simplifyCall       (const CallExpr* expr);
		ExprPtr simplifyBinary     (const BinaryExpr* expr);
	};
}
