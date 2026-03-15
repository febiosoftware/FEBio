#pragma once
#include "ast.h"

namespace febcode {
	ExprPtr simplify(const Expression* expr);
	inline ExprPtr simplify(const std::unique_ptr<Expression>& expr) { return simplify(expr.get()); }
}
