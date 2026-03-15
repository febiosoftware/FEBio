#include "simplifier.h"
using namespace febcode;

ExprPtr febcode::simplify(const Expression* expr)
{
	if (auto binary = dynamic_cast<const BinaryExpr*>(expr))
	{
		auto l = simplify(binary->left.get());
		auto r = simplify(binary->right.get());

		// If both sides are literals, we can compute the result at compile time
		if (isLiteral(l) && isLiteral(r))
		{
			Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
			Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;

			switch (binary->op)
			{
			case BinaryOp::Plus    : return Literal(lv + rv);
			case BinaryOp::Minus   : return Literal(lv - rv);
			case BinaryOp::Multiply: return Literal(lv * rv);
			case BinaryOp::Divide  : return Literal(lv / rv);
			case BinaryOp::Exponent:
			{
				if (isDouble(lv) && isDouble(rv))
				{
					double a = getDouble(lv);
					double b = getDouble(rv);
					return Literal(std::pow(a, b));
				}
				else
					throw std::runtime_error("Don't know how to exponeniate expression.");
			}
			default:
				throw std::runtime_error("Unsupported binary operator for simplification");
			}
		}

		// if only left operand is a literal and it's zero or one
		if (isLiteral(l))
		{
			Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
			if (isZero(lv))
			{
				switch (binary->op)
				{
				case BinaryOp::Plus    : return  r; // 0 + r = r
				case BinaryOp::Minus   : return -r; // 0 - r = -r
				case BinaryOp::Multiply: return Literal(0.0); // 0 * r = 0
				case BinaryOp::Divide  : return Literal(0.0); // 0 / r = 0
				}
			}

			if (isOne(lv))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return r; // 1 * r = r
				case BinaryOp::Exponent: return Literal(1.0); // 1 ^ r = 1
				}
			}

			if ((binary->op == BinaryOp::Multiply) && isBinary(r))
			{
				auto rbin = dynamic_cast<BinaryExpr*>(r.get());
				if (rbin->op == BinaryOp::Multiply)
				{
					auto a = simplify(rbin->left);
					auto b = simplify(rbin->right);

					if (isLiteral(a))
					{
						Value av = dynamic_cast<LiteralExpr*>(a.get())->value;
						return Binary(Literal(av * lv), BinaryOp::Multiply, b); // l * (a * b) --> (l * a) * b
					}
				}
			}
		}

		// if only right operand is a literal and it's zero or one
		if (isLiteral(r))
		{
			Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;
			if (isZero(rv))
			{
				switch (binary->op)
				{
				case BinaryOp::Plus    : return l; // l + 0 = l
				case BinaryOp::Minus   : return l; // l - 0 = l
				case BinaryOp::Multiply: return Literal(0.0); // l * 0 = 0
				case BinaryOp::Divide  : // l / 0 --> no can do!
					throw std::runtime_error("Division by zero in evaluation of derivative!");
				case BinaryOp::Exponent: return Literal(1.0); // l ^ 0 = 1
				}
			}

			if (isOne(rv))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return l; // l * 1 = l
				case BinaryOp::Divide  : return l; // l / 1 = l
				case BinaryOp::Exponent: return l; // l ^ 1 = l
				}
			}

			if (isInt(rv) && (getInt(rv) < 0))
			{
				if (binary->op == BinaryOp::Plus)
				{
					return simplify(Binary(l, BinaryOp::Minus, Literal(-getInt(rv))));
				}
				if (binary->op == BinaryOp::Minus)
				{
					return simplify(Binary(l, BinaryOp::Plus, Literal(-getInt(rv))));
				}
			}

			if (isDouble(rv) && (getDouble(rv) < 0))
			{
				if (binary->op == BinaryOp::Plus)
				{
					return simplify(Binary(l, BinaryOp::Minus, Literal(-getDouble(rv))));
				}
				if (binary->op == BinaryOp::Minus)
				{
					return simplify(Binary(l, BinaryOp::Plus, Literal(-getDouble(rv))));
				}
			}

			if (binary->op == BinaryOp::Multiply)
			{
				// move literal to the left for better chances of simplification
				return simplify(Binary(r, BinaryOp::Multiply, l));
			}
		}

		// if left and right are the same, we can simplify
		if (febcode::isEqual(l, r))
		{
			if (binary->op == BinaryOp::Plus    ) return simplify(Binary(Literal(2.0), BinaryOp::Multiply, l));
			if (binary->op == BinaryOp::Minus   ) return Literal(0.0);
			if (binary->op == BinaryOp::Multiply) return Binary(l, BinaryOp::Exponent, Literal(2.0));
		}

		if (isNegation(r))
		{
			if (binary->op == BinaryOp::Plus)
			{
				// a + (-b) --> a - b
				auto neg = dynamic_cast<UnaryExpr*>(r.get());
				return simplify(Binary(l, BinaryOp::Minus, neg->right));
			}
			if (binary->op == BinaryOp::Minus)
			{
				// a - (-b) --> a + b
				auto neg = dynamic_cast<UnaryExpr*>(r.get());
				return simplify(Binary(l, BinaryOp::Plus, neg->right));
			}
		}

		return Binary(l, binary->op, r);
	}
	if (auto unary = dynamic_cast<const UnaryExpr*>(expr))
	{
		if (unary->op == UnaryOp::Negate)
		{
			// absorb negative signs in number
			if (auto lit = dynamic_cast<LiteralExpr*>(unary->right.get()))
			{
				const Value& v = lit->value;
				if (isInt(v)) return Literal(-getInt(v));
				if (isDouble(v)) return Literal(-getDouble(v));
			}
			else if (auto innerUnary = dynamic_cast<const UnaryExpr*>(unary->right.get()))
			{
				if (innerUnary->op == UnaryOp::Negate)
				{
					// --x --> x
					return simplify(innerUnary->right.get());
				}
			}

			auto simplifiedRight = simplify(unary->right.get());
			if (isEqual(unary->right, simplifiedRight))
			{
				// no change after simplification, return original
				return copy_expression(expr);
			}
			else
				return simplify(Unary(UnaryOp::Negate, simplifiedRight));
		}
	}
	if (auto init = dynamic_cast<const InitializerExpr*>(expr))
	{
		std::vector<ExprPtr> simplifiedElements;
		for (const auto& elem : init->elements)
		{
			simplifiedElements.push_back(simplify(elem.get()));
		}
		return std::make_unique<InitializerExpr>(std::move(simplifiedElements));
	}

	return copy_expression(expr);
}
