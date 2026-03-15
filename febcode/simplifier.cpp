#include "simplifier.h"
using namespace febcode;

ExprPtr febcode::simplify(const Expression* expr)
{
	if (auto binary = dynamic_cast<const BinaryExpr*>(expr))
	{
		auto l = simplify(binary->left.get());
		auto r = simplify(binary->right.get());

		// If both sides are literals, we can compute the result at compile time
		if (dynamic_cast<LiteralExpr*>(l.get()) && dynamic_cast<LiteralExpr*>(r.get()))
		{
			Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
			Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;

			switch (binary->op)
			{
			case BinaryOp::Plus: return Literal(lv + rv);
			case BinaryOp::Minus: return Literal(lv - rv);
			case BinaryOp::Multiply: return Literal(lv * rv);
			case BinaryOp::Divide: return Literal(lv / rv);
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
		if (dynamic_cast<LiteralExpr*>(l.get()))
		{
			Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
			if (isZero(lv))
			{
				switch (binary->op)
				{
				case BinaryOp::Plus: return r; // 0 + r = r
				case BinaryOp::Minus: return -r; // 0 - r = -r
				case BinaryOp::Multiply: return Literal(0.0); // 0 * r = 0
				case BinaryOp::Divide: return Literal(0.0); // 0 / r = 0
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
		}

		// if only right operand is a literal and it's zero or one
		if (dynamic_cast<LiteralExpr*>(r.get()))
		{
			Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;
			if (isZero(rv))
			{
				switch (binary->op)
				{
				case BinaryOp::Plus: return l; // l + 0 = l
				case BinaryOp::Minus: return l; // l - 0 = l
				case BinaryOp::Multiply: return Literal(0.0); // l * 0 = 0
				case BinaryOp::Divide: // l / 0 --> no can do!
					throw std::runtime_error("Division by zero in evaluation of derivative!");
				case BinaryOp::Exponent: return Literal(1.0); // l ^ 0 = 1
				}
			}

			if (isOne(rv))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return l; // l * 1 = l
				case BinaryOp::Divide: return l; // l / 1 = l
				case BinaryOp::Exponent: return l; // l ^ 1 = l
				}
			}
		}

		return Binary(l, binary->op, r);
	}
	if (auto unary = dynamic_cast<const UnaryExpr*>(expr))
	{
		if (unary->op == UnaryOp::Negate)
		{
			if (auto lit = dynamic_cast<LiteralExpr*>(unary->right.get()))
			{
				const Value& v = lit->value;
				if (isInt(v)) return Literal(-getInt(v));
				if (isDouble(v)) return Literal(-getDouble(v));
			}
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
