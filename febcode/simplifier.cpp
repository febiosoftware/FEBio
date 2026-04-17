#include "simplifier.h"
using namespace febcode;

ExprPtr Simplifier::simplify(const Expression* expr)
{
	if (auto ctor   = dynamic_cast<const ConstructorExpr*>(expr)) return simplifyConstructor(ctor  );
	if (auto unary  = dynamic_cast<const UnaryExpr*      >(expr)) return simplifyUnary      (unary );
	if (auto init   = dynamic_cast<const InitExpr*       >(expr)) return simplifyInitializer(init  );
	if (auto call   = dynamic_cast<const CallExpr*       >(expr)) return simplifyCall       (call  );
	if (auto binary = dynamic_cast<const BinaryExpr*     >(expr)) return simplifyBinary     (binary);
	return copy_expression(expr);
}

ExprPtr Simplifier::simplifyConstructor(const ConstructorExpr* ctor)
{
	// simplify all arguments first
	std::vector<ExprPtr> simplifiedArgs;
	for (const auto& arg : ctor->args)
	{
		simplifiedArgs.push_back(simplify(arg.get()));
	}

	// if all args are numbers, we'll make this a literal too
	if (std::all_of(simplifiedArgs.begin(), simplifiedArgs.end(), [](const ExprPtr& arg) {
		auto lit = dynamic_cast<const LiteralExpr*>(arg.get());
		return lit && (isDouble(lit->value) || isInt(lit->value));
		}))
	{
		std::vector<double> v;
		for (const auto& arg : simplifiedArgs)
		{
			auto lit = dynamic_cast<const LiteralExpr*>(arg.get());
			if (isInt(lit->value))
				v.push_back(lit->value.i);
			else
				v.push_back(lit->value.d);
		}

		TypeKind type = ctor->valType->kind;
		if ((type == TypeKind::Vec2) && (v.size() == 1)) return Literal(vec2(v[0], v[0]));
		if ((type == TypeKind::Vec2) && (v.size() == 2)) return Literal(vec2(v[0], v[1]));
		if ((type == TypeKind::Vec3) && (v.size() == 1)) return Literal(vec3(v[0], v[0], v[0]));
		if ((type == TypeKind::Vec3) && (v.size() == 3)) return Literal(vec3(v[0], v[1], v[2]));
		if ((type == TypeKind::Mat2) && (v.size() == 1)) return Literal(mat2(v[0], 0.0, 0.0, v[0]));
		if ((type == TypeKind::Mat2) && (v.size() == 4)) return Literal(mat2(v[0], v[1], v[2], v[3]));
		if ((type == TypeKind::Mat3) && (v.size() == 1)) return Literal(mat3(v[0], 0.0, 0.0, 0.0, v[0], 0.0, 0.0, 0.0, v[0]));
		if ((type == TypeKind::Mat3) && (v.size() == 9)) return Literal(mat3(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]));
	}

	return std::make_unique<ConstructorExpr>(ctor->valType, std::move(simplifiedArgs));
}

ExprPtr Simplifier::simplifyUnary(const UnaryExpr* unary)
{
	if (unary->op == UnaryOp::Negate)
	{
		auto simplifiedRight = simplify(unary->right.get());

		// -0 = 0
		if (isZero(simplifiedRight)) {
			return simplifiedRight;
		}

		// absorb negative signs in number
		if (auto lit = dynamic_cast<LiteralExpr*>(simplifiedRight.get()))
		{
			const Value& v = lit->value;
			if (isInt(v)) return Literal(-getInt(v));
			if (isDouble(v)) return Literal(-getDouble(v));
		}
		else if (auto innerUnary = dynamic_cast<const UnaryExpr*>(simplifiedRight.get()))
		{
			if (innerUnary->op == UnaryOp::Negate)
			{
				// --x --> x
				return simplify(innerUnary->right.get());
			}
		}

		return Negate(simplifiedRight);
	}

	return copy_expression(unary);
}

ExprPtr Simplifier::simplifyInitializer(const InitExpr* init)
{
	std::vector<ExprPtr> simplifiedElements;
	for (const auto& elem : init->elements)
	{
		simplifiedElements.push_back(simplify(elem.get()));
	}
	return std::make_unique<InitExpr>(std::move(simplifiedElements));
}

ExprPtr Simplifier::simplifyCall(const CallExpr* call)
{
	std::vector<ExprPtr> simplifiedArgs;
	for (const auto& arg : call->arguments)
	{
		simplifiedArgs.push_back(simplify(arg.get()));
	}

	// handle some special cases for built-in functions
	std::string fncName = call->name;
	if ((fncName == "outer") && (simplifiedArgs.size() == 2))
	{
		auto& le = simplifiedArgs[0];
		auto& re = simplifiedArgs[1];

		if (isLiteral(le) && isZero(le))
		{
			Value lv = dynamic_cast<LiteralExpr*>(le.get())->value;
			if (isVec2(lv)) return Literal(mat2(0.0));
			else if (isVec3(lv)) return Literal(mat3(0.0));
		}

		if (isLiteral(re) && isZero(re))
		{
			Value rv = dynamic_cast<LiteralExpr*>(re.get())->value;
			if (isVec2(rv)) return Literal(mat2(0.0));
			else if (isVec3(rv)) return Literal(mat3(0.0));
		}
	}
	if ((fncName == "cross") && (simplifiedArgs.size() == 2))
	{
		auto& le = simplifiedArgs[0];
		auto& re = simplifiedArgs[1];

		if (isLiteral(le) && isZero(le))
		{
			return Literal(vec3(0.0, 0.0, 0.0));
		}

		if (isLiteral(re) && isZero(re))
		{
			return Literal(vec3(0.0, 0.0, 0.0));
		}
	}
	if ((fncName == "dot") && (simplifiedArgs.size() == 2))
	{
		auto& le = simplifiedArgs[0];
		auto& re = simplifiedArgs[1];

		if (isZero(le) || isZero(re))
		{
			return Literal(0.0);
		}
	}
	if ((fncName == "transpose") && (simplifiedArgs.size() == 1))
	{
		auto& e = simplifiedArgs[0];
		if (isLiteral(e))
		{
			Value v = dynamic_cast<LiteralExpr*>(e.get())->value;
			if (isMat2(v))
			{
				if (isSymmetric(v.mat2Value)) return Literal(v.mat2Value);
			}
			else if (isMat3(v))
			{
				if (isSymmetric(v.mat3Value)) return Literal(v.mat3Value);
			}
		}
	}

	return Call(fncName, std::move(simplifiedArgs));
}

ExprPtr Simplifier::simplifyBinary(const BinaryExpr* binary)
{
	auto l = simplify(binary->left.get());
	auto r = simplify(binary->right.get());

	// for scalar operations, we can compute the result here
	if (isScalar(l) && isScalar(r))
	{
		Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
		Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;

		if (isInt(lv) && isInt(rv))
		{
			int a = getInt(lv);
			int b = getInt(rv);

			switch (binary->op)
			{
			case BinaryOp::Plus    : return Literal(a + b);
			case BinaryOp::Minus   : return Literal(a - b);
			case BinaryOp::Multiply: return Literal(a * b);
			case BinaryOp::Divide  : return Literal(a / b);
			case BinaryOp::Exponent: return Literal(static_cast<int>(std::pow(a, b)));
			}
		}
		else
		{
			double a = (isInt(lv) ? static_cast<double>(getInt(lv)) : getDouble(lv));
			double b = (isInt(rv) ? static_cast<double>(getInt(rv)) : getDouble(rv));

			switch (binary->op)
			{
			case BinaryOp::Plus    : return Literal(a + b);
			case BinaryOp::Minus   : return Literal(a - b);
			case BinaryOp::Multiply: return Literal(a * b);
			case BinaryOp::Divide  : return Literal(a / b);
			case BinaryOp::Exponent: return Literal(std::pow(a, b));
			}
		}
	}

	// plus/minus apply to all types, so we can simplify them even if we don't know the exact type of the operands
	if (binary->op == BinaryOp::Plus)
	{
		if (isZero(l)) return r; // 0 + r = r
		if (isZero(r)) return l; // l + 0 = l
	}

	if (binary->op == BinaryOp::Minus)
	{
		if (isZero(l)) return Negate(r); // 0 - r = -r
		if (isZero(r)) return  l; // l - 0 = l
	}

	// if only left operand is a scalar and it's zero or one
	if (isLiteral(l))
	{
		Value lv = dynamic_cast<LiteralExpr*>(l.get())->value;
		if (isInt(lv) || isDouble(lv))
		{
			if (isZero(lv))
			{
				switch (binary->op)
				{
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

			if ((binary->op == BinaryOp::Multiply) && isBinary(r))
			{
				auto rbin = dynamic_cast<BinaryExpr*>(r.get());
				if (rbin->op == BinaryOp::Multiply)
				{
					auto a = simplify(rbin->left);
					auto b = simplify(rbin->right);

					if (isScalar(a))
					{
						ExprPtr newLeft = simplify(Binary(BinaryOp::Multiply, l, a)); // l * a
						return simplify(Binary(BinaryOp::Multiply, newLeft, b)); // l * (a * b) --> (l * a) * b
					}
				}
			}
		}

		if (isVec2(lv))
		{
			vec2 v = getVec2(lv);
			if (isZero(v) && r->valType)
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply:
				{
					if (isScalarType(r->valType))
						return Literal(vec2(0.0, 0.0)); // vec2(0) * scalar = vec2(0)
					if (r->valType->kind == TypeKind::Vec2)
						return Literal(0.0); // vec2(0) * vec2 = 0
				}
				break;
				}
			}
		}

		if (isVec3(lv))
		{
			vec3 v = getVec3(lv);
			if (isZero(v))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply:
					if (isScalarType(r->valType))
						return Literal(vec3(0.0, 0.0, 0.0)); // vec3(0) * scalar = vec3(0)
					if (r->valType->kind == TypeKind::Vec3)
						return Literal(0.0); // vec3(0) * vec3 = 0
				}
			}
		}

		if (isMat2(lv))
		{
			mat2 m = getMat2(lv);
			if (isIdentity(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return r; // I * r = r
				}
			}
		}

		if (isMat3(lv))
		{
			mat3 m = getMat3(lv);
			if (isIdentity(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return r; // I * r = r
				}
			}
		}
	}

	// if only right operand is a scalar and it's zero or one
	if (isLiteral(r))
	{
		Value rv = dynamic_cast<LiteralExpr*>(r.get())->value;
		if (isInt(rv) || isDouble(rv))
		{
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

			if (isInt(rv) && (getInt(rv) < 0))
			{
				if (binary->op == BinaryOp::Plus)
				{
					return simplify(Binary(BinaryOp::Minus, l, Literal(-getInt(rv))));
				}
				if (binary->op == BinaryOp::Minus)
				{
					return simplify(Binary(BinaryOp::Plus, l, Literal(-getInt(rv))));
				}
			}

			if (isDouble(rv) && (getDouble(rv) < 0))
			{
				if (binary->op == BinaryOp::Plus)
				{
					return simplify(Binary(BinaryOp::Minus, l, Literal(-getDouble(rv))));
				}
				if (binary->op == BinaryOp::Minus)
				{
					return simplify(Binary(BinaryOp::Plus, l, Literal(-getDouble(rv))));
				}
			}

			if (binary->op == BinaryOp::Multiply)
			{
				// move literal to the left for better chances of simplification
				return simplify(Binary(BinaryOp::Multiply, r, l));
			}
		}

		if (isVec3(rv))
		{
			vec3 v = getVec3(rv);
			if (isZero(v))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply:
					if (isScalarType(l->valType))
						return Literal(vec3(0.0, 0.0, 0.0)); // scalar * vec3(0) = vec3(0)
					if (l->valType->kind == TypeKind::Vec3)
						return Literal(0.0); // vec3 * vec3(0) = 0
				}
			}
		}

		if (isMat2(rv))
		{
			mat2 m = getMat2(rv);
			if (isZero(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: // l*0 = 0
				{
					// we need the signature for this operation to determine the result type
					BinaryOpSignature sig = prg.resolveBinaryOp(binary->op, l->valType, r->valType);
					return Zero(sig.resultType);
				}
				}
			}
			if (isIdentity(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return l; // l * I = l
				}
			}
		}

		if (isMat3(rv))
		{
			mat3 m = getMat3(rv);
			if (isZero(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return Zero(r->valType); // l * mat3(0) = mat3(0)
				}
			}
			if (isIdentity(m))
			{
				switch (binary->op)
				{
				case BinaryOp::Multiply: return l; // l * I = l
				}
			}
		}
	}

	// if left and right are the same, we can simplify
	if (febcode::isEqual(l, r))
	{
		if (binary->op == BinaryOp::Plus    ) return simplify(Binary(BinaryOp::Multiply, Literal(2.0), l));
		if (binary->op == BinaryOp::Minus   ) return Zero(l->valType);
		if ((binary->op == BinaryOp::Multiply) && isScalarType(l->valType)) return Binary(BinaryOp::Exponent, l, Literal(2.0));
	}

	if (isNegation(r))
	{
		if (binary->op == BinaryOp::Plus)
		{
			// a + (-b) --> a - b
			auto neg = dynamic_cast<UnaryExpr*>(r.get());
			return simplify(Binary(BinaryOp::Minus, l, neg->right));
		}
		if (binary->op == BinaryOp::Minus)
		{
			// a - (-b) --> a + b
			auto neg = dynamic_cast<UnaryExpr*>(r.get());
			return simplify(Binary(BinaryOp::Plus, l, neg->right));
		}
	}

	return Binary(binary->op, l, r);
}
