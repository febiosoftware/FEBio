#include "ast.h"
#include <iostream>

using namespace febcode;

ExprPtr febcode::Initializer(const std::vector<StructField>& fields)
{
	std::vector<ExprPtr> init;
	for (const auto& field : fields)
	{
		switch (field.first->kind)
		{
		case TypeKind::Bool  : init.push_back(Literal(Value(false))); break;
		case TypeKind::Int   : init.push_back(Literal(Value(0))); break;
		case TypeKind::Double: init.push_back(Literal(Value(0.0))); break;
		case TypeKind::Vec2  : init.push_back(Literal(Value(vec2(0., 0.)))); break;
		case TypeKind::Vec3  : init.push_back(Literal(Value(vec3(0., 0., 0.)))); break;
		case TypeKind::Array:
		{
			Value v;
			switch (field.first->elementType->kind)
			{
			case TypeKind::Bool: v = Value(false); break;
			case TypeKind::Int: v = Value(0); break;
			case TypeKind::Double: v = Value(0.0); break;
			case TypeKind::Vec2: v = Value(vec2(0., 0.)); break;
			case TypeKind::Vec3: v = Value(vec3(0., 0., 0.)); break;
			default:
				throw std::runtime_error("Unsupported array element type for initializer");
			};
			init.push_back(Initializer(field.first->arraySize, v));
			break;
		}
		case TypeKind::Struct:
		{
			std::vector<ExprPtr> structFields;
			for (const auto& subfield : field.first->fields)
			{
				Value v;
				switch (subfield.first->elementType->kind)
				{
				case TypeKind::Bool  : v = Value(false); break;
				case TypeKind::Int   : v = Value(0); break;
				case TypeKind::Double: v = Value(0.0); break;
				case TypeKind::Vec2  : v = Value(vec2(0., 0.)); break;
				case TypeKind::Vec3  : v = Value(vec3(0., 0., 0.)); break;
				default:
					throw std::runtime_error("Unsupported struct field type for initializer");
				};

				structFields.push_back(Literal(v));
			}
			init.push_back(std::make_unique<InitializerExpr>(std::move(structFields)));
		}
		break;
		default:
			throw std::runtime_error("Unsupported struct field type for initializer");
		}
	}
	return std::make_unique<InitializerExpr>(std::move(init));
}

ExprPtr febcode::copy_expression(const Expression* expr)
{
	if (expr == nullptr)  return nullptr;
	else if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
		return Literal(literal->value);
	}
	else if (auto variable = dynamic_cast<const VariableExpr*>(expr)) {
		return Variable(variable->name);
	}
	else if (auto unary = dynamic_cast<const UnaryExpr*>(expr)) {
		return Unary(unary->op, unary->right);
	}
	else if (auto binary = dynamic_cast<const BinaryExpr*>(expr)) {
		return Binary(binary->left, binary->op, binary->right);
	}
	else if (auto call = dynamic_cast<const CallExpr*>(expr))
	{
		return Call(call->callee, call->arguments);
	}
	else if (auto call = dynamic_cast<const InitializerExpr*>(expr))
	{
		std::vector<ExprPtr> copyArgs;
		for (auto& arg : call->elements)
		{
			copyArgs.emplace_back(copy_expression(arg.get()));
		}
		return std::make_unique<InitializerExpr>(std::move(copyArgs));
	}
	else if (auto assign = dynamic_cast<const AssignExpr*>(expr))
	{
		return Assign(assign->target, assign->value);
	}
	else if (auto member = dynamic_cast<const MemberExpr*>(expr))
	{
		return Member(member->object, member->property);
	}
	else
	{
		throw std::runtime_error("Unsupported expression type for copying");
	}
}

bool febcode::isEqual(const ExprPtr& l, const ExprPtr& r)
{
	if (!l && !r) return true;
	if (!l || !r) return false;

	if (auto litL = dynamic_cast<LiteralExpr*>(l.get()))
	{
		if (auto litR = dynamic_cast<LiteralExpr*>(r.get()))
		{
			return litL->value == litR->value;
		}
	}
	else if (auto varL = dynamic_cast<VariableExpr*>(l.get()))
	{
		if (auto varR = dynamic_cast<VariableExpr*>(r.get()))
		{
			return varL->name == varR->name;
		}
	}
	else if (auto binL = dynamic_cast<BinaryExpr*>(l.get()))
	{
		if (auto binR = dynamic_cast<BinaryExpr*>(r.get()))
		{
			return binL->op == binR->op &&
				isEqual(binL->left, binR->left) &&
				isEqual(binL->right, binR->right);
		}
	}
	else if (auto unL = dynamic_cast<UnaryExpr*>(l.get()))
	{
		if (auto unR = dynamic_cast<UnaryExpr*>(r.get()))
		{
			return unL->op == unR->op &&
				isEqual(unL->right, unR->right);
		}
	}
	else if (auto memberL = dynamic_cast<MemberExpr*>(l.get()))
	{
		if (auto memberR = dynamic_cast<MemberExpr*>(r.get()))
		{
			return isEqual(memberL->object, memberR->object) &&
				memberL->property == memberR->property;
		}
	}
	else if (auto indexL = dynamic_cast<IndexExpr*>(l.get()))
	{
		if (auto indexR = dynamic_cast<IndexExpr*>(r.get()))
		{
			return isEqual(indexL->object, indexR->object) &&
				isEqual(indexL->index, indexR->index);
		}
	}
	else if (auto callL = dynamic_cast<CallExpr*>(l.get()))
	{
		if (auto callR = dynamic_cast<CallExpr*>(r.get()))
		{
			return isEqual(callL->callee, callR->callee) && isEqual(callL->arguments, callR->arguments);
		}
	}

	return false;
}

bool febcode::isEqual(const std::vector<ExprPtr>& l, const std::vector<ExprPtr>& r)
{
	if (l.size() != r.size()) return false;
	for (int i = 0; i < l.size(); ++i)
	{
		if (!isEqual(l[i], r[i])) return false;
	}
	return true;
}

std::ostream& operator << (std::ostream& o, const febcode::Value& v)
{
	if      (isVoid  (v)) return o << "null";
	else if (isInt   (v)) return o << getInt(v);
	else if (isDouble(v)) return o << getDouble(v);
	else if (isBool  (v)) return o << (getBool(v) ? "true" : "false");
	else if (isArray(v))
	{
		const febcode::ArrayValue& arr = getArray(v);
		o << "[";
		for (size_t i = 0; i < arr.size(); ++i)
		{
			o << arr.elements[i];
			if (i != arr.size() - 1) o << ", ";
		}
		return o << "]";
	}
	else if (isStruct(v))
	{
		const febcode::StructValue& s = getStruct(v);
		o << "{";
		for (size_t i = 0; i < s.fields.size(); ++i)
		{
			o << s.fields[i];
			if (i != s.fields.size() - 1) o << ", ";
		}
		return o << "}";
	}
	else if (isVec2(v))
	{
		const febcode::vec2& vec = getVec2(v);
		return o << "vec2(" << vec.x << ", " << vec.y << ")";
	}
	else if (isVec3(v))
	{
		const febcode::vec3& vec = getVec3(v);
		return o << "vec3(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
	}
	else if (isMat2(v))
	{
		const febcode::mat2& mat = getMat2(v);
		return o << "mat2(" << mat.m[0][0] << ", " << mat.m[0][1] << ", " << mat.m[1][0] << ", " << mat.m[1][1] << ")";
	}
	else if (isMat3(v))
	{
		const febcode::mat3& mat = getMat3(v);
		return o << "mat3(" << mat.m[0][0] << ", " << mat.m[0][1] << ", " << mat.m[0][2] << ", " << mat.m[1][0] << ", " << mat.m[1][1] << ", " << mat.m[1][2] << ", " << mat.m[2][0] << ", " << mat.m[2][1] << ", " << mat.m[2][2] << ")";
	}
	else
		return o << "<unknown value>";
}

std::string ValueToString(const febcode::Value& v)
{
	std::string s;
	if      (isVoid  (v)) s = "null";
	else if (isBool  (v)) s = getBool(v) ? "true" : "false";
	else if (isInt   (v)) s = std::to_string(getInt (v));
	else if (isDouble(v)) s = std::to_string(getDouble(v));
	else if (isVec2  (v)) s = "vec2";
	else if (isVec3  (v)) s = "vec3";
	else if (isArray (v)) { auto& p = getArrayPtr (v); s = "[array:"  + std::to_string(p.use_count()) + "]"; }
	else if (isStruct(v)) { 
		const febcode::StructValue& o = febcode::getStruct(v);
		const std::string& name = o.type->name;
		auto& p = getStructPtr(v);
		s = "[" + name;
		s = "{";
		for (size_t i = 0; i < o.fields.size(); ++i)
		{
			s += ValueToString(o.fields[i]);
			if (i != o.fields.size() - 1) s += ", ";
		}
		s += "}";
		s += ":" + std::to_string(p.use_count()) + "]";
	}
	else s = "unknown";
	return s;
}

std::string ValueTypeToString(const febcode::Value& v)
{
	if      (isVoid  (v)) return "void";
	else if (isBool  (v)) return "bool";
	else if (isInt   (v)) return "int";
	else if (isDouble(v)) return "double";
	else if (isArray (v)) return "array";
	else if (isStruct(v)) return "struct";
	else if (isVec2  (v)) return "vec2";
	else if (isVec3  (v)) return "vec3";
	else if (isRef    (v)) return "ref";
	else return "<unknown type>";
}
