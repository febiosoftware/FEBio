#include "ast.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace febcode;

ExprPtr febcode::clone(const Expression* expr)
{
	ExprPtr cpy;

	if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
		cpy = std::make_unique<LiteralExpr>(literal->value);
	}
	else if (auto variable = dynamic_cast<const VariableExpr*>(expr)) {
		cpy = std::make_unique<VariableExpr>(variable->name);
	}
	else if (auto unary = dynamic_cast<const UnaryExpr*>(expr)) {
		cpy = std::make_unique<UnaryExpr>(unary->op, clone(unary->right.get()));
	}
	else if (auto binary = dynamic_cast<const BinaryExpr*>(expr)) {
		cpy = std::make_unique<BinaryExpr>(clone(binary->left.get()), binary->op, clone(binary->right.get()));
	}
	else if (auto call = dynamic_cast<const CallExpr*>(expr))
	{
		std::vector<ExprPtr> copyArgs;
		for (auto& arg : call->arguments)
		{
			copyArgs.emplace_back(clone(arg.get()));
		}
		cpy = std::make_unique<CallExpr>(call->name, std::move(copyArgs));
	}
	else if (auto call = dynamic_cast<const InitExpr*>(expr))
	{
		std::vector<ExprPtr> copyArgs;
		for (auto& arg : call->elements)
		{
			copyArgs.emplace_back(clone(arg.get()));
		}
		cpy = std::make_unique<InitExpr>(std::move(copyArgs));
	}
	else if (auto constructor = dynamic_cast<const ConstructorExpr*>(expr))
	{
		std::vector<ExprPtr> copyArgs;
		for (auto& arg : constructor->args)
		{
			copyArgs.emplace_back(clone(arg.get()));
		}
		cpy = std::make_unique<ConstructorExpr>(constructor->valType, std::move(copyArgs));
	}
	else if (auto assign = dynamic_cast<const AssignExpr*>(expr))
	{
		cpy = std::make_unique<AssignExpr>(clone(assign->target.get()), clone(assign->value.get()));
	}
	else if (auto member = dynamic_cast<const MemberExpr*>(expr))
	{
		cpy = std::make_unique<MemberExpr>(clone(member->object.get()), member->property);
	}
	else if (auto index = dynamic_cast<const IndexExpr*>(expr))
	{
		cpy = std::make_unique<IndexExpr>(clone(index->object.get()), clone(index->index.get()));
	}
	else
	{
		throw std::runtime_error("Unsupported expression type for copying");
	}

	if (cpy)
		cpy->valType = expr->valType; // Copy the value type as well

	return cpy;
}

bool febcode::isEqual(const ExprPtr& l, const ExprPtr& r)
{
	return isEqual(l.get(), r.get());
}

bool febcode::isEqual(const Expression* l, const Expression* r)
{
	if (!l && !r) return true;
	if (!l || !r) return false;

	if (auto litL = dynamic_cast<const LiteralExpr*>(l))
	{
		if (auto litR = dynamic_cast<const LiteralExpr*>(r))
		{
			return litL->value == litR->value;
		}
	}
	else if (auto varL = dynamic_cast<const VariableExpr*>(l))
	{
		if (auto varR = dynamic_cast<const VariableExpr*>(r))
		{
			return varL->name == varR->name;
		}
	}
	else if (auto binL = dynamic_cast<const BinaryExpr*>(l))
	{
		if (auto binR = dynamic_cast<const BinaryExpr*>(r))
		{
			return binL->op == binR->op &&
				isEqual(binL->left, binR->left) &&
				isEqual(binL->right, binR->right);
		}
	}
	else if (auto unL = dynamic_cast<const UnaryExpr*>(l))
	{
		if (auto unR = dynamic_cast<const UnaryExpr*>(r))
		{
			return unL->op == unR->op &&
				isEqual(unL->right, unR->right);
		}
	}
	else if (auto memberL = dynamic_cast<const MemberExpr*>(l))
	{
		if (auto memberR = dynamic_cast<const MemberExpr*>(r))
		{
			return isEqual(memberL->object, memberR->object) &&
				memberL->property == memberR->property;
		}
	}
	else if (auto indexL = dynamic_cast<const IndexExpr*>(l))
	{
		if (auto indexR = dynamic_cast<const IndexExpr*>(r))
		{
			return isEqual(indexL->object, indexR->object) &&
				isEqual(indexL->index, indexR->index);
		}
	}
	else if (auto callL = dynamic_cast<const CallExpr*>(l))
	{
		if (auto callR = dynamic_cast<const CallExpr*>(r))
		{
			return (callL->name == callR->name) && isEqual(callL->arguments, callR->arguments);
		}
	}
	else if (auto ctorL = dynamic_cast<const ConstructorExpr*>(l))
	{
		if (auto ctorR = dynamic_cast<const ConstructorExpr*>(r))
		{
			if (ctorL->valType != ctorR->valType) return false;
			assert(ctorL->args.size() == ctorR->args.size());
			if (ctorL->args.size() != ctorR->args.size()) return false;
			for (int i = 0; i < ctorL->args.size(); ++i)
			{
				if (!isEqual(ctorL->args[i], ctorR->args[i]))
					return false;
			}
			return true;
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
		return o << "mat2(" << mat(0,0) << ", " << mat(0,1) << ", " << mat(1,0) << ", " << mat(1,1) << ")";
	}
	else if (isMat3(v))
	{
		const febcode::mat3& mat = getMat3(v);
		return o << "mat3(" << mat(0,0) << ", " << mat(0,1) << ", " << mat(0,2) << ", " << mat(1,0) << ", " << mat(1,1) << ", " << mat(1,2) << ", " << mat(2,0) << ", " << mat(2,1) << ", " << mat(2,2) << ")";
	}
	else
		return o << "<unknown value>";
}

std::string to_nice_string(double d)
{
	std::ostringstream ss;
	ss << d;
	return ss.str();
}

std::string febcode::ValueToString(const febcode::Value& v)
{
	std::string s;
	if      (isVoid  (v)) s = "null";
	else if (isBool  (v)) s = getBool(v) ? "true" : "false";
	else if (isInt   (v)) s = std::to_string(getInt (v));
	else if (isDouble(v)) s = to_nice_string(getDouble(v));
	else if (isVec2(v))
	{
		s = "vec2(";
		vec2 q = getVec2(v);
		if (isZero(q))
			s += to_nice_string(0.0);
		else
			s += to_nice_string(q.x) + ", " + to_nice_string(q.y);
		s += ")";
	}
	else if (isVec3(v))
	{
		s = "vec3(";
		vec3 q = getVec3(v);
		if (isZero(q))
		{
			s += to_nice_string(0.0);
		}
		else
		{
			s += to_nice_string(q.x) + ", " + to_nice_string(q.y) + ", " + to_nice_string(q.z);
		}
		s += ")";
	}
	else if (isMat2(v))
	{
		s = "mat2(";

		const mat2& m = v.mat2Value;
		if (isZero(m))
		{
			s += to_nice_string(0.0);
		}
		else if (isIdentity(m))
		{
			s += to_nice_string(1.0);
		}
		else
		{
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j)
				{
					s += to_nice_string(m(i,j));
					if ((i != 1) || (j != 1)) s += ",";
				}
		}
		s += ")";
	}
	else if (isMat3(v))
	{
		s = "mat3(";

		const mat3& m = v.mat3Value;
		if (isZero(m))
		{
			s += to_nice_string(0.0);
		}
		else if (isIdentity(m))
		{
			s += to_nice_string(1.0);
		}
		else
		{
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
				{
					s += to_nice_string(m(i,j));
					if ((i != 2) || (j != 2)) s += ",";
				}
		}
		s += ")";
	}
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

std::string febcode::opToString(BinaryOp op)
{
	switch (op)
	{
	case BinaryOp::Plus:         return "+";
	case BinaryOp::Minus:        return "-";
	case BinaryOp::Multiply:     return "*";
	case BinaryOp::Divide:       return "/";
	case BinaryOp::Exponent:     return "**";
	case BinaryOp::Greater:      return ">";
	case BinaryOp::GreaterEqual: return ">=";
	case BinaryOp::Less:         return "<";
	case BinaryOp::LessEqual:    return "<=";
	case BinaryOp::EqualEqual:   return "==";
	case BinaryOp::NotEqual:     return "!=";
	case BinaryOp::AndAnd:       return "&&";
	case BinaryOp::OrOr:         return "||";
	default:
		return "unknown_op";
	}
}

[[noreturn]]
void febcode::error(SourceLocation loc, const std::string& msg)
{
	throw FatalError(msg, loc);
}

[[noreturn]]
void febcode::error(const ASTNode* node, const std::string& msg)
{
	throw FatalError(msg, (node ? node->location : SourceLocation()));
}

std::string febcode::StatementTypeToString(febcode::StatementType type)
{
	std::string stmtType;
	switch (type)
	{
	case StatementType::ExpressionStatement: stmtType = "ExpressionStatement"; break;
	case StatementType::VarDeclStatement   : stmtType = "VarDeclStatement"; break;
	case StatementType::ReturnStatement    : stmtType = "ReturnStatement"; break;
	case StatementType::BlockStatement     : stmtType = "BlockStatement"; break;
	case StatementType::IfStatement        : stmtType = "IfStatement"; break;
	case StatementType::WhileStatement     : stmtType = "WhileStatement"; break;
	case StatementType::ForStatement       : stmtType = "ForStatement"; break;
	case StatementType::FunctionStatement  : stmtType = "FunctionStatement"; break;
	case StatementType::StructStatement    : stmtType = "StructStatement"; break;
	default:
		assert(false);
		stmtType = "UnknownStatement"; break;
	}
	return stmtType;
}
