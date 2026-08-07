#include "modifier.h"
using namespace febcode;

ExprPtr Modifier::Call(const std::string& name, const std::vector<ExprPtr>& args) 
{
	std::vector<ExprPtr> copyArgs;
	std::vector<Type> argTypes;
	for (const auto& arg : args)
	{
		copyArgs.emplace_back(clone(arg.get()));
		argTypes.push_back(arg->valType);
	}


	Type returnType;
	// see if this is a function defined in the program
	int index = prg.resolveFunction(name, argTypes);
	if (index >= 0)
	{
		returnType = prg.functions[index].returnType;
	}
	else
	{
		// this could be a matrix variable. Unfortunately, we can't tell that here. 
		// If there are two args and they are of type int then we assume this is a matrix component extraction
		if ((args.size() == 2) && (argTypes[0] == prg.types.Int()) && (argTypes[1] == prg.types.Int()))
		{
			returnType = prg.types.Double(); // matrix components are always doubles
		}
		else
		{
			throw std::runtime_error("Undefined function: " + name);
		}
	}

	ExprPtr var = std::make_unique<VariableExpr>(name);

	ExprPtr c = std::make_unique<CallExpr>(std::move(var), std::move(copyArgs));
	c->valType = returnType;
	return c;
}

ExprPtr Modifier::Index(const ExprPtr& object, const ExprPtr& index)
{
	if (object->valType && object->valType->kind == TypeKind::Array)
	{
		Type elementType = object->valType->elementType;
		ExprPtr i = std::make_unique<IndexExpr>(clone(object.get()), clone(index.get()));
		i->valType = elementType;
		return i;
	}
	else
	{
		throw std::runtime_error("Cannot index non-array type");
	}
}

ExprPtr Modifier::Member(const ExprPtr& object, const std::string& property)
{
	if (object->valType && object->valType->kind == TypeKind::Struct)
	{
		const auto& fields = object->valType->fields;
		auto it = std::find_if(fields.begin(), fields.end(), [&](const auto& field) { return field.second == property; });
		if (it != fields.end())
		{
			Type fieldType = it->first;
			ExprPtr m = std::make_unique<MemberExpr>(clone(object.get()), property);
			m->valType = fieldType;
			return m;
		}
		else
		{
			throw std::runtime_error("Struct type does not have member: " + property);
		}
	}
	else if (object->valType && object->valType->kind == TypeKind::Vec2)
	{
		if (property == "x" || property == "y")
		{
			Type fieldType = prg.types.Double();
			ExprPtr m = std::make_unique<MemberExpr>(clone(object.get()), property);
			m->valType = fieldType;
			return m;
		}
		else
		{
			throw std::runtime_error("Vec2 type does not have member: " + property);
		}
	}
	else if (object->valType && object->valType->kind == TypeKind::Vec3)
	{
		if (property == "x" || property == "y" || property == "z")
		{
			Type fieldType = prg.types.Double();
			ExprPtr m = std::make_unique<MemberExpr>(clone(object.get()), property);
			m->valType = fieldType;
			return m;
		}
		else
		{
			throw std::runtime_error("Vec3 type does not have member: " + property);
		}
	}
	else
	{
		throw std::runtime_error("Cannot access member of non-struct type");
	}
}

ExprPtr Modifier::Initializer(Type type)
{
	if ((type == nullptr) || (type->kind != TypeKind::Array) || (type->elementType == nullptr))
	{
		throw std::runtime_error("Invalid type in Initializer");
	}

	Type elementType = type->elementType;
	size_t arraySize = type->arraySize;

	std::vector<ExprPtr> elements;
	for (int i = 0; i < arraySize; ++i)
	{
		elements.push_back(Zero(elementType));
	}

	ExprPtr init = std::make_unique<InitExpr>(std::move(elements));
	init->valType = type;
	return init;
}

ExprPtr Modifier::Constructor(Type type)
{
	if ((type == nullptr) || (type->kind != TypeKind::Struct))
	{
		throw std::runtime_error("Invalid type in Constructor");
	}

	std::vector<ExprPtr> init;
	for (const auto& field : type->fields)
	{
		init.push_back(Zero(field.first));
	}
	ExprPtr i = std::make_unique<ConstructorExpr>(type, std::move(init));
	return i;
}

ExprPtr Modifier::Zero(Type type)
{
	switch (type->kind)
	{
	case TypeKind::Bool  : return Literal(false);
	case TypeKind::Int   : return Literal(0);
	case TypeKind::Double: return Literal(0.0);
	case TypeKind::Vec2  : return Literal(vec2());
	case TypeKind::Vec3  : return Literal(vec3());
	case TypeKind::Mat2  : return Literal(mat2());
	case TypeKind::Mat3  : return Literal(mat3());
	case TypeKind::Array : return Initializer(type);
	case TypeKind::Struct: return Constructor(type);
	default:
		throw std::runtime_error("Unsupported type for Zero");
	}
}

ExprPtr Modifier::Unary(UnaryOp op, const Expression* arg)
{
	ExprPtr u = std::make_unique<UnaryExpr>(op, clone(arg));
	u->valType = arg->valType;
	return u;
}

ExprPtr Modifier::Binary(BinaryOp op, const Expression* left, const Expression* right)
{
	// get the signature for this operator and operand types (this throws if the operator is not defined for these types)
	BinaryOpSignature sig = prg.resolveBinaryOp(op, left->valType, right->valType);

	if (sig.resultType == nullptr)
	{
		throw std::runtime_error("Invalid binary operator: " + std::to_string(static_cast<int>(op)) + " for types " + std::to_string(static_cast<int>(left->valType->kind)) + " and " + std::to_string(static_cast<int>(right->valType->kind)));
	}

	ExprPtr b = std::make_unique<BinaryExpr>(std::move(clone(left)), op, std::move(clone(right)));
	b->valType = sig.resultType;
	return b;
}

ExprPtr Modifier::Component(const ExprPtr& expr, int component)
{
	if (auto lit = dynamic_cast<const LiteralExpr*>(expr.get()))
	{
		if (isVec2(lit->value)) return Literal(lit->value.vec2Value[component]);
		if (isVec3(lit->value)) return Literal(lit->value.vec3Value(component));
		if (isMat2(lit->value)) return Literal(lit->value.mat2Value(component / 2, component % 2));
		if (isMat3(lit->value)) return Literal(lit->value.mat3Value(component / 3, component % 3));
	}
	else if (auto var = dynamic_cast<const VariableExpr*>(expr.get()))
	{
		if (isVec2Type(var->valType) || isVec3Type(var->valType))
		{
			std::string compName;
			switch (component)
			{
			case 0: compName = "x"; break;
			case 1: compName = "y"; break;
			case 2: compName = "z"; break;
			default: throw std::runtime_error("Invalid component index in extractComponent");
			}

			auto varExpr = std::make_unique<VariableExpr>(var->name);
			varExpr->valType = var->valType;
			return Member(std::move(varExpr), compName);
		}
		else if (isMat2Type(var->valType))
		{
			auto varCopy = std::make_unique<VariableExpr>(var->name);
			varCopy->valType = var->valType;

			int i = component / 2;
			int j = component % 2;

			std::vector<ExprPtr> args;
			args.emplace_back(Literal(i));
			args.emplace_back(Literal(j));

			auto compExpr = std::make_unique<CallExpr>(std::move(varCopy), std::move(args));
			compExpr->valType = prg.types.Double();
			return compExpr;
		}
		else if (isMat3Type(var->valType))
		{
			auto varCopy = std::make_unique<VariableExpr>(var->name);
			varCopy->valType = var->valType;

			int i = component / 3;
			int j = component % 3;

			std::vector<ExprPtr> args;
			args.emplace_back(Literal(i));
			args.emplace_back(Literal(j));

			auto compExpr = std::make_unique<CallExpr>(std::move(varCopy), std::move(args));
			compExpr->valType = prg.types.Double();
			return compExpr;
		}
	}
	throw std::runtime_error("Cannot extract component from expression");
}
