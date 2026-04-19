#include "parser.h"
#include "resolver.h"
#include <iostream>
using namespace febcode;

std::unique_ptr<febcode::Statement> Parser::parseDeclaration() {
	if (match(TokenType::Input)) 
	{
		if (!isType()) {
			throw std::runtime_error("Expected type after 'in'.");
		}

		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			throw std::runtime_error("Unknown type name");

		if (!match(TokenType::Identifier)) {
			throw std::runtime_error("Expected identifier after type.");
		}

		std::string name = lexeme(previous());

		Type varType = type;

		std::vector<size_t> arraySizes;
		while (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				throw std::runtime_error("Expected array size after '['.");

			int size = std::stoul(lexeme(previous()));
			if (size == 0)
				throw std::runtime_error("Array size must be greater than zero.");

			if (!match(TokenType::RightBrack))
				throw std::runtime_error("Expected ']' after array size.");

			arraySizes.push_back(size);
		}

		for (int i = (int)arraySizes.size() - 1; i >= 0; --i)
		{
			varType = prg.types.getArrayType(varType, arraySizes[i]);
		}

		if (!match(TokenType::Semicolon)) {
			throw std::runtime_error("Expected ';' after input declaration.");
		}

		Var var{ name, {}, nullptr, true };
		return std::make_unique<VarDeclStmt>(varType, var);
	}
	else if (isType())
	{
		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			throw std::runtime_error("Unknown type name");

		if (match(TokenType::Identifier)) {
			std::string name = lexeme(previous());

			if (match(TokenType::LeftParen)) {

				if (checkType() || check(TokenType::RightParen))
				{
					// function declaration
					return parseFunctionDeclaration(type, name);
				}

				rewind();
			}

			return parseVarDeclaration(type, name);
		}
		else
			throw std::runtime_error("Expected identifier after type.");
	}
	if (match(TokenType::Struct)) return parseStructDeclaration();
	return parseStatement();
}

std::unique_ptr<febcode::Statement> Parser::parseBlockStatement() {
	auto block = std::make_unique<BlockStmt>();

	while (!check(TokenType::RightBrace) && !isAtEnd()) {
		auto decl = parseDeclaration();
		if (decl) {
			block->statements.push_back(std::move(decl));
		}
	}

	if (!match(TokenType::RightBrace)) {
		throw std::runtime_error("Expected '}' after block.");
	}

	return block;
}

std::unique_ptr<Statement> Parser::parseStatement() {
	if (match(TokenType::Return   )) return parseReturnStatement();
	if (match(TokenType::LeftBrace)) return parseBlockStatement();
	if (match(TokenType::If       )) return parseIfStatement();
	if (match(TokenType::While    )) return parseWhileStatement();
	if (match(TokenType::For      )) return parseForStatement();

	return parseExpressionStatement();
}

std::unique_ptr<Statement> Parser::parseExpressionStatement() {
	auto expr = parseExpression();
	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after expression.");
	}
	return std::make_unique<ExpressionStmt>(std::move(expr));
}

std::unique_ptr<Statement> Parser::parseVarDeclaration(Type type, const std::string& name) 
{
	if (type == prg.types.Void())
		throw std::runtime_error("Variables cannot be of type void.");

	std::vector<Var> vars;
	std::string varName = name;
	while (true)
	{
		Type varType = type;

		std::vector<size_t> arraySizes;
		while (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				throw std::runtime_error("Expected array size after '['.");

			int size = std::stoul(lexeme(previous()));
			if (size == 0)
				throw std::runtime_error("Array size must be greater than zero.");

			if (!match(TokenType::RightBrack))
				throw std::runtime_error("Expected ']' after array size.");

			arraySizes.push_back(size);
		}

		for (int i = (int)arraySizes.size() -1; i >= 0; --i)
		{
			varType = prg.types.getArrayType(varType, arraySizes[i]);
		}

		std::unique_ptr<Expression> initializer;
		if (match(TokenType::Equal)) {
			initializer = parseExpression();
		}
		else if (check(TokenType::LeftParen))
		{
			initializer = parseConstructor(type);
		}

		// make sure the name doesn't start with an underscore, which is reserved for internal use
		if (!varName.empty() && varName[0] == '_')
			throw std::runtime_error("Variable names cannot start with an underscore.");

		vars.push_back({varName, arraySizes, std::move(initializer)});

		if (!match(TokenType::Comma)) break;

		if (!match(TokenType::Identifier))
			throw std::runtime_error("Identifier expected after comma.");

		varName = std::string(lexeme(previous()));
	}

	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after variable declaration.");
	}

	return std::make_unique<VarDeclStmt>(type, vars);
}

std::unique_ptr<Statement> Parser::parseStructDeclaration() {
	if (!match(TokenType::Identifier)) {
		throw std::runtime_error("Expected struct name.");
	}
	std::string name(lexeme(previous()));
	if (!match(TokenType::LeftBrace)) {
		throw std::runtime_error("Expected '{' after struct name.");
	}
	std::vector<std::pair<Type, std::string>> fields;
	while (!check(TokenType::RightBrace) && !isAtEnd()) 
	{
		if (!isType()) {
			throw std::runtime_error("Expected field type in struct.");
		}

		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			throw std::runtime_error("Unknown type name");

		if (!match(TokenType::Identifier)) {
			throw std::runtime_error("Expected field name in struct.");
		}
		std::string fieldName(lexeme(previous()));

		if (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				throw std::runtime_error("Expected array size after '[' in struct field.");
			size_t arraySize = std::stoul(lexeme(previous()));
			if (arraySize == 0)
				throw std::runtime_error("Array size must be greater than zero.");
			if (!match(TokenType::RightBrack))
				throw std::runtime_error("Expected ']' after array size in struct field.");

			type = prg.types.getArrayType(type, arraySize);
		}

		if (!match(TokenType::Semicolon)) {
			throw std::runtime_error("Expected ';' after struct field declaration.");
		}

		fields.push_back({ type, fieldName });
	}
	if (!match(TokenType::RightBrace)) {
		throw std::runtime_error("Expected '}' after struct body.");
	}

	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after struct declaration.");
	}

	Type type = prg.types.defineStructType(name, fields);

	return std::make_unique<StructStmt>(name, type, std::move(fields));
}

std::unique_ptr<Statement> Parser::parseFunctionDeclaration(Type type, const std::string& name) 
{
	std::vector<std::pair<Type, std::string>> parameters;

	if (!check(TokenType::RightParen)) {
		do {
			if (!isType()) {
				throw std::runtime_error("Expect type.");
			}

			Type paramType = prg.types.getType(lexeme(previous()));
			if (paramType == nullptr)
				throw std::runtime_error("Unknown type name");

			if (!match(TokenType::Identifier)) {
				throw std::runtime_error("Expected parameter name.");
			}

			std::string param = lexeme(previous());

			if (match(TokenType::LeftBrack))
			{
				if (!match(TokenType::Integer))
					throw std::runtime_error("Expected array size after '[' in parameter.");
				size_t arraySize = std::stoul(lexeme(previous()));
				if (arraySize == 0)
					throw std::runtime_error("Array size must be greater than zero.");
				if (!match(TokenType::RightBrack))
					throw std::runtime_error("Expected ']' after array size in parameter.");
				paramType = prg.types.getArrayType(paramType, arraySize);
			}

			parameters.push_back({ paramType, param });
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen)) {
		throw std::runtime_error("Expected ')' after parameters.");
	}

	if (!match(TokenType::LeftBrace)) {
		throw std::runtime_error("Expected '{' before function body.");
	}

	auto body = parseBlockStatement();

	return std::make_unique<FunctionStmt>(
		name,
		type,
		std::move(parameters),
		std::move(body)
	);
}

std::unique_ptr<Statement> Parser::parseReturnStatement() {
	std::unique_ptr<Expression> value = nullptr;
	if (!check(TokenType::Semicolon)) {
		value = parseExpression();
	}

	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after return value.");
	}

	return std::make_unique<febcode::ReturnStmt>(std::move(value));
}

std::unique_ptr<febcode::Statement> Parser::parseIfStatement() {
	if (!match(TokenType::LeftParen)) {
		throw std::runtime_error("Expected '(' after 'if'.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::RightParen)) {
		throw std::runtime_error("Expected ')' after if condition.");
	}

	auto thenBranch = parseStatement(); // could be block or single statement

	std::unique_ptr<Statement> elseBranch = nullptr;
	if (match(TokenType::Else)) {
		elseBranch = parseStatement();
	}

	return std::make_unique<IfStmt>(std::move(condition),
		std::move(thenBranch),
		std::move(elseBranch));
}

std::unique_ptr<febcode::Statement> Parser::parseWhileStatement() {
	if (!match(TokenType::LeftParen)) {
		throw std::runtime_error("Expected '(' after 'while'.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::RightParen)) {
		throw std::runtime_error("Expected ')' after while condition.");
	}

	auto body = parseStatement(); // can be block or single statement

	return std::make_unique<WhileStmt>(
		std::move(condition),
		std::move(body));
}

std::unique_ptr<febcode::Statement> Parser::parseForStatement() {
	if (!match(TokenType::LeftParen)) {
		throw std::runtime_error("Expected '(' after 'for'.");
	}

	auto init = parseDeclaration(); // initializer (can be var decl, expression stmt, or empty)
	if (!init && !isVarDecl(init) && !isExprStmt(init)) {
		throw std::runtime_error("Expected variable declaration, expression statement, or ';' in for loop initializer.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after for initializer.");
	}

	auto increment = parseExpression();

	if (!match(TokenType::RightParen)) {
		throw std::runtime_error("Expected ')' after for increment.");
	}

	auto body = parseStatement(); // can be block or single statement

	return std::make_unique<ForStmt>(
		std::move(init),
		std::move(condition),
		std::move(increment),
		std::move(body));
}

std::unique_ptr<Expression> Parser::parseAssignment() {
	auto expr = parseOr(); // parse the left-hand side first

	// Check if this is an assignment
	if (match(TokenType::Equal)) {
		const Token& equalsToken = previous();

		auto value = parseAssignment();

		return std::make_unique<AssignExpr>(
			std::move(expr),
			std::move(value)
		);
	}

	return expr; // just an equality/expression if no '='
}

std::unique_ptr<Expression> Parser::parseOr()
{
	auto expr = parseAnd();

	while (match(TokenType::OrOr))
	{
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseAnd();
		expr = std::make_unique<BinaryExpr>(
			std::move(expr),
			op,
			std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseAnd()
{
	auto expr = parseEquality();

	while (match(TokenType::AndAnd))
	{
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseEquality();
		expr = std::make_unique<BinaryExpr>(
			std::move(expr),
			op,
			std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseEquality() {
	auto expr = parseComparison();

	while (match(TokenType::EqualEqual) || match(TokenType::NotEqual)) {
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseComparison();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseComparison() {
	auto expr = parseTerm();

	while (match(TokenType::Greater) || match(TokenType::GreaterEqual) ||
		match(TokenType::Less) || match(TokenType::LessEqual)) {
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseTerm();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseTerm() {
	auto expr = parseFactor();

	while (match(TokenType::Plus) || match(TokenType::Minus)) {
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseFactor();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseFactor() {
	auto expr = parseUnary();

	while (match(TokenType::Star) || match(TokenType::Slash)) {
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseUnary();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseUnary() {
	if (match(TokenType::Minus) || match(TokenType::Not)) {
		UnaryOp op = tokenToUnaryOp(previous());
		auto right = parseUnary();

		// absorbe negative sign for scalars.
		if (op == UnaryOp::Negate) {

			int n;
			if (isInt(right.get(), n)) return std::make_unique<LiteralExpr>(-n);
			double d;
			if (isDouble(right.get(), d)) return std::make_unique<LiteralExpr>(-d);
		}

		return std::make_unique<UnaryExpr>(op, std::move(right));
	}

	return parseExponent();
}

std::unique_ptr<Expression> Parser::parseExponent() {
	auto left = parseCall();

	if (match(TokenType::Exponent)) {
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseUnary();
		return std::make_unique<BinaryExpr>(
			std::move(left), op, std::move(right));
	}

	return left;
}

std::unique_ptr<Expression> Parser::parseCall() {

	if (isType())
	{
		std::string typeName = lexeme(previous());
		Type type = prg.types.getType(typeName);
		if (type == nullptr)
			throw std::runtime_error("Unknown type name in constructor expression.");

		return parseConstructor(type);
	}

	auto expr = parsePrimary();

	while (true) {
		if (match(TokenType::LeftParen)) {
			expr = finishCall(std::move(expr));
		}
		else if (match(TokenType::Dot))
		{
			if (!match(TokenType::Identifier)) {
				throw std::runtime_error("Expected property name after '.'.");
			}

			expr = std::make_unique<MemberExpr>(
				std::move(expr),
				std::string(previous().start, previous().length)
			);
		}
		else if (match(TokenType::LeftBrack)) {
			auto index = parseExpression();
			if (!match(TokenType::RightBrack))
				throw std::runtime_error("Expect ']' after index.");

			expr = std::make_unique<IndexExpr>(
				std::move(expr),
				std::move(index));
		}
		else {
			break;
		}
	}
	return expr;
}

std::unique_ptr<Expression> Parser::parseConstructor(Type type)
{
	std::vector<std::unique_ptr<Expression>> arguments;

	if (!match(TokenType::LeftParen))
		throw std::runtime_error("Expected '(' after type name in constructor expression.");

	if (!check(TokenType::RightParen)) {
		do {
			arguments.push_back(parseExpression());
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen))
		throw std::runtime_error("Expected ')' after arguments.");

	return std::make_unique<ConstructorExpr>(type, std::move(arguments));
}

std::unique_ptr<Expression> Parser::finishCall(std::unique_ptr<Expression> callee) 
{
	auto* varExpr = dynamic_cast<VariableExpr*>(callee.get());
	if (!varExpr) {
		throw std::runtime_error("Can only call functions by name.");
	}

	std::vector<std::unique_ptr<Expression>> arguments;

	if (!check(TokenType::RightParen)) {
		do {
			arguments.push_back(parseExpression());
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen))
		throw std::runtime_error("Expected ')' after arguments.");

	return std::make_unique<CallExpr>(
		varExpr->name,
		std::move(arguments)
	);
}

std::unique_ptr<Expression> Parser::parsePrimary() 
{
	if (match(TokenType::Integer)) {
		int value = std::stoi(std::string(previous().start, previous().length));
		return std::make_unique<LiteralExpr>(value);
	}

	if (match(TokenType::Double)) {
		double value = std::stod(std::string(previous().start, previous().length));
		return std::make_unique<LiteralExpr>(value);
	}

	if (match(TokenType::True))
		return std::make_unique<LiteralExpr>(true);
	
	if (match(TokenType::False))
		return std::make_unique<LiteralExpr>(false);

	if (match(TokenType::Identifier)) {
		std::string name(previous().start, previous().length);

		// check some special built-in constants like "PI"
		if (name == "PI") {
			return std::make_unique<LiteralExpr>(3.14159265358979323846);
		}

		// otherwise, it's a variable
		return std::make_unique<VariableExpr>(
			std::string(previous().start, previous().length));
	}

	if (match(TokenType::LeftParen)) {
		auto expr = parseExpression();
		if (!match(TokenType::RightParen)) {
			throw std::runtime_error("Expected ')' after expression.");
		}
		return expr;
	}

	if (match(TokenType::LeftBrace)) 
	{
		std::vector<std::unique_ptr<Expression>> elements;

		if (!check(TokenType::RightBrace)) {
			do {
				elements.push_back(parseExpression());
			} while (match(TokenType::Comma));
		}

		if (!match(TokenType::RightBrace))
			throw std::runtime_error("Expected '}' after initializer list.");

		return std::make_unique<InitExpr>(std::move(elements));
	}

	throw std::runtime_error("Expected expression.");
}

//---------------------------------------------------------------------
static int l = 0;
static void printTabs()
{
	for (int i = 0; i < l; ++i) std::cout << "    ";
}

static void printTabs(std::ostream& os)
{
	for (int i = 0; i < l; ++i) os << "    ";
}

std::ostream& operator << (std::ostream& o, const UnaryOp& op)
{
	switch (op)
	{
	case UnaryOp::Negate: return o << "-";
	case UnaryOp::Not: return o << "!";
	default: return o << "<unknown op>";
	}
}

std::ostream& operator << (std::ostream& o, const std::vector<std::string>& v)
{
	o << "[";
	for (size_t i = 0; i < v.size(); ++i)
	{
		o << "\"" << v[i] << "\"";
		if (i != v.size() - 1) o << ", ";
	}
	o << "]";
	return o;
}

std::ostream& operator << (std::ostream& o, const BinaryOp& op)
{
	switch (op)
	{
	case BinaryOp::Plus        : return o << "+";
	case BinaryOp::Minus       : return o << "-";
	case BinaryOp::Multiply    : return o << "*";
	case BinaryOp::Divide      : return o << "/";
	case BinaryOp::Exponent    : return o << "**";
	case BinaryOp::EqualEqual  : return o << "==";
	case BinaryOp::NotEqual    : return o << "!=";
	case BinaryOp::Greater     : return o << ">";
	case BinaryOp::GreaterEqual: return o << ">=";
	case BinaryOp::Less        : return o << "<";
	case BinaryOp::LessEqual   : return o << "<=";
	case BinaryOp::AndAnd      : return o << "&&";
	case BinaryOp::OrOr        : return o << "||";
	default: return o << "<unknown op>";
	}
}

static void printExpr(const Expression* e);

static void printLiteralExpr(const LiteralExpr* e)
{
	std::cout << ValueToString(e->value);
}

static void printVariableExpr(const VariableExpr* e)
{
	std::cout << e->name;
}

static void printMemberExpr(const MemberExpr* e)
{
	std::cout << "MemberExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "property: " << e->property << "\n";
	l--; printTabs(); std::cout << "}";
}

static void printAssignmentExpr(const AssignExpr* e)
{
	std::cout << "AssignExpr {\n"; l++;
	printTabs(); std::cout << "target: "; printExpr(e->target.get()); std::cout << ",\n";
	printTabs(); std::cout << "value: "; printExpr(e->value.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printUnaryExpr(const UnaryExpr* e)
{
	std::cout << "UnaryExpr {\n"; l++;
	printTabs(); std::cout << "op: " << e->op << ",\n";
	printTabs(); std::cout << "right: "; printExpr(e->right.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printBinaryExpr(const BinaryExpr* e)
{
	std::cout << "BinaryExpr {\n"; l++;
	printTabs(); std::cout << "op: " << e->op << ",\n";
	printTabs(); std::cout << "left: "; printExpr(e->left.get()); std::cout << ",\n";
	printTabs(); std::cout << "right: "; printExpr(e->right.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printExprList(const std::vector<std::unique_ptr<Expression>>& a)
{
	std::cout << "[\n"; l++;
	for (size_t i = 0; i < a.size(); ++i)
	{
		printTabs(); printExpr(a[i].get());
		if (i != a.size() - 1) std::cout << ",\n";
		else std::cout << "\n";
	}
	l--; printTabs(); std::cout << "]";
}

static void printCallExpr(const CallExpr* e)
{
	std::cout << "CallExpr {\n"; l++;
	printTabs(); std::cout << "callee: " << e->name << ",\n";
	printTabs(); std::cout << "args: "; printExprList(e->arguments); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printInitializerExpr(const InitExpr* e)
{
	std::cout << "InitializerExpr {\n"; l++;
	printTabs(); std::cout << "elements: "; printExprList(e->elements); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printIndexExpr(const IndexExpr* e)
{
	std::cout << "IndexExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "index: "; printExpr(e->index.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printConstructorExpr(const ConstructorExpr* e)
{
	std::cout << "ConstructorExpr {\n"; l++;
	printTabs(); std::cout << "type: " << TypeToString(e->valType) << ",\n";
	printTabs(); std::cout << "args: "; printExprList(e->args); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printExpr(const Expression* e)
{
	if      (auto l = dynamic_cast<const LiteralExpr*      >(e)) printLiteralExpr      (l);
	else if (auto v = dynamic_cast<const VariableExpr*     >(e)) printVariableExpr     (v);
	else if (auto m = dynamic_cast<const MemberExpr*       >(e)) printMemberExpr       (m);
	else if (auto a = dynamic_cast<const AssignExpr*       >(e)) printAssignmentExpr   (a);
	else if (auto u = dynamic_cast<const UnaryExpr*        >(e)) printUnaryExpr        (u);
	else if (auto b = dynamic_cast<const BinaryExpr*       >(e)) printBinaryExpr       (b);
	else if (auto c = dynamic_cast<const CallExpr*         >(e)) printCallExpr         (c);
	else if (auto c = dynamic_cast<const InitExpr*         >(e)) printInitializerExpr  (c);
	else if (auto c = dynamic_cast<const IndexExpr*        >(e)) printIndexExpr        (c);
	else if (auto c = dynamic_cast<const ConstructorExpr*  >(e)) printConstructorExpr  (c);
	else if (e == nullptr)
		std::cout << "null";
	else
		std::cout << "(Unknown Expression)";
}

static void printExpressionStmt(const ExpressionStmt* s)
{
	std::cout << "ExpressionStmt: {\n"; l++;
	printTabs(); std::cout << "expr: "; printExpr(s->expr.get());
	std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printStatement(const Statement* stmt);

static std::ostream& operator << (std::ostream& o, TypeKind type)
{
	switch (type)
	{
	case TypeKind::Void  : return o << "void";
	case TypeKind::Bool  : return o << "bool";
	case TypeKind::Int   : return o << "int";
	case TypeKind::Double: return o << "double";
	case TypeKind::Vec2  : return o << "vec2";
	case TypeKind::Vec3  : return o << "vec3";
	case TypeKind::Mat2  : return o << "mat2";
	case TypeKind::Mat3  : return o << "mat3";
	default:
		return o << "<unknown type>";
	}
}

static void printVarDeclStmt(const VarDeclStmt* s)
{
	std::cout << "VarDeclStmt: {\n"; l++;
	printTabs(); std::cout << "type: " << TypeToString(s->type) << ",\n";
	printTabs(); std::cout << "vars: [\n"; l++;
	for (const auto& var : s->vars) {
		printTabs(); std::cout << "{ name: " << var.name;
		if (!var.arraySizes.empty())
		{
			std::cout << ", size: [";
			for (size_t i = 0; i < var.arraySizes.size(); ++i)
			{
				std::cout << var.arraySizes[i];
				if (i != var.arraySizes.size() - 1) std::cout << "][";
			}
			std::cout << "]";
		}
		std::cout << ", initializer: ";
		printExpr(var.initializer.get());
		std::cout << " },\n";
	}
	l--; printTabs(); std::cout << "],\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printReturnStmt(const ReturnStmt* s)
{
	std::cout << "ReturnStmt: {\n"; l++;
	printTabs(); std::cout << "value: "; printExpr(s->value.get());
	std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printBlockStmt(const BlockStmt* s)
{
	std::cout << "BlockStmt: [\n"; l++;
	size_t n = s->statements.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& stmt = s->statements[i];
		printStatement(stmt.get());
		if (i != n - 1) std::cout << ",\n";
	}	std::cout << "\n";
	l--;
	printTabs(); std::cout << "]";
}

static void printIfStmt(const IfStmt* s)
{
	std::cout << "IfStmt: {\n"; l++;
	printTabs(); std::cout << "condition: "; printExpr(s->condition.get()); std::cout << ",\n";
	printTabs(); std::cout << "thenBranch: {\n"; l++; printStatement(s->thenBranch.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
	if (s->elseBranch)
	{
		printTabs(); std::cout << "elseBranch: {\n"; l++;
		printStatement(s->elseBranch.get());
		std::cout << "\n";
		l--; printTabs(); std::cout << "}\n";
	}
	l--;
	printTabs(); std::cout << "}";
}

static void printWhileStmt(const WhileStmt* s)
{
	std::cout << "WhileStmt: {\n"; l++;
	printTabs(); std::cout << "condition: "; printExpr(s->condition.get()); std::cout << ",\n";
	printTabs(); std::cout << "body: {\n"; l++; printStatement(s->body.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
}

static void printForStmt(const ForStmt* s)
{
	std::cout << "ForStmt: {\n"; l++;
	printTabs(); std::cout << "initializer: "; printStatement(s->initializer.get()); std::cout << ",\n";
	printTabs(); std::cout << "condition: "; printExpr(s->condition.get()); std::cout << ",\n";
	printTabs(); std::cout << "increment: "; printExpr(s->increment.get()); std::cout << ",\n";
	printTabs(); std::cout << "body: {\n"; l++; printStatement(s->body.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
	l--; printTabs(); std::cout << "}";
}

static void printFunctionStmt(const FunctionStmt* s)
{
	std::cout << "FunctionStmt: {\n"; l++;
	printTabs(); std::cout << "type: " << TypeToString(s->returnType) << ",\n";
	printTabs(); std::cout << "name: " << s->name << ",\n";
	printTabs(); std::cout << "params: [\n"; l++;
	for (const auto& param : s->params)
	{
		printTabs(); std::cout << "{ type: " << TypeToString(param.first) << ", name: " << param.second << " },\n";
	}
	l--; printTabs(); std::cout << "],\n";
	printTabs(); std::cout << "body: {\n"; l++; printStatement(s->body.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
	l--; printTabs(); std::cout << "}";
}

static void printStructStmt(const StructStmt* s)
{
	std::cout << "StructStmt: {\n"; l++;
	printTabs(); std::cout << "name: " << s->name << ",\n";
	printTabs(); std::cout << "fields: [\n"; l++;
	for (const auto& field : s->fields)
	{
		printTabs(); std::cout << "{ name: " << field.second << ", type: " << TypeToString(field.first) << " },\n";
	}
	l--; printTabs(); std::cout << "]\n";
	l--; printTabs(); std::cout << "}";
}

static void printStatement(const Statement* stmt)
{
	printTabs();
	if      (auto e = dynamic_cast<const ExpressionStmt*>(stmt)) printExpressionStmt(e);
	else if (auto v = dynamic_cast<const VarDeclStmt*   >(stmt)) printVarDeclStmt   (v);
	else if (auto r = dynamic_cast<const ReturnStmt*    >(stmt)) printReturnStmt    (r);
	else if (auto b = dynamic_cast<const BlockStmt*     >(stmt)) printBlockStmt     (b);
	else if (auto i = dynamic_cast<const IfStmt*        >(stmt)) printIfStmt        (i);
	else if (auto w = dynamic_cast<const WhileStmt*     >(stmt)) printWhileStmt     (w);
	else if (auto l = dynamic_cast<const ForStmt*       >(stmt)) printForStmt       (l);
	else if (auto f = dynamic_cast<const FunctionStmt*  >(stmt)) printFunctionStmt  (f);
	else if (auto s = dynamic_cast<const StructStmt*    >(stmt)) printStructStmt    (s);
	else
		std::cout << "(Unknown Statement)";
}

void febcode::printAST(const AST& ast)
{
	l = 0;
	size_t n = ast.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto stmt = ast[i];
		printStatement(stmt);
		if (i != n - 1) std::cout << ",\n";
	}

	std::cout << std::endl;
}

void febcode::ParseSource(Program& prg, const std::string& source)
{
	Scanner scanner(source);
	std::vector<Token> tokens = scanner.scanTokens();

	Parser parser(prg);
	parser.parse(tokens);

	Resolver resolver(prg);
	resolver.resolve();
}

//---------------------------------------------------------------------
// pretty print functions for debugging purposes

static void prettyPrintLiteralExpr(std::ostream& os, const LiteralExpr& expr)
{
	os << ValueToString(expr.value);
}

static void prettyPrintVariableExpr(std::ostream& os, const VariableExpr& expr)
{
	os << expr.name;
}

static void prettyPrintMemberExpr(std::ostream& os, const MemberExpr& expr)
{
	prettyPrintExpression(os, *expr.object);
	os << "." << expr.property;
}

static void prettyPrintAssignmentExpr(std::ostream& os, const AssignExpr& expr)
{
	prettyPrintExpression(os, *expr.target);
	os << " = ";
	prettyPrintExpression(os, *expr.value);
}

static void prettyPrintUnaryExpr(std::ostream& os, const UnaryExpr& expr)
{
	os << expr.op;
	os << "(";
	prettyPrintExpression(os, *expr.right);
	os << ")";
}

static void prettyPrintBinaryExpr(std::ostream& os, const BinaryExpr& expr)
{
	bool leftBinary  = dynamic_cast<const BinaryExpr*>(expr.left.get()) != nullptr;
	bool rightBinary = dynamic_cast<const BinaryExpr*>(expr.right.get()) != nullptr;

	if (leftBinary) os << "(";
	prettyPrintExpression(os, *expr.left);
	if (leftBinary) os << ")";
	os << expr.op;
	if (rightBinary) os << "(";
	prettyPrintExpression(os, *expr.right);
	if (rightBinary) os << ")";
}

static void prettyPrintCallExpr(std::ostream& os, const CallExpr& expr)
{
	os << expr.name;
	os << "(";
	size_t n = expr.arguments.size();
	for (size_t i = 0; i < n; ++i)
	{
		prettyPrintExpression(os, *expr.arguments[i]);
		if (i != n - 1) os << ", ";
	}
	os << ")";
}

static void prettyPrintInitializerExpr(std::ostream& os, const InitExpr& expr)
{
	os << "{ ";
	size_t n = expr.elements.size();
	for (size_t i = 0; i < n; ++i)
	{
		prettyPrintExpression(os, *expr.elements[i]);
		if (i != n - 1) os << ", ";
	}
	os << " }";
}

static void prettyPrintIndexExpr(std::ostream& os, const IndexExpr& expr)
{
	prettyPrintExpression(os, *expr.object);
	os << "[";
	prettyPrintExpression(os, *expr.index);
	os << "]";
}

static void prettyPrintConstructorExpr(std::ostream& os, const ConstructorExpr& expr)
{
	os << TypeToString(expr.valType);
	os << "(";
	size_t n = expr.args.size();
	for (size_t i = 0; i < n; ++i)
	{
		prettyPrintExpression(os, *expr.args[i]);
		if (i != n - 1) os << ", ";
	}
	os << ")";
}

void febcode::prettyPrintExpression(std::ostream& os, const Expression& expr)
{
	if      (auto l = dynamic_cast<const LiteralExpr*      >(&expr)) prettyPrintLiteralExpr      (os, *l);
	else if (auto v = dynamic_cast<const VariableExpr*     >(&expr)) prettyPrintVariableExpr     (os, *v);
	else if (auto m = dynamic_cast<const MemberExpr*       >(&expr)) prettyPrintMemberExpr       (os, *m);
	else if (auto a = dynamic_cast<const AssignExpr*       >(&expr)) prettyPrintAssignmentExpr   (os, *a);
	else if (auto u = dynamic_cast<const UnaryExpr*        >(&expr)) prettyPrintUnaryExpr        (os, *u);
	else if (auto b = dynamic_cast<const BinaryExpr*       >(&expr)) prettyPrintBinaryExpr       (os, *b);
	else if (auto c = dynamic_cast<const CallExpr*         >(&expr)) prettyPrintCallExpr         (os, *c);
	else if (auto c = dynamic_cast<const InitExpr*         >(&expr)) prettyPrintInitializerExpr  (os, *c);
	else if (auto c = dynamic_cast<const IndexExpr*        >(&expr)) prettyPrintIndexExpr        (os, *c);
	else if (auto c = dynamic_cast<const ConstructorExpr*  >(&expr)) prettyPrintConstructorExpr  (os, *c);
	else
		os << "(Unknown Expression)";
}

static void prettyPrintStatement(std::ostream& os, const Statement& stmt);

static void prettyPrintExpressionStmt(std::ostream& os, const ExpressionStmt& stmt)
{
	prettyPrintExpression(os, *stmt.expr);
	os << ";";
}

static void prettyPrintVarDeclStmt(std::ostream& os, const VarDeclStmt& stmt)
{
	size_t n = stmt.vars.size();
	os << TypeToString(stmt.type) << " ";
	for (size_t i = 0; i < n; ++i)
	{
		const auto& var = stmt.vars[i];
		os << var.name;
		if (!var.arraySizes.empty())
		{
			for (size_t j = 0; j < var.arraySizes.size(); ++j)
			{
				os << "[";
				os << var.arraySizes[j];
				os << "]";
			}
		}
		if (var.initializer)
		{
			os << " = ";
			prettyPrintExpression(os, *var.initializer);
		}
		if (i != n - 1) os << ", ";
	}
	os << ";";
}

static void prettyPrintReturnStmt(std::ostream& os, const ReturnStmt& stmt)
{
	os << "return";
	if (stmt.value)
	{
		os << " ";
		prettyPrintExpression(os, *stmt.value);
	}
	os << ";";
}

static void prettyPrintBlockStmt(std::ostream& os, const BlockStmt& stmt)
{
	os << "{\n"; l++;
	size_t n = stmt.statements.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& s = stmt.statements[i];
		prettyPrintStatement(os, *s);
		if (i != n - 1) os << "\n";
	}
	os << "\n"; l--;
	printTabs(os); os << "}";
}

static void prettyPrintIfStmt(std::ostream& os, const IfStmt& stmt)
{
	os << "if (";
	prettyPrintExpression(os, *stmt.condition);
	os << ")\n";
	printTabs(os); prettyPrintStatement(os, *stmt.thenBranch);
	if (stmt.elseBranch)
	{
		os << "\n";
		printTabs(os); os << "else\n";
		printTabs(os); prettyPrintStatement(os, *stmt.elseBranch);
	}
}

static void prettyPrintWhileStmt(std::ostream& os, const WhileStmt& stmt)
{
	os << "while (";
	prettyPrintExpression(os, *stmt.condition);
	os << ")\n";
	printTabs(os); prettyPrintStatement(os, *stmt.body);
}

static void prettyPrintFunctionStmt(std::ostream& os, const FunctionStmt& stmt)
{
	os << TypeToString(stmt.returnType) << " " << stmt.name << "(";
	size_t n = stmt.params.size();
	for (size_t i = 0; i < n; ++i)
	{
		const auto& param = stmt.params[i];
		os << TypeToString(param.first) << " " << param.second;
		if (i != n - 1) std::cout << ", ";
	}
	os << ")\n";
	printTabs(os); prettyPrintStatement(os, *stmt.body);
}

static void prettyPrintStructStmt(std::ostream& os, const StructStmt& stmt)
{
	os << "struct " << stmt.name << " {\n";
	for (auto& field : stmt.fields)
	{
		os << "    " << TypeToString(field.first) << " " << field.second << ";\n";
	}
	os << "};";
}

static void prettyPrintStatement(std::ostream& os, const febcode::Statement& stmt)
{
	printTabs(os);
	if      (auto e = dynamic_cast<const ExpressionStmt*>(&stmt)) prettyPrintExpressionStmt(os, *e);
	else if (auto v = dynamic_cast<const VarDeclStmt*   >(&stmt)) prettyPrintVarDeclStmt   (os, *v);
	else if (auto r = dynamic_cast<const ReturnStmt*    >(&stmt)) prettyPrintReturnStmt    (os, *r);
	else if (auto b = dynamic_cast<const BlockStmt*     >(&stmt)) prettyPrintBlockStmt     (os, *b);
	else if (auto i = dynamic_cast<const IfStmt*        >(&stmt)) prettyPrintIfStmt        (os, *i);
	else if (auto w = dynamic_cast<const WhileStmt*     >(&stmt)) prettyPrintWhileStmt     (os, *w);
	else if (auto f = dynamic_cast<const FunctionStmt*  >(&stmt)) prettyPrintFunctionStmt  (os, *f);
	else if (auto s = dynamic_cast<const StructStmt*    >(&stmt)) prettyPrintStructStmt    (os, *s);
	else
		std::cout << "(Unknown Statement)";
}

void febcode::prettyPrintAST(std::ostream& os, const AST& ast)
{
	l = 0;
	size_t n = ast.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& stmt = *ast[i];
		prettyPrintStatement(os, stmt);
		if (i != n - 1) os << "\n";
	}
}

void febcode::prettyPrintAST(const AST& ast)
{
	prettyPrintAST(std::cout, ast);
}
