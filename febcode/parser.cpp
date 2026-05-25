#include "parser.h"
#include "resolver.h"
#include <iostream>
using namespace febcode;

std::unique_ptr<febcode::Statement> Parser::parseDeclaration() {

	SourceLocation currentLoc = peek().loc;

	if (match(TokenType::Input)) 
	{
		// make sure we are at the global scope, since inputs can only be declared at the global level
		if (scopeDepth > 0) {
			error(currentLoc, "Input declarations can only be made at the global scope.");
		}

		if (!isType()) {
			error(currentLoc, "Expected type after 'in'.");
		}

		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			error(currentLoc, "Unknown type name");

		if (!match(TokenType::Identifier)) {
			error(currentLoc, "Expected identifier after type.");
		}

		std::string name = lexeme(previous());

		Type varType = type;

		std::vector<size_t> arraySizes;
		while (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				error(currentLoc, "Expected array size after '['.");

			int size = std::stoul(lexeme(previous()));
			if (size == 0)
				error(currentLoc, "Array size must be greater than zero.");

			if (!match(TokenType::RightBrack))
				error(currentLoc, "Expected ']' after array size.");

			arraySizes.push_back(size);
		}

		for (int i = (int)arraySizes.size() - 1; i >= 0; --i)
		{
			varType = prg.types.getArrayType(varType, arraySizes[i]);
		}

		if (!match(TokenType::Semicolon)) {
			error(currentLoc, "Expected ';' after input declaration.");
		}

		VarPtr var = std::make_unique<Var>(Var{ name, varType, nullptr });
		auto stmt = std::make_unique<VarDeclStmt>(type, std::move(var), true);
		stmt->location = currentLoc;
		return stmt;
	}
	else if (isType())
	{
		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			error(currentLoc, "Unknown type name");

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
			error(currentLoc, "Expected identifier after type.");
	}
	if (match(TokenType::Struct)) return parseStructDeclaration();
	return parseStatement();
}

std::unique_ptr<febcode::Statement> Parser::parseBlockStatement() 
{
	scopeDepth++;

	Token startToken = previous(); assert(startToken.type == TokenType::LeftBrace);

	auto block = std::make_unique<BlockStmt>();
	block->location = startToken.loc;

	while (!check(TokenType::RightBrace) && !isAtEnd()) {
		auto decl = parseDeclaration();
		if (decl) {
			block->statements.push_back(std::move(decl));
		}
	}

	if (!match(TokenType::RightBrace)) {
		error(currentLoc, "Expected '}' after block.");
	}

	scopeDepth--;

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

	SourceLocation currentLoc = peek().loc;
	auto expr = parseExpression();
	if (!match(TokenType::Semicolon)) {
		error(currentLoc, "Expected ';' after expression.");
	}
	auto stmt = std::make_unique<ExpressionStmt>(std::move(expr));
	stmt->location = currentLoc;
	return stmt;
}

std::unique_ptr<Statement> Parser::parseVarDeclaration(Type type, const std::string& name) 
{
	Token startToken = previous(); assert(startToken.type == TokenType::Identifier);

	if (type == prg.types.Void())
		error(currentLoc, "Variables cannot be of type void.");

	std::vector<VarPtr> vars;
	std::string varName = name;
	while (true)
	{
		std::vector<size_t> arraySizes;
		while (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				error(currentLoc, "Expected array size after '['.");

			int size = std::stoul(lexeme(previous()));
			if (size == 0)
				error(currentLoc, "Array size must be greater than zero.");

			if (!match(TokenType::RightBrack))
				error(currentLoc, "Expected ']' after array size.");

			arraySizes.push_back(size);
		}

		Type varType = type;
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
			error(currentLoc, "Variable names cannot start with an underscore.");

		vars.push_back(std::make_unique<Var>(Var{ varName, varType, std::move(initializer) }));

		if (!match(TokenType::Comma)) break;

		if (!match(TokenType::Identifier))
			error(currentLoc, "Identifier expected after comma.");

		varName = std::string(lexeme(previous()));
	}

	if (!match(TokenType::Semicolon)) {
		error(currentLoc, "Expected ';' after variable declaration.");
	}

	auto stmt = std::make_unique<VarDeclStmt>(type, std::move(vars));
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<Statement> Parser::parseStructDeclaration() {

	Token startToken = previous(); assert(startToken.type == TokenType::Struct);

	if (!match(TokenType::Identifier)) {
		error(currentLoc, "Expected struct name.");
	}
	std::string name(lexeme(previous()));
	if (!match(TokenType::LeftBrace)) {
		error(currentLoc, "Expected '{' after struct name.");
	}
	std::vector<std::pair<Type, std::string>> fields;
	while (!check(TokenType::RightBrace) && !isAtEnd()) 
	{
		if (!isType()) {
			error(currentLoc, "Expected field type in struct.");
		}

		Type type = prg.types.getType(lexeme(previous()));
		if (type == nullptr)
			error(currentLoc, "Unknown type name.");

		if (!match(TokenType::Identifier)) {
			error(currentLoc, "Expected field name in struct.");
		}
		std::string fieldName(lexeme(previous()));

		if (match(TokenType::LeftBrack))
		{
			if (!match(TokenType::Integer))
				error(currentLoc, "Expected array size after '[' in struct field.");
			size_t arraySize = std::stoul(lexeme(previous()));
			if (arraySize == 0)
				error(currentLoc, "Array size must be greater than zero.");
			if (!match(TokenType::RightBrack))
				error(currentLoc, "Expected ']' after array size in struct field.");

			type = prg.types.getArrayType(type, arraySize);
		}

		if (!match(TokenType::Semicolon)) {
			error(currentLoc, "Expected ';' after struct field declaration.");
		}

		fields.push_back({ type, fieldName });
	}
	if (!match(TokenType::RightBrace)) {
		error(currentLoc, "Expected '}' after struct body.");
	}

	if (!match(TokenType::Semicolon)) {
		error(currentLoc, "Expected ';' after struct declaration.");
	}

	Type type = prg.types.defineStructType(name, fields);

	auto stmt = std::make_unique<StructStmt>(name, type, std::move(fields));
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<Statement> Parser::parseFunctionDeclaration(Type type, const std::string& name) 
{
	Token startToken = previous(); assert(startToken.type == TokenType::LeftParen);

	std::vector<VarPtr> parameters;

	if (!check(TokenType::RightParen)) {
		do {
			if (!isType()) {
				error(currentLoc, "Expect type.");
			}

			Type paramType = prg.types.getType(lexeme(previous()));
			if (paramType == nullptr)
				error(currentLoc, "Unknown type name.");

			if (!match(TokenType::Identifier)) {
				error(currentLoc, "Expected parameter name.");
			}

			std::string param = lexeme(previous());

			if (match(TokenType::LeftBrack))
			{
				if (!match(TokenType::Integer))
					error(currentLoc, "Expected array size after '[' in parameter.");
				size_t arraySize = std::stoul(lexeme(previous()));
				if (arraySize == 0)
					error(currentLoc, "Array size must be greater than zero.");
				if (!match(TokenType::RightBrack))
					error(currentLoc, "Expected ']' after array size in parameter.");
				paramType = prg.types.getArrayType(paramType, arraySize);
			}

			parameters.push_back(std::make_unique<Var>(Var{ param, paramType, nullptr }));
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen)) {
		error(currentLoc, "Expected ')' after parameters.");
	}

	if (!match(TokenType::LeftBrace)) {
		error(currentLoc, "Expected '{' before function body.");
	}

	auto body = parseBlockStatement();

	auto stmt = std::make_unique<FunctionStmt>(
		name,
		type,
		std::move(parameters),
		std::move(body)
	);
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<Statement> Parser::parseReturnStatement() {

	Token startToken = previous(); assert(startToken.type == TokenType::Return);

	std::unique_ptr<Expression> value = nullptr;
	if (!check(TokenType::Semicolon)) {
		value = parseExpression();
	}

	if (!match(TokenType::Semicolon)) {
		error(currentLoc, "Expected ';' after return value.");
	}

	std::unique_ptr<Statement> ret = std::make_unique<febcode::ReturnStmt>(std::move(value));
	ret->location = startToken.loc;
	return ret;
}

std::unique_ptr<febcode::Statement> Parser::parseIfStatement() {

	Token startToken = previous(); assert(startToken.type == TokenType::If);

	if (!match(TokenType::LeftParen)) {
		error(currentLoc, "Expected '(' after 'if'.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::RightParen)) {
		error(currentLoc, "Expected ')' after if condition.");
	}

	auto thenBranch = parseStatement(); // could be block or single statement

	std::unique_ptr<Statement> elseBranch = nullptr;
	if (match(TokenType::Else)) {
		elseBranch = parseStatement();
	}

	auto stmt = std::make_unique<IfStmt>(std::move(condition),
		std::move(thenBranch),
		std::move(elseBranch));
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<febcode::Statement> Parser::parseWhileStatement() {

	Token startToken = previous(); assert(startToken.type == TokenType::While);

	if (!match(TokenType::LeftParen)) {
		error(currentLoc, "Expected '(' after 'while'.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::RightParen)) {
		error(currentLoc, "Expected ')' after while condition.");
	}

	auto body = parseStatement(); // can be block or single statement

	auto stmt = std::make_unique<WhileStmt>(
		std::move(condition),
		std::move(body));
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<febcode::Statement> Parser::parseForStatement() {

	Token startToken = previous(); assert(startToken.type == TokenType::For);

	if (!match(TokenType::LeftParen)) {
		error(currentLoc, "Expected '(' after 'for'.");
	}

	auto init = parseDeclaration(); // initializer (can be var decl, expression stmt, or empty)
	if (!init && !isVarDecl(init) && !isExprStmt(init)) {
		error(currentLoc, "Expected variable declaration, expression statement, or ';' in for loop initializer.");
	}

	auto condition = parseExpression();

	if (!match(TokenType::Semicolon)) {
		error(currentLoc, "Expected ';' after for initializer.");
	}

	auto increment = parseExpression();

	if (!match(TokenType::RightParen)) {
		error(currentLoc, "Expected ')' after for increment.");
	}

	auto body = parseStatement(); // can be block or single statement

	auto stmt = std::make_unique<ForStmt>(
		std::move(init),
		std::move(condition),
		std::move(increment),
		std::move(body));
	stmt->location = startToken.loc;
	return stmt;
}

std::unique_ptr<Expression> Parser::parseAssignment() {
	auto expr = parseOr(); // parse the left-hand side first

	// Check if this is an assignment
	if (match(TokenType::Equal)) {
		SourceLocation opLoc = previous().loc;
		auto value = parseAssignment();
		auto assign = std::make_unique<AssignExpr>(
			std::move(expr),
			std::move(value)
		);
		assign->location = opLoc;
		return assign;
	}

	return expr; // just an equality/expression if no '='
}

std::unique_ptr<Expression> Parser::parseOr()
{
	auto expr = parseAnd();

	while (match(TokenType::OrOr))
	{
		BinaryOp op = tokenToBinaryOp(previous());
		SourceLocation opLoc = previous().loc;
		auto right = parseAnd();
		expr = std::make_unique<BinaryExpr>(
			std::move(expr),
			op,
			std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseAnd()
{
	auto expr = parseEquality();

	while (match(TokenType::AndAnd))
	{
		BinaryOp op = tokenToBinaryOp(previous());
		SourceLocation opLoc = previous().loc;
		auto right = parseEquality();
		expr = std::make_unique<BinaryExpr>(
			std::move(expr),
			op,
			std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseEquality() {
	auto expr = parseComparison();

	while (match(TokenType::EqualEqual) || match(TokenType::NotEqual)) {
		BinaryOp op = tokenToBinaryOp(previous());
		SourceLocation opLoc = previous().loc;
		auto right = parseComparison();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseComparison() {
	auto expr = parseTerm();

	while (match(TokenType::Greater) || match(TokenType::GreaterEqual) ||
		match(TokenType::Less) || match(TokenType::LessEqual)) {
		BinaryOp op = tokenToBinaryOp(previous());
		SourceLocation opLoc = previous().loc;
		auto right = parseTerm();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseTerm() {
	auto expr = parseFactor();

	while (match(TokenType::Plus) || match(TokenType::Minus)) {
		SourceLocation opLoc = previous().loc;
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseFactor();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseFactor() {
	auto expr = parseUnary();

	while (match(TokenType::Star) || match(TokenType::Slash)) {
		SourceLocation opLoc = previous().loc;
		BinaryOp op = tokenToBinaryOp(previous());
		auto right = parseUnary();
		expr = std::make_unique<BinaryExpr>(std::move(expr), op, std::move(right));
		expr->location = opLoc;
	}

	return expr;
}

std::unique_ptr<Expression> Parser::parseUnary() {
	if (match(TokenType::Minus) || match(TokenType::Not)) {
		UnaryOp op = tokenToUnaryOp(previous());
		SourceLocation opLoc = previous().loc;

		auto right = parseUnary();

		// absorbe negative sign for scalars.
		if (op == UnaryOp::Negate) {

			int n;
			if (isInt(right.get(), n))
			{
				auto expr = std::make_unique<LiteralExpr>(-n);
				expr->location = opLoc;
				return expr;
			}

			double d;
			if (isDouble(right.get(), d))
			{
				auto expr = std::make_unique<LiteralExpr>(-d);
				expr->location = opLoc;
				return expr;
			}
		}

		auto expr = std::make_unique<UnaryExpr>(op, std::move(right));
		expr->location = opLoc;
		return expr;
	}

	return parseExponent();
}

std::unique_ptr<Expression> Parser::parseExponent() {
	auto left = parseCall();

	if (match(TokenType::Exponent)) {
		BinaryOp op = tokenToBinaryOp(previous());
		SourceLocation opLoc = previous().loc;

		auto right = parseUnary();
		auto expr = std::make_unique<BinaryExpr>(
			std::move(left), op, std::move(right));
		expr->location = opLoc;
		return expr;
	}

	return left;
}

std::unique_ptr<Expression> Parser::parseCall() {

	if (isType())
	{
		std::string typeName = lexeme(previous());
		Type type = prg.types.getType(typeName);
		if (type == nullptr)
			error(currentLoc, "Unknown type name in constructor expression.");

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
				error(currentLoc, "Expected property name after '.'.");
			}

			expr = std::make_unique<MemberExpr>(
				std::move(expr),
				std::string(previous().start, previous().length)
			);
		}
		else if (match(TokenType::LeftBrack)) {
			auto index = parseExpression();
			if (!match(TokenType::RightBrack))
				error(currentLoc, "Expect ']' after index.");

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
		error(currentLoc, "Expected '(' after type name in constructor expression.");

	if (!check(TokenType::RightParen)) {
		do {
			arguments.push_back(parseExpression());
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen))
		error(currentLoc, "Expected ')' after arguments.");

	return std::make_unique<ConstructorExpr>(type, std::move(arguments));
}

std::unique_ptr<Expression> Parser::finishCall(std::unique_ptr<Expression> callee) 
{
	auto* varExpr = dynamic_cast<VariableExpr*>(callee.get());
	if (!varExpr) {
		error(currentLoc, "Can only call functions by name.");
	}

	std::vector<std::unique_ptr<Expression>> arguments;

	if (!check(TokenType::RightParen)) {
		do {
			arguments.push_back(parseExpression());
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen))
		error(currentLoc, "Expected ')' after arguments.");

	return std::make_unique<CallExpr>(
		varExpr->name,
		std::move(arguments)
	);
}

std::unique_ptr<Expression> Parser::parsePrimary() 
{
	std::unique_ptr<Expression> expr;

	if (match(TokenType::Integer)) {
		int value = std::stoi(std::string(previous().start, previous().length));
		expr = std::make_unique<LiteralExpr>(value);
		expr->location = previous().loc;
		return expr;
	}

	if (match(TokenType::Double)) {
		double value = std::stod(std::string(previous().start, previous().length));
		expr = std::make_unique<LiteralExpr>(value);
		expr->location = previous().loc;
		return expr;
	}

	if (match(TokenType::True))
	{
		expr = std::make_unique<LiteralExpr>(true);
		expr->location = previous().loc;
		return expr;
	}
	
	if (match(TokenType::False))
	{
		expr = std::make_unique<LiteralExpr>(false);
		expr->location = previous().loc;
		return expr;
	}

	if (match(TokenType::Identifier)) {
		std::string name(previous().start, previous().length);

		// check some special built-in constants like "PI"
		if (name == "PI") {
			expr = std::make_unique<LiteralExpr>(3.14159265358979323846);
			expr->location = previous().loc;
			return expr;
		}

		// otherwise, it's a variable
		auto expr = std::make_unique<VariableExpr>(
			std::string(previous().start, previous().length));
		expr->location = previous().loc;
		return expr;
	}

	if (match(TokenType::LeftParen)) {
		auto expr = parseExpression();
		if (!match(TokenType::RightParen)) {
			error(currentLoc, "Expected ')' after expression.");
		}
		return expr;
	}

	if (match(TokenType::LeftBrace)) 
	{
		std::vector<std::unique_ptr<Expression>> elements;

		SourceLocation opLoc = previous().loc;

		if (!check(TokenType::RightBrace)) {
			do {
				elements.push_back(parseExpression());
			} while (match(TokenType::Comma));
		}

		if (!match(TokenType::RightBrace))
			error(currentLoc, "Expected '}' after initializer list.");

		auto expr = std::make_unique<InitExpr>(std::move(elements));
		expr->location = opLoc;
		return expr;
	}

	if (match(TokenType::Error))
	{
		error(previous().loc, "Unexpected token.");
	}

	error(currentLoc, "Expected expression.");
}

//---------------------------------------------------------------------
static int l = 0;
static void printTabs(std::ostream& os)
{
	for (int i = 0; i < l; ++i) os << "    ";
}

static std::ostream& operator << (std::ostream& o, const SourceLocation& loc)
{
	return o << "Line " << loc.line << ", Column " << loc.column;
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

static void printExpr(std::ostream& os, const Expression* e);

static void printLiteralExpr(std::ostream& os, const LiteralExpr* e)
{
	os << ValueToString(e->value);
}

static void printVariableExpr(std::ostream& os, const VariableExpr* e)
{
	os << e->name;
}

static void printMemberExpr(std::ostream& os, const MemberExpr* e)
{
	os << "MemberExpr {\n"; l++;
	printTabs(os); os << "object: "; printExpr(os, e->object.get()); os << ",\n";
	printTabs(os); os << "property: " << e->property << "\n";
	l--; printTabs(os); os << "}";
}

static void printAssignmentExpr(std::ostream& os, const AssignExpr* e)
{
	os << "AssignExpr {\n"; l++;
	printTabs(os); os << "target: "; printExpr(os, e->target.get()); os << ",\n";
	printTabs(os); os << "value: "; printExpr(os, e->value.get()); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printUnaryExpr(std::ostream& os, const UnaryExpr* e)
{
	os << "UnaryExpr {\n"; l++;
	printTabs(os); os << "op: " << e->op << ",\n";
	printTabs(os); os << "right: "; printExpr(os, e->right.get()); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printBinaryExpr(std::ostream& os, const BinaryExpr* e)
{
	os << "BinaryExpr {\n"; l++;
	printTabs(os); os << "op: " << e->op << ",\n";
	printTabs(os); os << "left: "; printExpr(os, e->left.get()); os << ",\n";
	printTabs(os); os << "right: "; printExpr(os, e->right.get()); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printExprList(std::ostream& os, const std::vector<std::unique_ptr<Expression>>& a)
{
	os << "[\n"; l++;
	for (size_t i = 0; i < a.size(); ++i)
	{
		printTabs(os); printExpr(os, a[i].get());
		if (i != a.size() - 1) os << ",\n";
		else os << "\n";
	}
	l--; printTabs(os); os << "]";
}

static void printCallExpr(std::ostream& os, const CallExpr* e)
{
	os << "CallExpr {\n"; l++;
	printTabs(os); os << "callee: " << e->name << ",\n";
	printTabs(os); os << "args: "; printExprList(os, e->arguments); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printInitializerExpr(std::ostream& os, const InitExpr* e)
{
	os << "InitializerExpr {\n"; l++;
	printTabs(os); os << "elements: "; printExprList(os, e->elements); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printIndexExpr(std::ostream& os, const IndexExpr* e)
{
	os << "IndexExpr {\n"; l++;
	printTabs(os); os << "object: "; printExpr(os, e->object.get()); os << ",\n";
	printTabs(os); os << "index: "; printExpr(os, e->index.get()); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printConstructorExpr(std::ostream& os, const ConstructorExpr* e)
{
	os << "ConstructorExpr {\n"; l++;
	printTabs(os); os << "type: " << TypeToString(e->valType) << ",\n";
	printTabs(os); os << "args: "; printExprList(os, e->args); os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printExpr(std::ostream& os, const Expression* e)
{
	if      (auto l = dynamic_cast<const LiteralExpr*      >(e)) printLiteralExpr      (os, l);
	else if (auto v = dynamic_cast<const VariableExpr*     >(e)) printVariableExpr     (os, v);
	else if (auto m = dynamic_cast<const MemberExpr*       >(e)) printMemberExpr       (os, m);
	else if (auto a = dynamic_cast<const AssignExpr*       >(e)) printAssignmentExpr   (os, a);
	else if (auto u = dynamic_cast<const UnaryExpr*        >(e)) printUnaryExpr        (os, u);
	else if (auto b = dynamic_cast<const BinaryExpr*       >(e)) printBinaryExpr       (os, b);
	else if (auto c = dynamic_cast<const CallExpr*         >(e)) printCallExpr         (os, c);
	else if (auto c = dynamic_cast<const InitExpr*         >(e)) printInitializerExpr  (os, c);
	else if (auto c = dynamic_cast<const IndexExpr*        >(e)) printIndexExpr        (os, c);
	else if (auto c = dynamic_cast<const ConstructorExpr*  >(e)) printConstructorExpr  (os, c);
	else if (e == nullptr)
		os << "null";
	else
		os << "(Unknown Expression)";
}

static void printExpressionStmt(std::ostream& os, const ExpressionStmt* s)
{
	os << "ExpressionStmt: {\n"; l++;
	printTabs(os); os << "expr: "; printExpr(os, s->expr.get());
	os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printStatement(std::ostream& os, const Statement* stmt);

static std::ostream& operator << (std::ostream& o, TypeKind type)
{
	o << febcode::TypeKindToString(type);
	return o;
}

static void printVarDeclStmt(std::ostream& os, const VarDeclStmt* s)
{
	os << "VarDeclStmt: {\n"; l++;
	printTabs(os); os << "type: " << TypeToString(s->baseType) << ",\n";
	printTabs(os); os << "vars: [\n"; l++;
	for (const auto& var : s->vars) {
		printTabs(os); os << "{ name: " << var->name;
		if (isArrayType(var->type))
		{
			os << ", size: ";
			Type type = var->type;
			while (isArrayType(type))
			{
				os << "[" << type->arraySize <<	"]";
				type = type->elementType;	
			}
		}
		os << ", initializer: ";
		printExpr(os, var->initializer.get());
		os << " },\n";
	}
	l--; printTabs(os); os << "],\n";
	l--;
	printTabs(os); os << "}";
}

static void printReturnStmt(std::ostream& os, const ReturnStmt* s)
{
	os << "ReturnStmt: {\n"; l++;
	printTabs(os); os << "value: "; printExpr(os, s->value.get());
	os << "\n";
	l--;
	printTabs(os); os << "}";
}

static void printBlockStmt(std::ostream& os, const BlockStmt* s)
{
	os << "BlockStmt: [\n"; l++;
	size_t n = s->statements.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& stmt = s->statements[i];
		printStatement(os, stmt.get());
		if (i != n - 1) os << ",\n";
	}	os << "\n";
	l--;
	printTabs(os); os << "]";
}

static void printIfStmt(std::ostream& os, const IfStmt* s)
{
	os << "IfStmt: {\n"; l++;
	printTabs(os); os << "condition: "; printExpr(os, s->condition.get()); os << ",\n";
	printTabs(os); os << "thenBranch: {\n"; l++; printStatement(os, s->thenBranch.get());
	os << "\n";
	l--; printTabs(os); os << "}\n";
	if (s->elseBranch)
	{
		printTabs(os); os << "elseBranch: {\n"; l++;
		printStatement(os, s->elseBranch.get());
		os << "\n";
		l--; printTabs(os); os << "}\n";
	}
	l--;
	printTabs(os); os << "}";
}

static void printWhileStmt(std::ostream& os, const WhileStmt* s)
{
	os << "WhileStmt: {\n"; l++;
	printTabs(os); os << "condition: "; printExpr(os, s->condition.get()); os << ",\n";
	printTabs(os); os << "body: {\n"; l++; printStatement(os, s->body.get());
	os << "\n";
	l--; printTabs(os); os << "}\n";
}

static void printForStmt(std::ostream& os, const ForStmt* s)
{
	os << "ForStmt: {\n"; l++;
	printTabs(os); os << "initializer: "; printStatement(os, s->initializer.get()); os << ",\n";
	printTabs(os); os << "condition: "; printExpr(os, s->condition.get()); os << ",\n";
	printTabs(os); os << "increment: "; printExpr(os, s->increment.get()); os << ",\n";
	printTabs(os); os << "body: {\n"; l++; printStatement(os, s->body.get());
	os << "\n";
	l--; printTabs(os); os << "}\n";
	l--; printTabs(os); os << "}";
}

static void printFunctionStmt(std::ostream& os, const FunctionStmt* s)
{
	os << "FunctionStmt: {\n"; l++;
	printTabs(os); os << "type: " << TypeToString(s->returnType) << ",\n";
	printTabs(os); os << "name: " << s->name << ",\n";
	printTabs(os); os << "params: [\n"; l++;
	for (const auto& param : s->params)
	{
		printTabs(os); os << "{ type: " << TypeToString(param->type) << ", name: " << param->name << " },\n";
	}
	l--; printTabs(os); os << "],\n";
	printTabs(os); os << "body: {\n"; l++; printStatement(os, s->body.get());
	os << "\n";
	l--; printTabs(os); os << "}\n";
	l--; printTabs(os); os << "}";
}

static void printStructStmt(std::ostream& os, const StructStmt* s)
{
	os << "StructStmt: {\n"; l++;
	printTabs(os); os << "name: " << s->name << ",\n";
	printTabs(os); os << "fields: [\n"; l++;
	for (const auto& field : s->fields)
	{
		printTabs(os); os << "{ name: " << field.second << ", type: " << TypeToString(field.first) << " },\n";
	}
	l--; printTabs(os); os << "]\n";
	l--; printTabs(os); os << "}";
}

static void printStatement(std::ostream& os, const Statement* stmt)
{
	printTabs(os); os << stmt->location; os << ": ";
	if      (auto e = dynamic_cast<const ExpressionStmt*>(stmt)) printExpressionStmt(os, e);
	else if (auto v = dynamic_cast<const VarDeclStmt*   >(stmt)) printVarDeclStmt   (os, v);
	else if (auto r = dynamic_cast<const ReturnStmt*    >(stmt)) printReturnStmt    (os, r);
	else if (auto b = dynamic_cast<const BlockStmt*     >(stmt)) printBlockStmt     (os, b);
	else if (auto i = dynamic_cast<const IfStmt*        >(stmt)) printIfStmt        (os, i);
	else if (auto w = dynamic_cast<const WhileStmt*     >(stmt)) printWhileStmt     (os, w);
	else if (auto l = dynamic_cast<const ForStmt*       >(stmt)) printForStmt       (os, l);
	else if (auto f = dynamic_cast<const FunctionStmt*  >(stmt)) printFunctionStmt  (os, f);
	else if (auto s = dynamic_cast<const StructStmt*    >(stmt)) printStructStmt    (os, s);
	else
		os << "(Unknown Statement)";
}

void febcode::printAST(std::ostream& os, const AST& ast)
{
	l = 0;
	size_t n = ast.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto stmt = ast[i];
		printStatement(os, stmt);
		if (i != n - 1) os << ",\n";
	}

	os << std::endl;
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
	if (stmt.input)
		os << "in ";
	os << TypeToString(stmt.baseType) << " ";
	for (size_t i = 0; i < n; ++i)
	{
		const auto& var = stmt.vars[i];
		os << var->name;

		Type type = var->type;
		while (isArrayType(type))
		{
			os << "[" << type->arraySize << "]";
			type = type->elementType;
		}
		if (var->initializer)
		{
			os << " = ";
			prettyPrintExpression(os, *var->initializer);
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
		os << TypeToString(param->type) << " " << param->name;
		if (i != n - 1) os << ", ";
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
		os << "(Unknown Statement)";
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
