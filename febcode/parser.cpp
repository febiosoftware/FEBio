#include "parser.h"
#include <iostream>
using namespace febcode;

void Parser::synchronize() {
	advance(); // skip the erroneous token

	while (!isAtEnd()) {
		// Stop if the previous token was a semicolon (likely end of statement)
		if (previous().type == TokenType::Semicolon) return;

		// Check if the current token could be the start of a new statement/declaration
		switch (peek().type) {
		case TokenType::Type:
		case TokenType::Return:
		case TokenType::If:
		case TokenType::While:
		case TokenType::Struct:
			return;
		default:
			break;
		}

		advance(); // keep skipping
	}
}

std::unique_ptr<febcode::Statement> Parser::parseDeclaration() {
	try {
		if (isType())
		{
			Type type = prg.types.getType(lexeme(previous()));

			if (match(TokenType::Identifier)) {
				std::string name = lexeme(previous());

				if (match(TokenType::LeftParen)) {
					// function declaration
					return parseFunctionDeclaration(type, name);
				}

				return parseVarDeclaration(type, name);
			}
			else
				throw std::runtime_error("Expected identifier after type.");
		}
		if (match(TokenType::Struct)) return parseStructDeclaration();
		return parseStatement();
	}
	catch (const std::runtime_error& e) {
		std::cerr << "ERROR: Parse error in parseDeclaration: " << e.what() << "\n";
		synchronize();  // skip tokens until a safe point
		return nullptr; // skip malformed declaration
	}
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
	if (match(TokenType::Return)) return parseReturnStatement();
	if (match(TokenType::LeftBrace)) return parseBlockStatement();
	if (match(TokenType::If)) return parseIfStatement();
	if (match(TokenType::While)) return parseWhileStatement();

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
	if (type == prg.types.getVoidType())
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

		vars.push_back({varType, varName, std::move(initializer)});

		if (!match(TokenType::Comma)) break;

		if (!match(TokenType::Identifier))
			throw std::runtime_error("Identifier expected after comma.");

		varName = std::string(lexeme(previous()));
	}

	if (!match(TokenType::Semicolon)) {
		throw std::runtime_error("Expected ';' after variable declaration.");
	}

	return std::make_unique<VarDeclStmt>(vars);
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

std::unique_ptr<Expression> Parser::parseAssignment() {
	auto expr = parseOr(); // parse the left-hand side first

	// Check if this is an assignment
	if (match(TokenType::Equal)) {
		const Token& equalsToken = previous();

		auto value = parseAssignment();

		// Variable assignment
		if (auto var = dynamic_cast<VariableExpr*>(expr.get()))
		{
			return std::make_unique<AssignExpr>(
				var->name,
				std::move(value));
		}

		// Member assignment
		if (auto member = dynamic_cast<MemberExpr*>(expr.get()))
		{
			return std::make_unique<SetMemberExpr>(
				std::move(member->object),
				member->property,
				std::move(value)
			);
		}

		// index assignment
		if (auto index = dynamic_cast<IndexExpr*>(expr.get()))
		{
			return std::make_unique<SetIndexExpr>(
				std::move(index->object),
				std::move(index->index),
				std::move(value)
			);
		}

		// Error: invalid assignment target
		throw std::runtime_error("Invalid assignment target.");
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

std::unique_ptr<Expression> Parser::finishCall(std::unique_ptr<Expression> callee) {
	std::vector<std::unique_ptr<Expression>> arguments;

	if (!check(TokenType::RightParen)) {
		do {
			arguments.push_back(parseExpression());
		} while (match(TokenType::Comma));
	}

	if (!match(TokenType::RightParen))
		throw std::runtime_error("Expected ')' after arguments.");

	return std::make_unique<CallExpr>(
		std::move(callee),
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
		return std::make_unique<VariableExpr>(
			std::string(previous().start, previous().length));
	}

	if (match(TokenType::String)) {
		return std::make_unique<LiteralExpr>(std::string(previous().start, previous().length));
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

		return std::make_unique<InitializerExpr>(std::move(elements));
	}

	throw std::runtime_error("Expected expression.");
}

//---------------------------------------------------------------------
static int l = 0;
static void printTabs()
{
	for (int i = 0; i < l; ++i) std::cout << "    ";
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

static void printExpr(const Program& prg, const Expression* e);

static void printLiteralExpr(const Program& prg, const LiteralExpr* e)
{
	std::cout << ValueToString(e->value);
}

static void printMemberExpr(const Program& prg, const MemberExpr* e)
{
	std::cout << "MemberExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(prg, e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "property: " << e->property << "\n";
	l--; printTabs(); std::cout << "}";
}

static void printSetMemberExpr(const Program& prg, const SetMemberExpr* e)
{
	std::cout << "SetMemberExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(prg, e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "property: " << e->property << "\n";
	printTabs(); std::cout << "value: "; printExpr(prg, e->value.get()); std::cout << "\n";
	l--; printTabs(); std::cout << "}";
}


static void printVariableExpr(const Program& prg, const VariableExpr* e)
{
	std::cout << e->name;
}

static void printAssignmentExpr(const Program& prg, const AssignExpr* e)
{
	std::cout << "{ name: " << e->name << ", value: "; 
	printExpr(prg, e->value.get()); std::cout << " }";
}

static void printUnaryExpr(const Program& prg, const UnaryExpr* e)
{
	std::cout << "UnaryExpr {\n"; l++;
	printTabs(); std::cout << "op: " << e->op << ",\n";
	printTabs(); std::cout << "right: "; printExpr(prg, e->right.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printBinaryExpr(const Program& prg, const BinaryExpr* e)
{
	std::cout << "BinaryExpr {\n"; l++;
	printTabs(); std::cout << "op: " << e->op << ",\n";
	printTabs(); std::cout << "left: "; printExpr(prg, e->left.get()); std::cout << ",\n";
	printTabs(); std::cout << "right: "; printExpr(prg, e->right.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printExprList(const Program& prg, const std::vector<std::unique_ptr<Expression>>& a)
{
	std::cout << "[\n"; l++;
	for (size_t i = 0; i < a.size(); ++i)
	{
		printTabs(); printExpr(prg, a[i].get());
		if (i != a.size() - 1) std::cout << ",\n";
		else std::cout << "\n";
	}
	l--; printTabs(); std::cout << "]";
}

static void printCallExpr(const Program& prg, const CallExpr* e)
{
	std::cout << "CallExpr {\n"; l++;
	printTabs(); std::cout << "callee: "; printExpr(prg, e->callee.get()); std::cout << ",\n";
	printTabs(); std::cout << "args: "; printExprList(prg, e->arguments); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printInitializerExpr(const Program& prg, const InitializerExpr* e)
{
	std::cout << "InitializerExpr {\n"; l++;
	printTabs(); std::cout << "elements: "; printExprList(prg, e->elements); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printSetIndexExpr(const Program& prg, const SetIndexExpr* e)
{
	std::cout << "SetIndexExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(prg, e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "index: "; printExpr(prg, e->index.get()); std::cout << ",\n";
	printTabs(); std::cout << "value: "; printExpr(prg, e->value.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printIndexExpr(const Program& prg, const IndexExpr* e)
{
	std::cout << "IndexExpr {\n"; l++;
	printTabs(); std::cout << "object: "; printExpr(prg, e->object.get()); std::cout << ",\n";
	printTabs(); std::cout << "index: "; printExpr(prg, e->index.get()); std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printExpr(const Program& prg, const Expression* e)
{
	if      (auto l = dynamic_cast<const LiteralExpr*      >(e)) printLiteralExpr      (prg, l);
	else if (auto m = dynamic_cast<const MemberExpr*       >(e)) printMemberExpr       (prg, m);
	else if (auto s = dynamic_cast<const SetMemberExpr*    >(e)) printSetMemberExpr    (prg, s);
	else if (auto v = dynamic_cast<const VariableExpr*     >(e)) printVariableExpr     (prg, v);
	else if (auto a = dynamic_cast<const AssignExpr*       >(e)) printAssignmentExpr   (prg, a);
	else if (auto u = dynamic_cast<const UnaryExpr*        >(e)) printUnaryExpr        (prg, u);
	else if (auto b = dynamic_cast<const BinaryExpr*       >(e)) printBinaryExpr       (prg, b);
	else if (auto c = dynamic_cast<const CallExpr*         >(e)) printCallExpr         (prg, c);
	else if (auto c = dynamic_cast<const InitializerExpr*  >(e)) printInitializerExpr  (prg, c);
	else if (auto c = dynamic_cast<const SetIndexExpr*     >(e)) printSetIndexExpr     (prg, c);
	else if (auto c = dynamic_cast<const IndexExpr*        >(e)) printIndexExpr        (prg, c);
	else if (e == nullptr)
		std::cout << "null";
	else
		std::cout << "(Unknown Expression)";
}

static void printExpressionStmt(const Program& prg, const ExpressionStmt* s)
{
	std::cout << "ExpressionStmt: {\n"; l++;
	printTabs(); std::cout << "expr: "; printExpr(prg, s->expr.get());
	std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printStatement(const Program& prg, const Statement* stmt);

static std::ostream& operator << (std::ostream& o, TypeKind type)
{
	switch (type)
	{
	case TypeKind::Void  : return o << "void";
	case TypeKind::Bool  : return o << "bool";
	case TypeKind::Int   : return o << "int";
	case TypeKind::Double: return o << "double";
	case TypeKind::String: return o << "string";
	default:
		return o << "<unknown type>";
	}
}

static void printVarDeclStmt(const Program& prg, const VarDeclStmt* s)
{
	std::cout << "VarDeclStmt: {\n"; l++;
	printTabs(); std::cout << "vars: [\n"; l++;
	for (const auto& var : s->vars) {
		printTabs(); std::cout << "{ type: " << TypeToString(var.type) << ", name: " << var.name << ", initializer: ";
		printExpr(prg, var.initializer.get());
		std::cout << " },\n";
	}
	l--; printTabs(); std::cout << "],\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printReturnStmt(const Program& prg, const ReturnStmt* s)
{
	std::cout << "ReturnStmt: {\n"; l++;
	printTabs(); std::cout << "value: "; printExpr(prg, s->value.get());
	std::cout << "\n";
	l--;
	printTabs(); std::cout << "}";
}

static void printBlockStmt(const Program& prg, const BlockStmt* s)
{
	std::cout << "BlockStmt: [\n"; l++;
	size_t n = s->statements.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& stmt = s->statements[i];
		printStatement(prg, stmt.get());
		if (i != n - 1) std::cout << ",\n";
	}	std::cout << "\n";
	l--;
	printTabs(); std::cout << "]";
}

static void printIfStmt(const Program& prg, const IfStmt* s)
{
	std::cout << "IfStmt: {\n"; l++;
	printTabs(); std::cout << "condition: "; printExpr(prg, s->condition.get()); std::cout << ",\n";
	printTabs(); std::cout << "thenBranch: {\n"; l++; printStatement(prg, s->thenBranch.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
	if (s->elseBranch)
	{
		printTabs(); std::cout << "elseBranch: {\n"; l++;
		printStatement(prg, s->elseBranch.get());
		std::cout << "\n";
		l--; printTabs(); std::cout << "}\n";
	}
	l--;
	printTabs(); std::cout << "}";
}

static void printWhileStmt(const Program& prg, const WhileStmt* s)
{
	std::cout << "WhileStmt: {\n"; l++;
	printTabs(); std::cout << "condition: "; printExpr(prg, s->condition.get()); std::cout << ",\n";
	printTabs(); std::cout << "body: {\n"; l++; printStatement(prg, s->body.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
}

static void printFunctionStmt(const Program& prg, const FunctionStmt* s)
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
	printTabs(); std::cout << "body: {\n"; l++; printStatement(prg, s->body.get());
	std::cout << "\n";
	l--; printTabs(); std::cout << "}\n";
	l--; printTabs(); std::cout << "}";
}

static void printStructStmt(const Program& prg, const StructStmt* s)
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

static void printStatement(const Program& prg, const Statement* stmt)
{
	printTabs();
	if      (auto e = dynamic_cast<const ExpressionStmt*>(stmt)) printExpressionStmt(prg, e);
	else if (auto v = dynamic_cast<const VarDeclStmt*   >(stmt)) printVarDeclStmt   (prg, v);
	else if (auto r = dynamic_cast<const ReturnStmt*    >(stmt)) printReturnStmt    (prg, r);
	else if (auto b = dynamic_cast<const BlockStmt*     >(stmt)) printBlockStmt     (prg, b);
	else if (auto i = dynamic_cast<const IfStmt*        >(stmt)) printIfStmt        (prg, i);
	else if (auto w = dynamic_cast<const WhileStmt*     >(stmt)) printWhileStmt     (prg, w);
	else if (auto f = dynamic_cast<const FunctionStmt*  >(stmt)) printFunctionStmt  (prg, f);
	else if (auto s = dynamic_cast<const StructStmt*    >(stmt)) printStructStmt    (prg, s);
	else
		std::cout << "(Unknown Statement)";
}

void febcode::printAST(const Program& prg)
{
	l = 0;
	size_t n = prg.ast.statements.size();
	for (size_t i = 0; i < n; ++i)
	{
		auto& stmt = prg.ast.statements[i];
		printStatement(prg, stmt.get());
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
}
