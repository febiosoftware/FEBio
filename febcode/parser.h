#pragma once
#include "scanner.h"
#include "ast.h"
#include "program.h"
#include <stdexcept>
#include <ostream>
#include <cstring>

namespace febcode {

	class Parser {
	public:
		Parser(Program& program) : prg(program), current(0), m_ast(nullptr) {}

		void parse(const std::vector<Token>& tokens) {
			m_ast = prg.ast.get();
			this->tokens = tokens;
			current = 0;
			while (!isAtEnd()) {
				m_ast->root.statements.push_back(parseDeclaration());
			}
			m_ast = nullptr;
		}

	private:
		Program& prg;
		std::vector<Token> tokens;
		size_t current;

		AST* m_ast;

		bool isAtEnd() const { return peek().type == TokenType::EndOfFile; }
		const Token& peek() const { return tokens[current]; }
		const Token& previous() const { return tokens[current - 1]; }

		std::string lexeme(const Token& t) const {
			return std::string(t.start, t.length);
		}

		const Token& advance() {
			if (!isAtEnd()) current++;
			return previous();
		}

		void rewind() {
			if (current > 0) current--;
		}

		bool check(TokenType type) const {
			if (isAtEnd()) return false;
			return peek().type == type;
		}

		bool match(TokenType type) {
			if (check(type)) {
				advance();
				return true;
			}
			return false;
		}

		bool isType()
		{
			if (check(TokenType::Type)) return match(TokenType::Type);
			if (check(TokenType::Identifier)) {
				// Check if it's a user-defined struct type
				std::string name = lexeme(peek());
				Type type = prg.types.getStructType(name);
				if (type != nullptr) {
					advance(); // consume the identifier
					return true;
				}
			}
			return false;
		}

		bool checkType()
		{
			if (check(TokenType::Type)) return true;
			if (check(TokenType::Identifier)) {
				// Check if it's a user-defined struct type
				std::string name = lexeme(peek());
				Type type = prg.types.getStructType(name);
				if (type != nullptr) {
					return true;
				}
			}
			return false;
		}

		// --- Parsing functions (placeholders) ---
		std::unique_ptr<Statement> parseDeclaration();
		std::unique_ptr<Statement> parseStatement();
		std::unique_ptr<Statement> parseVarDeclaration(Type type, const std::string& name);
		std::unique_ptr<Statement> parseFunctionDeclaration(Type type, const std::string& name);
		std::unique_ptr<Statement> parseStructDeclaration();
		std::unique_ptr<Statement> parseReturnStatement();
		std::unique_ptr<Statement> parseIfStatement();
		std::unique_ptr<Statement> parseWhileStatement();
		std::unique_ptr<Statement> parseForStatement();
		std::unique_ptr<Statement> parseBlockStatement();
		std::unique_ptr<Statement> parseExpressionStatement();

		std::unique_ptr<Expression> parseExpression() { return parseAssignment(); }
		std::unique_ptr<Expression> parseAssignment();
		std::unique_ptr<Expression> parseOr();
		std::unique_ptr<Expression> parseAnd();
		std::unique_ptr<Expression> parseEquality();
		std::unique_ptr<Expression> parseComparison();
		std::unique_ptr<Expression> parseTerm();
		std::unique_ptr<Expression> parseFactor();
		std::unique_ptr<Expression> parseExponent();
		std::unique_ptr<Expression> parseUnary();
		std::unique_ptr<Expression> parseCall();
		std::unique_ptr<Expression> parseConstructor(Type type);
		std::unique_ptr<Expression> parsePrimary();

		std::unique_ptr<Expression> finishCall(std::unique_ptr<Expression> callee);

		// --- Helper functions to map tokens to enums ---
		BinaryOp tokenToBinaryOp(const Token& t) {
			switch (t.type) {
			case TokenType::Plus: return BinaryOp::Plus;
			case TokenType::Minus: return BinaryOp::Minus;
			case TokenType::Star: return BinaryOp::Multiply;
			case TokenType::Slash: return BinaryOp::Divide;
			case TokenType::Exponent: return BinaryOp::Exponent;
			case TokenType::EqualEqual: return BinaryOp::EqualEqual;
			case TokenType::NotEqual: return BinaryOp::NotEqual;
			case TokenType::Greater: return BinaryOp::Greater;
			case TokenType::GreaterEqual: return BinaryOp::GreaterEqual;
			case TokenType::Less: return BinaryOp::Less;
			case TokenType::LessEqual: return BinaryOp::LessEqual;
			case TokenType::AndAnd: return BinaryOp::AndAnd;
			case TokenType::OrOr: return BinaryOp::OrOr;
			default: throw std::runtime_error("Invalid binary operator token");
			}
		}

		UnaryOp tokenToUnaryOp(const Token& t) {
			switch (t.type) {
			case TokenType::Minus: return UnaryOp::Negate;
			case TokenType::Not: return UnaryOp::Not;
			default: throw std::runtime_error("Invalid unary operator token");
			}
		}
	};

	void printAST(const AST& ast);

	void prettyPrintAST(std::ostream& os, const AST& ast);
	void prettyPrintExpression(std::ostream& os, const Expression& expr);

	void ParseSource(Program& prg, const std::string& source);

} // namespace febcode

