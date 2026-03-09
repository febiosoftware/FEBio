#pragma once
#include <string>
#include <vector>
namespace febcode {

	enum class TokenType {
		// Single-character tokens
		LeftParen, RightParen,
		LeftBrace, RightBrace,
		LeftBrack, RightBrack,
		Plus, Minus, Star, Slash, Exponent,
		Less, Greater,
		Equal, Not,
		Semicolon,
		Colon,
		Comma,
		Dot,

		// Multi-character tokens
		NotEqual, EqualEqual,
		GreaterEqual, LessEqual,
		AndAnd, OrOr,

		// Literals
		Identifier,
		Integer,
		Double,
		True,
		False,
		String,

		// Keywords
		Return,
		If,
		Else,
		While,
		Struct,

		// for type definitions
		Type,

		// Special
		EndOfFile,
		Error,
	};

	struct Token {
		TokenType type;
		const char* start;
		int length;
		int line;
	};

	class Scanner {
	public:
		Scanner(const std::string& source)
			: source(source)
		{
			start = source.c_str();
			current = start;
			line = 1;
		}

		std::vector<Token> scanTokens() {
			std::vector<Token> tokens;
			while (true) {
				Token token = nextToken();
				tokens.push_back(token);
				if (token.type == TokenType::EndOfFile) break;
			}
			return tokens;
		}

		Token nextToken() {
			skipWhitespace();
			start = current;

			if (isAtEnd()) return makeToken(TokenType::EndOfFile);

			char c = advance();

			// multi-character tokens
			switch (c)
			{
			case '!':
				if (peek() == '=') {
					advance();
					return makeToken(TokenType::NotEqual);
				}
				break;
			case '=':
				if (peek() == '=') {
					advance();
					return makeToken(TokenType::EqualEqual);
				}
				break;
			case '<':
				if (peek() == '=') {
					advance();
					return makeToken(TokenType::LessEqual);
				}
				break;
			case '>':
				if (peek() == '=') {
					advance();
					return makeToken(TokenType::GreaterEqual);
				}
				break;
			case '&':
				if (peek() == '&') {
					advance();
					return makeToken(TokenType::AndAnd);
				}
				break;
			case '|':
				if (peek() == '|') {
					advance();
					return makeToken(TokenType::OrOr);
				}
				break;
			case '*':
				if (peek() == '*') {
					advance();
					return makeToken(TokenType::Exponent);
				}
				break;
			case '/':
				if (peek() == '/') {
					// Comment until end of line
					while (peek() != '\n' && !isAtEnd()) advance();
					return nextToken(); // skip comment token
				}
				break;
			}

			// Single-character tokens
			switch (c) {
			case '(': return makeToken(TokenType::LeftParen);
			case ')': return makeToken(TokenType::RightParen);
			case '{': return makeToken(TokenType::LeftBrace);
			case '}': return makeToken(TokenType::RightBrace);
			case '[': return makeToken(TokenType::LeftBrack);
			case ']': return makeToken(TokenType::RightBrack);
			case '+': return makeToken(TokenType::Plus);
			case '-': return makeToken(TokenType::Minus);
			case '*': return makeToken(TokenType::Star);
			case '/': return makeToken(TokenType::Slash);
			case '=': return makeToken(TokenType::Equal);
			case ';': return makeToken(TokenType::Semicolon);
			case ':': return makeToken(TokenType::Colon);
			case ',': return makeToken(TokenType::Comma);
			case '.': return makeToken(TokenType::Dot);
			case '<': return makeToken(TokenType::Less);
			case '>': return makeToken(TokenType::Greater);
			case '!': return makeToken(TokenType::Not);
			case '\'': return string('\'');
			case '\"': return string('\"');
			}

			if (std::isdigit(c)) return number();
			if (isAlpha(c)) return identifier();

			return errorToken("Unexpected character.");
		}

	private:
		bool isAtEnd() const {
			return *current == '\0';
		}

		char advance() {
			return *current++;
		}

		void skipWhitespace() {
			while (true) {
				char c = *current;
				switch (c) {
				case ' ':
				case '\r':
				case '\t':
					current++;
					break;
				case '\n':
					line++;
					current++;
					break;
				default:
					return;
				}
			}
		}

		Token makeToken(TokenType type) const {
			return Token{ type, start, int(current - start), line };
		}

		Token errorToken(const char* message) const {
			return Token{ TokenType::Error, message, int(std::strlen(message)), line };
		}

		Token number() {
			while (std::isdigit(peek())) advance();

			// fractional part
			if (peek() == '.') {
				advance(); // consume '.'
				while (std::isdigit(peek())) advance();
				if ((peek() == 'e') || (peek() == 'E'))
				{
					advance(); // consume 'e' or 'E'
					if ((peek() == '+') || (peek() == '-')) advance();
					while (std::isdigit(peek())) advance();
				}
				return makeToken(TokenType::Double);
			}
			return makeToken(TokenType::Integer);
		}

		Token identifier() {
			while (isAlphaNumeric(peek())) advance();

			int length = int(current - start);

			// Keyword checks
			if (length == 2 && std::memcmp(start, "if", 2) == 0)
				return makeToken(TokenType::If);

			if (length == 4 && std::memcmp(start, "else", 4) == 0)
				return makeToken(TokenType::Else);

			if (length == 6 && std::memcmp(start, "return", 6) == 0)
				return makeToken(TokenType::Return);

			if (length == 5 && std::memcmp(start, "while", 5) == 0)
				return makeToken(TokenType::While);

			if (length == 6 && std::memcmp(start, "struct", 6) == 0)
				return makeToken(TokenType::Struct);

			if (length == 4 && std::memcmp(start, "true", 4) == 0)
				return makeToken(TokenType::True);

			if (length == 5 && std::memcmp(start, "false", 5) == 0)
				return makeToken(TokenType::False);

			// type checks
			if (length == 4 && std::memcmp(start, "void", 4) == 0)
				return makeToken(TokenType::Type);
			if (length == 4 && std::memcmp(start, "bool", 4) == 0)
				return makeToken(TokenType::Type);
			if (length == 3 && std::memcmp(start, "int", 3) == 0)
				return makeToken(TokenType::Type);
			if (length == 6 && std::memcmp(start, "double", 6) == 0)
				return makeToken(TokenType::Type);
			if (length == 6 && std::memcmp(start, "string", 6) == 0)
				return makeToken(TokenType::Type);
			if (length == 4 && std::memcmp(start, "vec2", 4) == 0)
				return makeToken(TokenType::Type);
			if (length == 4 && std::memcmp(start, "vec3", 4) == 0)
				return makeToken(TokenType::Type);

			// Fallback: regular identifier
			return makeToken(TokenType::Identifier);
		}

		Token string(char quote)
		{
			while (peek() != quote && !isAtEnd()) {
				if (peek() == '\n') line++;
				advance();
			}
			if (isAtEnd()) return errorToken("Unterminated string.");
			advance(); // consume closing '

			// don't include quotes characters
			return Token{ TokenType::String, start + 1, int(current - start) - 2, line };
		}

		char peek() const {
			return *current;
		}

		char peekNext() const {
			if (isAtEnd()) return '\0';
			return current[1];
		}

		static bool isAlpha(char c) {
			return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
		}

		static bool isAlphaNumeric(char c) {
			return isAlpha(c) || std::isdigit(static_cast<unsigned char>(c));
		}

	private:
		const std::string& source;
		const char* start;
		const char* current;
		int line;
	};
}
