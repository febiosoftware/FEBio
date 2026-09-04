#include "scanner.h"

using namespace febcode;

const char* febcode::TokenTypeToString(TokenType type)
{
	switch (type)
	{
	case TokenType::LeftParen   : return "LeftParen";
	case TokenType::RightParen  : return "RightParen";
	case TokenType::LeftBrace   : return "LeftBrace";
	case TokenType::RightBrace  : return "RightBrace";
	case TokenType::LeftBrack   : return "LeftBrack";
	case TokenType::RightBrack  : return "RightBrack";
	case TokenType::Plus        : return "Plus";
	case TokenType::Minus       : return "Minus";
	case TokenType::Star        : return "Star";
	case TokenType::Slash       : return "Slash";
	case TokenType::Exponent    : return "Exponent";
	case TokenType::Less        : return "Less";
	case TokenType::Greater     : return "Greater";
	case TokenType::Equal       : return "Equal";
	case TokenType::Not         : return "Not";
	case TokenType::Semicolon   : return "Semicolon";
	case TokenType::Colon       : return "Colon";
	case TokenType::Comma       : return "Comma";
	case TokenType::Dot         : return "Dot";
	case TokenType::NotEqual    : return "NotEqual";
	case TokenType::EqualEqual  : return "EqualEqual";
	case TokenType::GreaterEqual: return "GreaterEqual";
	case TokenType::LessEqual   : return "LessEqual";
	case TokenType::AndAnd      : return "AndAnd";
	case TokenType::OrOr        : return "OrOr";
	case TokenType::Identifier  : return "Identifier";
	case TokenType::Integer     : return "Integer";
	case TokenType::Double      : return "Double";
	case TokenType::True        : return "True";
	case TokenType::False       : return "False";
	case TokenType::Return      : return "Return";
	case TokenType::If          : return "If";
	case TokenType::Else        : return "Else";
	case TokenType::While       : return "While";
	case TokenType::For         : return "For";
	case TokenType::Struct      : return "Struct";
	case TokenType::Input       : return "Input";
	case TokenType::Type        : return "Type";
	case TokenType::EndOfFile   : return "EndOfFile";
	case TokenType::Error       : return "Error";
	default:
		return "<unknown token>";
	}
}
