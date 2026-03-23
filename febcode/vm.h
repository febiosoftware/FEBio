#pragma once
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include "compiler.h"
#include <assert.h>

namespace febcode
{
	class VM
	{
	public:
		enum { MAX_CALL_DEPTH = 8 };

	public:
		VM() : m_program(nullptr) {}

		void setProgram(const Program& program)
		{
			m_program = &program;
			globalStackSize = m_program->globalStackSize;
			m_stack.resize(globalStackSize + program.maxStackSize);
			stackTop = globalStackSize; // stack starts after global region
		}

		VM(const Program& program)
		{
			setProgram(program);
		}

		Value run()
		{
			if (m_program == nullptr) return Value();
			callFunction(0, 0);              // call "main"
			return execute();
		}

		void setDebugMode(bool b) { m_debug = b; }

		bool stackEmpty() const { return stackTop == globalStackSize; }
		size_t stackSize() const { return stackTop - globalStackSize; }

		Value getGlobal(size_t n)
		{
			if (n >= globalStackSize)
				throw std::runtime_error("Invalid global index: " + std::to_string(n));

			const Program::Global& glob = m_program->globals[n];
			int slot = glob.slot;
			switch (glob.type->kind)
			{
			case TypeKind::Bool  : return getBoolAt  (slot);
			case TypeKind::Int   : return getIntAt   (slot);
			case TypeKind::Double: return getDoubleAt(slot);
			case TypeKind::Vec2  : return getVec2At  (slot);
			case TypeKind::Vec3  : return getVec3At  (slot);
			default:
				return m_stack[slot];
			};
		}

		Value getGlobal(const std::string& name)
		{
			auto it = m_program->globalIndices.find(name);
#ifndef NDEBUG
			if (it == m_program->globalIndices.end())
				throw std::runtime_error("Undefined global variable: " + name);
#endif
			return getGlobal(it->second);
		}

		void setGlobal(int i, const Value& v)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;

			if      (isStruct(v)) m_stack[slot] = copyStruct(getStruct(v));
			else if (isArray (v)) m_stack[slot] = copyArray(getArray(v));
			else
			{
				switch (glob.type->kind)
				{
				case TypeKind::Bool  : setBoolAt  (slot, getBool  (v)); break;
				case TypeKind::Int   : setIntAt   (slot, getInt   (v)); break;
				case TypeKind::Double: setDoubleAt(slot, getDouble(v)); break;
				case TypeKind::Vec2  : setVec2At  (slot, getVec2  (v)); break;
				case TypeKind::Vec3  : setVec3At  (slot, getVec3  (v)); break;
					default:
					throw std::runtime_error("Unsupported global variable type for assignment.");
				}
			}
		}

		void setGlobal(int n, std::initializer_list<Value> values)
		{
#ifndef NDEBUG
			if ((n<0) || (n >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(n));
#endif

			const Program::Global& glob = m_program->globals[n];
			int slot = (int)glob.slot;

			Value& v = m_stack[slot];
			if (isArray(v))
			{
				ArrayValue& arr = getArray(v);

				if (arr.size() != values.size())
					throw std::runtime_error("Array initializer has incorrect number of elements.");

				Type elemType = m_program->types.getBuiltinType(*values.begin());
				std::copy(values.begin(), values.end(), arr.elements.begin());
			}
			else if (isStruct(v))
			{
				StructValue& obj = getStruct(v);
				if (obj.fields.size() != values.size())
					throw std::runtime_error("Array initializer has incorrect number of elements.");

				for (size_t i = 0; i < values.size(); ++i)
				{
					Type fieldType = obj.type->fields[i].first;
					Type valType = m_program->types.getBuiltinType(values.begin()[i]);
					if (fieldType != valType)
						throw std::runtime_error("Type mismatch in struct initializer for field: " + obj.type->fields[i].second);
					obj.fields[i] = values.begin()[i];
				}
			}
			else
				throw std::runtime_error("Initializer lists can only be assigned to arrays and structs.");
		}

	private:

		struct CallFrame
		{
			int functionIndex;
			size_t ip;
			size_t base;
		};

		// ===== Execution Loop =====

		Value execute();
	
		// ===== Call Handling =====

		void callFunction(int fnIndex, int argCount);

		// ===== Helpers =====

		CallFrame& currentFrame()
		{
			return m_frames[frameCount - 1];
		}

		uint8_t readByte()
		{
			return m_program->code[currentFrame().ip++];
		}

		uint16_t readUint16()
		{
			uint16_t high = readByte();
			uint16_t low = readByte();
			return (high << 8) | low;
		}

		void push(const Value& v)
		{
#ifndef NDEBUG
			if (stackTop >= m_stack.size())
				throw std::runtime_error("Stack overflow.");
#endif
			m_stack[stackTop++] = v;
		}

		void pushVoid()
		{
			push(Value());
		}

		void pushBool(bool b)
		{
			push(b);
		}

		void pushInt(int n)
		{
			push(n);
		}

		void pushDouble(double d)
		{
			push(d);
		}

		void pushVec2(const vec2& v)
		{
			push(v.x);
			push(v.y);
		}

		void pushVec3(const vec3& v)
		{
			push(v.x);
			push(v.y);
			push(v.z);
		}

		void pushArray(const ArrayValue& arr)
		{
			for (int i = 0; i < arr.size(); ++i)
			{
				const Value& v = arr.elements[i];
				switch (arr.type->elementType->kind)
				{
				case TypeKind::Bool  : pushBool  (v.b); break;
				case TypeKind::Int   : pushInt   (v.i); break;
				case TypeKind::Double: pushDouble(v.d); break;
				case TypeKind::Vec2  : pushVec2  (v.vec2Value); break;
				case TypeKind::Vec3  : pushVec3  (v.vec3Value); break;
				default:
					throw std::runtime_error("Unsupported array element type for push.");
				};
			}
		}

		void pushStruct(const StructValuePtr& arr)
		{
			push(arr);
		}

		const Value& pop()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[--stackTop];
		}

		bool popBool()
		{
			return pop().b;
		}

		int popInt()
		{
			return pop().i;
		}

		double popDouble()
		{
			return pop().d;
		}

		vec2 popVec2()
		{
			double y = pop().d;
			double x = pop().d;
			return vec2(x, y);
		}

		vec3 popVec3()
		{
			double z = pop().d;
			double y = pop().d;
			double x = pop().d;
			return vec3(x, y, z);
		}

		void popValues(size_t count)
		{
			if (stackTop < globalStackSize + count)
				throw std::runtime_error("Stack underflow.");
			stackTop -= count;
		}

		ArrayValue popArray(Type type)
		{
			assert(type->kind == TypeKind::Array);
			ArrayValue arr;
			arr.type = type;
			arr.elements.resize(type->arraySize);
			for (int i = (int)arr.size() - 1; i >= 0; --i)
			{
				switch (type->elementType->kind)
				{
				case TypeKind::Bool  : arr.elements[i] = popBool(); break;
				case TypeKind::Int   : arr.elements[i] = popInt(); break;
				case TypeKind::Double: arr.elements[i] = popDouble(); break;
				case TypeKind::Vec2  : arr.elements[i] = popVec2(); break;
				case TypeKind::Vec3  : arr.elements[i] = popVec3(); break;
				case TypeKind::Array : arr.elements[i] = std::make_shared<ArrayValue>(popArray(type->elementType)); break;
				case TypeKind::Struct: arr.elements[i] = std::make_shared<StructValue>(popStruct()); break;
				default:
					throw std::runtime_error("Unsupported array element type for pop.");
				};
			}
			return arr;
		}

		const StructValue& popStruct()
		{
			return getStruct(pop());
		}

		const Ref& popRef()
		{
			return pop().ref;
		}

		Value& peek()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[stackTop-1];
		}

		bool peekBool()
		{
			return peek().b;
		}

		int peekInt()
		{
			return peek().i;
		}

		double peekDouble()
		{
			return peek().d;
		}

		vec2 peekVec2()
		{
			return vec2(m_stack[stackTop - 2].d, m_stack[stackTop - 1].d);
		}

		vec3 peekVec3()
		{
			return vec3(m_stack[stackTop - 3].d, m_stack[stackTop - 2].d, m_stack[stackTop - 1].d);
		}

		ArrayValue peekArray(Type type)
		{
			assert(type->kind == TypeKind::Array);
			ArrayValue arr;
			arr.type = type;
			arr.elements.resize(type->size());
			int c = (int)(stackTop - type->size());
			for (int i = 0; i < arr.size(); ++i)
			{
				switch (type->elementType->kind)
				{
				case TypeKind::Bool  : arr.elements[i] = m_stack[c].b; break;
				case TypeKind::Int   : arr.elements[i] = m_stack[c].i; break;
				case TypeKind::Double: arr.elements[i] = m_stack[c].d; break;
				case TypeKind::Vec2  : arr.elements[i] = vec2(m_stack[c].d, m_stack[c + 1].d); break;
				case TypeKind::Vec3  : arr.elements[i] = vec3(m_stack[c].d, m_stack[c + 1].d, m_stack[c + 2].d); break;
				default:
					throw std::runtime_error("Unsupported array element type for peek.");
				};
				c += (int)type->elementType->size();
			}
			return arr;
		}

		void setBoolAt(int slot, bool b)
		{
			m_stack[slot] = b;
		}

		void setIntAt(int slot, int n)
		{
			m_stack[slot] = n;
		}

		void setDoubleAt(int slot, double d)
		{
			m_stack[slot] = d;
		}

		void setVec2At(int slot, const vec2& v)
		{
			m_stack[slot  ] = v.x;
			m_stack[slot+1] = v.y;
		}

		void setVec3At(int slot, const vec3& v)
		{
			m_stack[slot    ] = v.x;
			m_stack[slot + 1] = v.y;
			m_stack[slot + 2] = v.z;
		}

		void setArrayAt(int slot, const ArrayValue& arr)
		{
			for (int i = 0; i < arr.size(); ++i)
			{
				const Value& v = arr.elements[i];
				switch (arr.type->elementType->kind)
				{
				case TypeKind::Bool  : setBoolAt  (slot, v.b); break;
				case TypeKind::Int   : setIntAt   (slot, v.i); break;
				case TypeKind::Double: setDoubleAt(slot, v.d); break;
				case TypeKind::Vec2  : setVec2At  (slot, v.vec2Value); break;
				case TypeKind::Vec3  : setVec3At  (slot, v.vec3Value); break;
				default:
					throw std::runtime_error("Unsupported array element type for setArrayAt.");
				};
				slot += (int)arr.type->elementType->size();
			}
		}

		bool getBoolAt(int slot)
		{
			return m_stack[slot].b;
		}

		int getIntAt(int slot)
		{
			return m_stack[slot].i;
		}

		double getDoubleAt(int slot)
		{
			return m_stack[slot].d;
		}

		vec2 getVec2At(int slot)
		{
			return vec2(m_stack[slot].d, m_stack[slot + 1].d);
		}

		vec3 getVec3At(int slot)
		{
			return vec3(m_stack[slot].d, m_stack[slot + 1].d, m_stack[slot + 2].d);
		}

		ArrayValue getArrayAt(int slot, Type type)
		{
			assert(type->kind == TypeKind::Array);
			ArrayValue arr;
			arr.type = type;
			arr.elements.resize(type->arraySize);

			int c = slot;
			for (int i = 0; i < arr.size(); ++i)
			{
				Value& elem = arr.elements[i];
				switch (type->elementType->kind)
				{
				case TypeKind::Bool  : elem = getBoolAt  (c); break;
				case TypeKind::Int   : elem = getIntAt   (c); break;
				case TypeKind::Double: elem = getDoubleAt(c); break;
				case TypeKind::Vec2  : elem = getVec2At  (c); break;
				case TypeKind::Vec3  : elem = getVec3At  (c); break;
				default:
					throw std::runtime_error("Unsupported type in getArrayAt");
				}
				c += (int)type->elementType->size();
			}
			return arr;
		}

		void copy(int dest, int src, int size)
		{
			for (int i = 0; i < size; ++i)
				m_stack[dest + i] = m_stack[src + i];

			if (dest + size > stackTop)
				stackTop = dest + size;
		}

		bool isTruthy(const Value& v)
		{
			switch (v.index)
			{
				case ValueIndex::VOID  : return false;
				case ValueIndex::BOOL  : return getBool(v);
				case ValueIndex::INT   : return getInt(v) != 0;
				case ValueIndex::DOUBLE: return getDouble(v) != 0.0;
				case ValueIndex::ARRAY : return true;
				case ValueIndex::STRUCT: return true;
			}

			return false;
		}

		std::string toString(const Value& v)
		{
			switch (v.index)
			{
			case ValueIndex::VOID  : return "void";
			case ValueIndex::BOOL  : return getBool(v) ? "true" : "false";
			case ValueIndex::INT   : return std::to_string(getInt(v));
			case ValueIndex::DOUBLE: return std::to_string(getDouble(v));
			}
			return "";
		}

	private:
		const Program* m_program;
		size_t globalStackSize = 0;

		std::vector<Value> m_stack;
		size_t stackTop = 0;

		CallFrame m_frames[MAX_CALL_DEPTH];
		size_t frameCount = 0;

		bool m_debug = false;
	};

	Value runScript(const std::string& script, bool initModules = false);
}
