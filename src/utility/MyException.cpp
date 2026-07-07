//
// MyException.cpp
//
#include "MyException.h"
using namespace utility;

MyException::MyException()
{
}
MyException::MyException(string msg)
{
	errMsg = msg;
}

MyException::~MyException() throw ()
{
}

string MyException::getMessage()
{
	return errMsg;
}

string MyException::what_message(string msg)
{
	string str;

	str  = exception::what();
	str += " " + msg;

	return str;
}
