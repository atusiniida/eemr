//
// MyException.h
//
#ifndef MyException_H
#define MyException_H

#include <iostream>
#include <string>
#include <exception>
using namespace std;

namespace utility{

class MyException : public exception {
public:
	MyException();
	MyException(string msg);
	virtual ~MyException() throw();

private:
	static long long serialVersionUID;// = 1L;
	string errMsg;

public:
	string getMessage();
	string what_message(string);
};
} // namespace utility
#endif // MyException_H
