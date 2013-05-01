
#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include<string>
using namespace std;
/**
 * Class error of filename
 */
class WrongFile
{
	string name; ///> name of file
public:
	/**
	 * Constructor
	 * @param name File's name
	 */
	WrongFile(string name) {
		this->name = name;
	}
	/**
	 *
	 * @return file's name
	 */
	string getName() const { return name;}
};

/**
 * Class error with edge
 */
class InvalidEdge
{
	int origin; ///> edge's origin
	int destination; ///> edge's destination
public:
	/**
	 * Constructor
	 * @param o origin node
	 * @param d destination node
	 */
	InvalidEdge(int o, int d) :	origin(o), destination(d) {
	}
	/**
	 *
	 * @return edge's destination
	 */
	int getDestination() const {
		return destination;
	}
	/**
	 *
	 * @return edge's origin
	 */
	int getOrigin() const {
		return origin;
	}
};


#endif /* EXCEPTIONS_H_ */
