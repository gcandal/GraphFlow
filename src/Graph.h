#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <queue>
#include <list>
#include <time.h>
#include <climits>
#include <cmath>
#include <iostream>
#include <Windows.h>
#include <winbase.h>
#include <string>
#include <sstream>
#include <fstream>
#include <time.h>
#include "graphviewer.h"
#include "Exceptions.h"

using namespace std;

class Edge;
class Graph;

const int INT_INFINITY = INT_MAX;

/*
 * ================================================================================================
 * Class Vertex
 * ================================================================================================
 */

/**
 * Represents a node of a graph
 */
class Vertex {
	int info; ///> node's information
	vector<Edge  > adj; ///> edges which have this node as origin
	bool visited; ///> checks if was visited or not (used for
	int indegree; ///> how many edges have this node as destiny
	double dist; ///> distance to a given node
public:
	Vertex(int in);

	int in; ///> maximum flow that can get in
	int out; ///> maximum flow that can get out
	/**
	 * Adds an edge with this node as origin
	 * @param dest destiny node
	 * @param w edge's maximum capacity
	 */
	void addEdge(Vertex *dest, double w);
	/**
	 * Adds an edge with this node as origin
	 * @param dest destiny node
	 * @param w edge's maximum capacity
	 * @param f edge's flow
	 */
	void addEdge(Vertex *dest, double w, double f);
	/**
	 * Adds an edge with this node as origin
	 * @param dest destiny node
	 * @param w edge's maximum capacity
	 * @param f edge's flow
	 * @param hide hidden edge or not
	 */
	void addEdge(Vertex *dest, double w, double f, bool hide);
	/**
	 * Removes an edge, if it exists
	 * @param d edge's destiny node
	 * @return True if edge exists
	 */
	bool removeEdgeTo(Vertex *d);
	/**
	 * Hides an edge, if it exists
	 * @param d edge's destiny node
	 * @return True if edge exists
	 */
	bool hideEdgeTo(Vertex *d);
	/**
	 * Adds flow a given edge, if it exists
	 * @param d edge's destiny node
	 * @return True if edge exists
	 */
	bool addFlowTo(Vertex *d, const double flow);
	/**
	 * Changes an edge's flow, if it exists
	 * @param d edge's destiny node
	 * @return True if edge exists
	 */
	void updateEdgeFlow(unsigned int index, float f);
	/**
	 * @return Node's information
	 */
	int getInfo() const;
	/**
	 * @return Distance to a given node
	 */
	int getDist() const;
	/**
	 * @return Number of edges with this node as destiny
	 */
	int getIndegree() const;
	/**
	 * @return Set of edges with this node as origin
	 */
	vector<Edge > getAdj() const;
	/**
	 * @return Node before this one in the path which connects it to a given node
	 */
	Vertex* getPath() const;

	Vertex* path; ///> node before this one in the path which connects it to a given node


	friend class Graph;
};


struct vertex_greater_than {
	bool operator()(Vertex * a, Vertex * b) const {
		return a->getDist() > b->getDist();
	}
};



/* ================================================================================================
 * Class Edge
 * ================================================================================================
 */

/*
 * Represents an edge connection two nodes.
 */
class Edge {
	Vertex * dest; ///> destiny node
	Vertex * orig; ///> origin node
	double weight; ///> edge's capacity
	double flow; ///> edge's current flow
	bool hidden; ///> keeps track if the edge is hidden or not (used in Dinic's algorithm)
public:
	/**
	 * @param d destiny node
	 * @param w edge's capacity
	 * @param f edge's flow
	 */
	Edge(Vertex *d, double w, double f=0);
	/**
	 * @param d destiny node
	 * @param w edge's capacity
	 * @param f edge's flow
	 * @param rev hidden edge or not
	 */
	Edge(Vertex *d, double w, double f, bool rev);
	/**
	 * @return Edge's flow
	 */
	double getFlow() const;
	/**
	 * @return Edge's weight
	 */
	double getWeight() const;
	/**
	 * @return Edge's origin node
	 */
	Vertex * getOri();
	/**
	 * @return Edge's desitny node
	 */
	Vertex * getDest();

	friend class Graph;
	friend class Vertex;
};



/* ================================================================================================
 * Class Graph
 * ================================================================================================
 */

class Graph {
	vector<Vertex *> vertexSet; ///> set of the nodes in the graph
	/**
	 * Prints a set of edges
	 * @param e Edges to be printed
	 */
	void printEdges(vector<Edge > e);
	/**
	 * Builds a residual Graph, used on Dinic's and Fold-Fulkerson's algorithms
	 * @return Residual graph
	 */
	Graph getGr();
	/**
	 * @param s source node
	 * @param t target node
	 * @return Flow of a given edge
	 */
	int getEdgeFlow(const int &s, const int &t);
public:
	/**
	 * @return All edges in the graph
	 */
	vector<Edge > getEdges();
	/**
	 * @param in New node's information
	 * @return True on success (false if already existed)
	 */
	bool addVertex(const int &in);
	/**
	 * @param sourc Origin node
	 * @param dest Destiny node
	 * @param w Edge's weight
	 * @param f Edge's flow
	 * @return True on success (false if already existed)
	 */
	bool addEdge(const int &sourc, const int &dest, double w,double f=0);
	/**
	 * @param sourc Origin node
	 * @param dest Destiny node
	 * @param w Edge's weight
	 * @param f Edge's flow
	 * @param hidden If its an hidden edge or not
	 * @return True on success (false if already existed)
	 */
	bool addEdge(const int &sourc, const int &dest, double w,double f,bool hidden);
	/**
	 * @param Information within the node to be removed
	 * @return True on success (false if doesn't exist)
	 */
	bool removeVertex(const int &in);
	/**
	 * @param sourc Edge's origin node information
	 * @param sourc Edge's destiny node information
	 * @return True on success (false if doesn't exist)
	 */
	bool removeEdge(const int &sourc, const int &dest);
	/**
	 * @param sourc Edge's origin node information
	 * @param sourc Edge's destiny node information
	 * @return True on success (false if doesn't exist)
	 */
	bool hideEdge(const int &sourc, const int &dest);
	/**
	 * @param sourc Edge's origin node information
	 * @param sourc Edge's destiny node information
	 * @param flow Flow to be added
	 * @return True on success (false if doesn't exist)
	 */
	bool addFlow(const int &sourc, const int &dest, const double flow);
	/**
	 * Changes the flow of all edges in the graph to 0
	 */
	void resetEdgeFlow();
	/**
	 * @return All nodes in the graph
	 */
	vector<Vertex * > getVertexSet() const;
	/**
	 * @return Number of nodes in the graph
	 */
	int getNumVertex() const;
	/**
	 * @param v information in the node to be returned
	 * @return Node with a given information
	 */
	Vertex* getVertex(const int &v) const;
	/**
	 * @return Nodes in the graph which aren't destiny for any edge
	 */
	vector<Vertex*> getSources() const;
	/**
	 * @param origin Start node
	 * @param dest End node
	 * @return Information of the nodes which make a path between two other nodes
	 */
	vector<int> getPath(const int &origin, const int &dest);
	/**
	 * Calculates the shortest distance between a node in the graph with all others
	 * @param v Information of the node which serves as referential
	 */
	void unweightedShortestPath(const int &v);
	/**
	 * Runs Dinic's algorithm for maximum flow in the graph,
	 * altering it's edges flow for the final one
	 * @param s Source node's information
	 * @param t Sink node's information
	 * @return Graph's maximum flow
	 */
	int Dinic(const int &s, const int &t);
	/**
	 * Runs Ford-Fulkerson's algorithm for maximum flow in the graph,
	 * altering it's edges flow for the final one
	 * @param s Source node's information
	 * @param t Sink node's information
	 * @return Graph's maximum flow
	 */
	int FordFulk(const int &s, const int &t);
	/**
	 * Converts this graph into GraphViewer (graphical API) and show it
	 */
	void plotGraph();
	/**
	 * @param Filename/path.
	 * Reads a graph from a file
	 */
	void readFile(string name);
	/**
	 * @param Filename/path.
	 * Writes the graph to a file
	 */
	void writeFile(string name);
	/**
	 * @param str String to be checked
	 * @return True if all chars are digits
	 */
	bool isNumber(string str);

	int menuEdges();
	void menuFord();
	void menuDinic();

};


#endif
