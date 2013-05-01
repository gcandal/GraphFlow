#ifndef FORDFULKERSON_H_
#define FORDFULKERSON_H_

#include "graphviewer.h"
#include "Exceptions.h"
#include <iostream>
#include <time.h>
#include <vector>
#include <queue>
#include <fstream>
#include <climits>
#include <string>
#include <cstdio>
#include <sstream>


using namespace std;
/**
 * Ford-Fulkerson Algorithm - time optimized
 */
class FordFulkerson {
	static int const WHITE1 = 0;  ///> needed for breadth-first search - means not visited
	static int const GRAY1 = 1; ///> needed for breadth-first search - means in queue
	static int const BLACK1 = 2; ///> needed for breadth-first search - means visited
	static int const oo = INT_MAX; ///> representation of infinity
	int nodes; ///> number of nodes
	int edges; ///> number of edges
	vector< vector<int> > capacity; ///> capacity matrix
	vector< vector<int> > flow; ///> flow matrix
	vector<int> color; ///> needed for breadth-first search
	vector<int> pred; ///> stores augmenting path
	queue<int> q; ///> queue for breadth-first search
	string file_name; ///> name of file with nodes/edges info
	GraphViewer *g; ///> graph api to show flow/capacity
	int n_edges_g; ///> number of edges in g
	GraphViewer *gr; ///> graph api to show the residual graph
	int n_edges_gr; ///> number of edges in gr
	/**
	 *
	 * @param x number 1
	 * @param y number 2
	 * @return minimum of x and y
	 */
	int min (int x, int y);
	/**
	 * Puts in queue/changes color to gray
	 * @param x node
	 */
	void enqueue(int x);
	/**
	 * Changes color to black
	 * @return first in queue
	 */
	int dequeue();
	/**
	 * Breath first-search - finds path from start to target
	 * @param start origin node
	 * @param target destination node
	 * @return if achieved the target
	 */
	bool bfs(int start, int target);
	/**
	 * Reads and fulfills nodes'/edges' info
	 */
	void read_input_file();
public:
	/**
	 * Constructor
	 * @param file name of file with nodes'/edges' info
	 */
	FordFulkerson(string file);
	/**
	 * Calculates the max flow
	 * @param source origin node
	 * @param sink destination node
	 * @return max flow
	 */
	int max_flow (int source, int sink);
	/**
	 * Prints flow/capacity info and updates the graphs api
	 */
	void print();
	/**
	 * Calculates bottlenecks - say it is bottleneck if capacity out <= 0.7 * capacity in
	 * @return vector with potential nodes that are bottlenecks
	 */
	vector<int> getBottlenecks() const;
	/**
	 * Modifies edge capacity
	 * @param origin origin node
	 * @param destination destination node
	 * @param capac edge capacity
	 */
	void setEdgeCapacity(int origin, int destination, int capac);
	/**
	* Saves the current graph format
	* @param name name of destination file
	*/
	void savetoFile(string name);
	/**
	 * Ford-Fulkerson time optimized menu
	 */
	void menu();
	/**
	 * Menu to edit edge
	 * @return 0 upon success
	 */
	int menuSetEdge();
	/**
	 * Checks if str is a number
	 * @param str str to check
	 * @return
	 */
	bool isNumber(string str);
	/**
	 * Closes graphic api
	 */
	void closeGraphs();

};




#endif /* FORDFULKERSON_H_ */
