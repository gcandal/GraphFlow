#include "Graph.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <ctime>
#include "graphviewer.h"

using namespace std;

template <class T>
void plotGraph(Graph<T> graph) {
	stringstream ss;
	GraphViewer *gv = new GraphViewer(600, 600, true);
	gv->createWindow(600, 600);

	vector<Edge<int> > e=graph.getEdges();
	vector<Vertex<int> *> v=graph.getVertexSet();

	for(unsigned int i=0; i<v.size(); i++) {
		gv->addNode(v[i]->getInfo());
		if(v[i]->in>v[i]->out)
			gv->setVertexColor(v[i]->getInfo(),"red");
	}

	for(unsigned int i=0; i<e.size(); i++) {
		ss.str("");
		ss << e[i].getFlow() << "/" << e[i].getWeight();
		gv->addEdge(i,e[i].getOri()->getInfo(),e[i].getDest()->getInfo(),EdgeType::DIRECTED);
		gv->setEdgeLabel(i,ss.str());
		gv->setEdgeThickness(i,e[i].getFlow()*0.5);

	}

	gv->rearrange();
}

Graph<int> createTestFlowGraph()
		{
	Graph<int> myGraph;

	for(int i = 1; i < 7; i++)
		myGraph.addVertex(i);

	myGraph.addEdge(1, 2, 3, 0);
	myGraph.addEdge(1, 3, 2, 0);
	myGraph.addEdge(2, 5, 4, 0);
	myGraph.addEdge(2, 4, 3, 0);
	myGraph.addEdge(2, 3, 1, 0);
	myGraph.addEdge(3, 5, 2, 0);
	myGraph.addEdge(4, 6, 2, 0);
	myGraph.addEdge(5, 6, 3, 0);

	return myGraph;
		}

Graph<int> readFile(string name) {
	Graph<int> myGraph;
	stringstream ss;
	ifstream f(name);
	int s,t,w;

	while(!f.eof()) {
		f>>s;
		f>>t;
		f>>w;
		myGraph.addVertex(s);
		myGraph.addVertex(t);
		myGraph.addEdge(s, t, w, 0);
	}


	return myGraph;
}

int main()
{
	time_t h1,h2;
	Graph<int> graph = readFile("nosc3.txt");
	cout << graph.getNumVertex() << endl;
	/*time(&h1);
	cout << graph.FordFulk(0,1) << endl;
	time(&h2);
	cout << h2-h1 << endl;
	graph = readFile("nosc3.txt");
	time(&h1);
	cout << graph.Dinic(0,1) << endl;
	time(&h2);
	cout << h2-h1 << endl;*/
	//plotGraph(graph);
	//cin >> h1;
	return 0;
}
