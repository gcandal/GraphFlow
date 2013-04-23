/*
 * Graph.h
 */
#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <queue>
#include <list>
#include <climits>
#include <cmath>
#include <iostream>
#include <Windows.h>
#include <winbase.h>
using namespace std;

template <class T> class Edge;
template <class T> class Graph;

const int INT_INFINITY = INT_MAX;

/*
 * ================================================================================================
 * Class Vertex
 * ================================================================================================
 */
template <class T>
class Vertex {
	T info;
	vector<Edge<T>  > adj;
	bool visited;
	int indegree;
	double dist;
public:
	Vertex(T in);
	friend class Graph<T>;
	int in;
	int out;
	void addEdge(Vertex<T> *dest, double w);
	void addEdge(Vertex<T> *dest, double w, double f);
	void addEdge(Vertex<T> *dest, double w, double f, bool hide);
	bool removeEdgeTo(Vertex<T> *d);
	bool hideEdgeTo(Vertex<T> *d);
	bool addFlowTo(Vertex<T> *d, const double flow);

	T getInfo() const;

	int getDist() const;
	int getIndegree() const;
	vector<Edge<T> > getAdj() const;

	Vertex<T>* getPath() const;

	Vertex* path;

	void updateEdgeFlow(unsigned int index, float f);
};


template <class T>
struct vertex_greater_than {
	bool operator()(Vertex<T> * a, Vertex<T> * b) const {
		return a->getDist() > b->getDist();
	}
};


template <class T>
bool Vertex<T>::removeEdgeTo(Vertex<T> *d) {
	d->indegree--;
	typename vector<Edge<T> >::iterator it= adj.begin();
	typename vector<Edge<T> >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			adj.erase(it);
			return true;
		}
		else it++;
	}
	return false;
}

template <class T>
bool Vertex<T>::hideEdgeTo(Vertex<T> *d) {
	typename vector<Edge<T> >::iterator it= adj.begin();
	typename vector<Edge<T> >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			it->hidden=true;
			return true;
		}
		else it++;
	}
	return false;
}

template <class T>
bool Vertex<T>::addFlowTo(Vertex<T> *d, const double flow) {
	typename vector<Edge<T> >::iterator it= adj.begin();
	typename vector<Edge<T> >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			it->flow+=flow;
			return true;
		}
		else it++;
	}
	return false;
}

//atualizado pelo exercício 5
template <class T>
Vertex<T>::Vertex(T in): info(in), visited(false), indegree(0), dist(0), in(0), out(0) {
	path = NULL;
}


template <class T>
void Vertex<T>::addEdge(Vertex<T> *dest, double w) {
	Edge<T> edgeD(dest,w);
	edgeD.orig = this;
	dest->in++;
	this->out++;
	adj.push_back(edgeD);
}

template <class T>
void Vertex<T>::addEdge(Vertex<T> *dest, double w, double f)
{
	Edge<T> edgeD(dest, w, f);
	edgeD.orig = this;
	dest->in++;
	this->out++;
	adj.push_back(edgeD);
}

template <class T>
void Vertex<T>::addEdge(Vertex<T> *dest, double w, double f, bool hide)
{
	Edge<T> edgeD(dest, w, f);
	edgeD.orig = this;
	//edgeD.hidden=hide;
	adj.push_back(edgeD);
}

template <class T>
T Vertex<T>::getInfo() const {
	return this->info;
}

template <class T>
int Vertex<T>::getDist() const {
	return this->dist;
}

template <class T>
vector<Edge<T> > Vertex<T>::getAdj() const {
	return this->adj;
}

template <class T>
Vertex<T>* Vertex<T>::getPath() const {
	return this->path;
}

template <class T>
int Vertex<T>::getIndegree() const {
	return this->indegree;
}

template <class T>
void Vertex<T>::updateEdgeFlow(unsigned int index, float f)
{
	if (index >= adj.size())
		return;
	adj[index].flow = f;
}

/* ================================================================================================
 * Class Edge
 * ================================================================================================
 */
template <class T>
class Edge {
	Vertex<T> * dest;
	Vertex<T> * orig;
	double weight;
	double flow;
	bool hidden;
public:
	Edge(Vertex<T> *d, double w, double f=0);
	Edge(Vertex<T> *d, double w, double f, bool rev);
	double getFlow() const;
	double getWeight() const;
	Vertex<T> * getOri();
	Vertex<T> * getDest();

	friend class Graph<T>;
	friend class Vertex<T>;
};

template <class T>
Edge<T>::Edge(Vertex<T> *d, double w, double f): dest(d), weight(w), flow(f), hidden(false){}

template <class T>
Edge<T>::Edge(Vertex<T> *d, double w, double f, bool rev): dest(d), weight(w), flow(f),hidden(rev){}

template <class T>
double Edge<T>::getFlow() const {
	return flow;
}

template <class T>
double Edge<T>::getWeight() const {
	return weight;
}

template <class T>
Vertex<T> * Edge<T>::getOri() {
	return orig;
}

template <class T>
Vertex<T> * Edge<T>::getDest() {
	return dest;
}

/* ================================================================================================
 * Class Graph
 * ================================================================================================
 */
template <class T>
class Graph {
	vector<Vertex<T> *> vertexSet;
	void getPathTo(Vertex<T> *origin, list<T> &res);
	void printEdges(vector<Edge<T> > e);
	Graph<T> getGr();
	int getEdgeFlow(const T &s, const T &t);
public:
	vector<Edge<T> > getEdges();

	bool addVertex(const T &in);
	bool addEdge(const T &sourc, const T &dest, double w,double f=0);
	bool addEdge(const T &sourc, const T &dest, double w,double f,bool hidden);
	bool removeVertex(const T &in);
	bool removeEdge(const T &sourc, const T &dest);
	bool hideEdge(const T &sourc, const T &dest);
	bool addFlow(const T &sourc, const T &dest, const double flow);
	vector<Vertex<T> * > getVertexSet() const;
	int getNumVertex() const;

	Vertex<T>* getVertex(const T &v) const;
	vector<Vertex<T>*> getSources() const;
	vector<T> topologicalOrder();
	vector<T> getPath(const T &origin, const T &dest);
	void unweightedShortestPath(const T &v);

	int Dinic(const T &s, const T &t);
	int FordFulk(const T &s, const T &t);
};


template <class T>
int Graph<T>::getNumVertex() const {
	return vertexSet.size();
}
template <class T>
vector<Vertex<T> * > Graph<T>::getVertexSet() const {
	return vertexSet;
}

template <class T>
bool Graph<T>::addVertex(const T &in) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	for (; it!=ite; it++)
		if ((*it)->info == in) return false;
	Vertex<T> *v1 = new Vertex<T>(in);
	vertexSet.push_back(v1);
	return true;
}

template <class T>
bool Graph<T>::removeVertex(const T &in) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	for (; it!=ite; it++) {
		if ((*it)->info == in) {
			Vertex<T> * v= *it;
			vertexSet.erase(it);
			typename vector<Vertex<T>*>::iterator it1= vertexSet.begin();
			typename vector<Vertex<T>*>::iterator it1e= vertexSet.end();
			for (; it1!=it1e; it1++) {
				(*it1)->removeEdgeTo(v);
			}

			typename vector<Edge<T> >::iterator itAdj= v->adj.begin();
			typename vector<Edge<T> >::iterator itAdje= v->adj.end();
			for (; itAdj!=itAdje; itAdj++) {
				itAdj->dest->indegree--;
			}
			delete v;
			return true;
		}
	}
	return false;
}


template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w, double f) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex<T> *vS, *vD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == sourc )
		{ vS=*it; found++;}
		if ( (*it)->info == dest )
		{ vD=*it; found++;}
		it ++;
	}
	if (found!=2) return false;
	vD->indegree++;
	vS->addEdge(vD,w,f);

	return true;
}

template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w, double f, bool hidden) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex<T> *vS, *vD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == sourc )
		{ vS=*it; found++;}
		if ( (*it)->info == dest )
		{ vD=*it; found++;}
		it ++;
	}
	if (found!=2) return false;
	vD->indegree++;
	vS->addEdge(vD,w,f,hidden);

	return true;
}



template <class T>
bool Graph<T>::removeEdge(const T &sourc, const T &dest) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex<T> *vS, *vD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == sourc )
		{ vS=*it; found++;}
		if ( (*it)->info == dest )
		{ vD=*it; found++;}
		it ++;
	}
	if (found!=2)
		return false;

	vD->indegree--;

	return vS->removeEdgeTo(vD);
}

template <class T>
bool Graph<T>::hideEdge(const T &sourc, const T &dest) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex<T> *vS, *vD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == sourc )
		{ vS=*it; found++;}
		if ( (*it)->info == dest )
		{ vD=*it; found++;}
		it ++;
	}
	if (found!=2)
		return false;

	vD->indegree--;

	return vS->hideEdgeTo(vD);
}

template <class T>
bool Graph<T>::addFlow(const T &sourc, const T &dest, const double flow) {
	typename vector<Vertex<T>*>::iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex<T> *vS, *vD;
	while (found!=2 && it!=ite ) {
		if ( (*it)->info == sourc )
		{ vS=*it; found++;}
		if ( (*it)->info == dest )
		{ vD=*it; found++;}
		it ++;
	}
	if (found!=2)
		return false;

	vD->indegree--;

	return vS->addFlowTo(vD,flow);
}



template <class T>
Vertex<T>* Graph<T>::getVertex(const T &v) const {
	for(unsigned int i = 0; i < vertexSet.size(); i++)
		if (vertexSet[i]->info == v) return vertexSet[i];
	return NULL;
}

template<class T>
vector<Vertex<T>*> Graph<T>::getSources() const {
	vector< Vertex<T>* > buffer;
	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		if( vertexSet[i]->indegree == 0 ) buffer.push_back( vertexSet[i] );
	}
	return buffer;
}

template<class T>
vector<T> Graph<T>::getPath(const T &origin, const T &dest){

	list<T> buffer;
	Vertex<T>* v = getVertex(dest);

	buffer.push_front(v->info);
	while ( v->path != NULL &&  v->path->info != origin) {
		v = v->path;
		buffer.push_front(v->info);
	}
	if( v->path != NULL )
		buffer.push_front(v->path->info);


	vector<T> res;
	while( !buffer.empty() ) {
		res.push_back( buffer.front() );
		buffer.pop_front();
	}
	return res;
}



template<class T>
void Graph<T>::unweightedShortestPath(const T &s) {

	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		vertexSet[i]->path = NULL;
		vertexSet[i]->dist = INT_INFINITY;
	}

	Vertex<T>* v = getVertex(s);
	v->dist = 0;
	queue< Vertex<T>* > q;
	q.push(v);

	while( !q.empty() ) {
		v = q.front(); q.pop();
		for(unsigned int i = 0; i < v->adj.size(); i++) {
			if(!v->adj[i].hidden) {
				Vertex<T>* w = v->adj[i].dest;
				if( w->dist == INT_INFINITY ) {
					w->dist = v->dist + 1;
					w->path = v;
					q.push(w);
				}
			}
		}
	}
}

void printSquareArray(int ** arr, unsigned int size)
{
	for(unsigned int k = 0; k < size; k++)
	{
		if(k == 0)
		{
			cout <<  "   ";
			for(unsigned int i = 0; i < size; i++)
				cout <<  " " << i+1 << " ";
			cout << endl;
		}

		for(unsigned int i = 0; i < size; i++)
		{
			if(i == 0)
				cout <<  " " << k+1 << " ";

			if(arr[k][i] == INT_INFINITY)
				cout << " - ";
			else
				cout <<  " " << arr[k][i] << " ";
		}

		cout << endl;
	}
}

template <class T>
int Graph<T>::FordFulk(const T &s, const T &t) {
	int min=INT_INFINITY, cand;
	Graph<T> gr;
	vector<T> res;
	SYSTEMTIME t1, t2;
	double a[]={0,0,0,0};
	while(1) {
		gr=getGr();
		GetSystemTime(&t1);
		gr.unweightedShortestPath(s);
		GetSystemTime(&t2);
		a[0]+=t2.wMilliseconds-t1.wMilliseconds;
		res=gr.getPath(s,t);
		GetSystemTime(&t1);
		a[1]+=t1.wMilliseconds-t2.wMilliseconds;
		for(unsigned int i=0; i<res.size()-1; i++) {
			cand=gr.getEdgeFlow(res[i],res[i+1]);
			if(cand<min)
				min=cand;
		}
		GetSystemTime(&t2);
		a[2]=t2.wMilliseconds-t1.wMilliseconds;/*
		for(unsigned int i=0; i<res.size()-1; i++)
			if(!addFlow(res[i],res[i+1],min))
				addFlow(res[i+1],res[i],-min);*/
		for(unsigned int i=0; i<res.size()-1; i++)
			if(getEdgeFlow(res[i],res[i+1])==0)
				addFlow(res[i+1],res[i],-min);
			else addFlow(res[i+1],res[i],min);
		GetSystemTime(&t1);
		if(*res.rbegin()!=t || res.size()<2)
			break;
		a[3]=t1.wMilliseconds-t2.wMilliseconds;
		min=INT_INFINITY;
	}

	int flow=0;

	vector<Edge<T> > e;
	e=getEdges();
	for(unsigned int i=0; i<e.size(); i++)
		if(e[i].dest->getInfo()==t)
			flow+=e[i].flow;

	cout << "unweightedShortestPath" << a[0] << endl;
	cout << "getPath" << a[1] << endl;
	cout << "getEdgeFlow" << a[2] << endl;
	cout << "addFlow" << a[3] << endl;

	return flow;
}

template <class T>
int Graph<T>::Dinic(const T &s, const T &t) {
	int num=getNumVertex()-1;
	int min=INT_INFINITY, cand;
	Graph<T> gr;
	vector<T> res;
	vector<Edge<T> > er;
	vector<Edge<T> > e;
	int a=0;
	while(1) {
		gr=getGr();
		gr.unweightedShortestPath(s);

		if(gr.getVertex(t)->dist==INT_INFINITY)
			break;

		er=gr.getEdges();
		for(unsigned int i=0; i<er.size(); i++)
			if(er[i].orig->getDist()+1!=er[i].getDest()->getDist())
				gr.hideEdge(er[i].orig->getInfo(),er[i].getDest()->getInfo());

		res=gr.getPath(s,t);

		for(unsigned int i=0; i<res.size()-1; i++) {
			cand=gr.getEdgeFlow(res[i],res[i+1]);
			if(cand<min)
				min=cand;
		}

		for(unsigned int i=0; i<res.size()-1; i++)
			if(!addFlow(res[i],res[i+1],min))
				addFlow(res[i+1],res[i],-min);

		e=getEdges();

		if(getVertex(t)->dist==num)
			break;

		min=INT_INFINITY;
		a++;
	}

	int flow=0;

	e=getEdges();
	for(unsigned int i=0; i<e.size(); i++)
		if(e[i].dest->getInfo()==t)
			flow+=e[i].flow;

	return flow;
}

template<class T>
vector<Edge<T> > Graph<T>::getEdges() {
	vector<Vertex<T> *> v = getVertexSet();
	vector<Edge<T> > e;
	vector<Edge<T> > tmp;
	for(unsigned int i=0; i<v.size(); i++) {
		tmp=v[i]->getAdj();
		e.insert(e.end(),tmp.begin(),tmp.end());
	}

	return e;
}

template<class T>
int Graph<T>::getEdgeFlow(const T &s, const T &t) {
	vector<Vertex<T> *> v = getVertexSet();
	vector<Edge<T> > tmp;
	for(unsigned int i=0; i<v.size(); i++) {
		if(v[i]->getInfo()==s)
			tmp=v[i]->getAdj();
		for(unsigned int a=0; a<tmp.size(); a++)
			if(tmp[a].dest->getInfo()==t)
				return tmp[a].weight;
	}

	return 0;
}

template<class T>
Graph<T> Graph<T>::getGr() {
	vector<Edge<T> > e=getEdges();
	Graph<T> ret;

	for (unsigned int i = 0; i < this->vertexSet.size(); i++)
		ret.addVertex(this->vertexSet[i]->info);

	for (unsigned int a = 0; a < e.size(); a++) {
		if(e[a].flow<e[a].weight)
			ret.addEdge(e[a].orig->info, e[a].dest->info, e[a].weight-e[a].flow, 0);

		if(e[a].flow>0)
			ret.addEdge(e[a].dest->info, e[a].orig->info, e[a].flow, 0, true);

	}

	return ret;
}

template<class T>
void Graph<T>::printEdges(vector<Edge<T> > e) {
	for(unsigned int i=0; i<e.size(); i++)
		cout << e[i].hidden << "  " << e[i].getFlow() << "/" << e[i].getWeight() << "  " << e[i].orig->getInfo() << " " << e[i].getDest()->getInfo() << "\n";
}

#endif /* GRAPH_H_ */
