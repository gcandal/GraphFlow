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
using namespace std;

template <class T> class Edge;
template <class T> class Graph;

const int NOT_VISITED = 0;
const int BEING_VISITED = 1;
const int DONE_VISITED = 2;
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
	bool processing;
	int indegree;
	double dist;
	int set;

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
	void setInfo(T info);

	int getDist() const;
	int getIndegree() const;
	vector<Edge<T> > getAdj() const;

	Vertex<T>* getPath() const;

	bool operator<(const Vertex<T> vertex);

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
	d->indegree--; //adicionado do exercicio 5
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
Vertex<T>::Vertex(T in): info(in), visited(false), processing(false), indegree(0), dist(0), in(0), out(0) {
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
void Vertex<T>::setInfo(T info) {
	this->info = info;
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
public:
	bool hidden;
	Edge(Vertex<T> *d, double w, double f=0);
	Edge(Vertex<T> *d, double w, double f, bool rev);
	double getFlow() const;
	double getWeight() const;
	Vertex<T> *getDest() const;
	Vertex<T> *getOri() const;
	bool operator<(const Edge<T> &other) const;

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
Vertex<T>* Edge<T>::getDest() const {
	return dest;
}

template <class T>
Vertex<T>* Edge<T>::getOri() const {
	return orig;
}

template <class T>
bool Edge<T>::operator<(const Edge<T> &other) const {
	return this->weight < other.weight;
}

template <class T>
struct edge_greater_than {
	bool operator()(Edge<T> a, Edge<T>  b) const {
		return a.getWeight() > b.getWeight();
	}
};



/* ================================================================================================
 * Class Graph
 * ================================================================================================
 */
template <class T>
class Graph {
	vector<Vertex<T> *> vertexSet;
	void dfs(Vertex<T> *v, vector<T> &res) const;
	bool dfs(Vertex<T> *v, vector<T> &res, const T &s, int &min) const;
	//exercicio 5
	int numCycles;
	void dfsVisit(Vertex<T> *v);
	void dfsVisit();
	void getPathTo(Vertex<T> *origin, list<T> &res);

	//exercicio 6
	int ** W;   //weight
	int ** P;   //path


	void printEdges(vector<Edge<T> > e);
	Graph<T> getGr();
	vector<T> dfs(const T &s, const T &t, int &min) const;
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
	vector<T> dfs() const;
	vector<T> bfs(Vertex<T> *v) const;
	int maxNewChildren(Vertex<T> *v, T &inf) const;
	vector<Vertex<T> * > getVertexSet() const;
	int getNumVertex() const;

	//exercicio 5
	Vertex<T>* getVertex(const T &v) const;
	void resetIndegrees();
	vector<Vertex<T>*> getSources() const;
	int getNumCycles();
	vector<T> topologicalOrder();
	vector<T> getPath(const T &origin, const T &dest);
	void unweightedShortestPath(const T &v);
	bool isDAG();

	//exercicio 6
	void bellmanFordShortestPath(const T &s);
	void dijkstraShortestPath(const T &s);
	void floydWarshallShortestPath();
	int edgeCost(int vOrigIndex, int vDestIndex);
	vector<T> getfloydWarshallPath(const T &origin, const T &dest);
	void getfloydWarshallPathAux(int index1, int index2, vector<T> & res);

	//exercicio 8
	Graph<T> clone();
	void resetEdgeFlow();
	int Dinic(const T &s, const T &t);
	int FordFulk(const T &s, const T &t);
	//vector<Vertex<T>*> calculateFordFulkerson(T source);
	//float calculateFordFulkerson(Vertex<T>* current, Vertex<T>* parent, float min, Graph<T>* gf, Graph<T>* gr, vector<Vertex<T>*> visited);
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
int Graph<T>::getNumCycles() {
	numCycles = 0;
	dfsVisit();
	return this->numCycles;
}

template <class T>
bool Graph<T>::isDAG() {
	return (getNumCycles() == 0);
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
vector<T> Graph<T>::dfs() const {
	typename vector<Vertex<T>*>::const_iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::const_iterator ite= vertexSet.end();
	for (; it !=ite; it++)
		(*it)->visited=false;
	vector<T> res;
	it=vertexSet.begin();
	for (; it !=ite; it++)
		if ( (*it)->visited==false )
			dfs(*it,res);
	return res;
}

template <class T>
void Graph<T>::dfs(Vertex<T> *v,vector<T> &res) const {
	v->visited = true;
	res.push_back(v->info);
	typename vector<Edge<T> >::iterator it= (v->adj).begin();
	typename vector<Edge<T> >::iterator ite= (v->adj).end();
	for (; it !=ite; it++)
		if ( it->dest->visited == false && !it->hidden ){
			//cout << "ok ";
			dfs(it->dest, res);
		}
}

template <class T>
vector<T> Graph<T>::dfs(const T &s, const T &t, int &min) const {
	vector<Vertex<T>*> v;
	v.push_back(getVertex(s));
	typename vector<Vertex<T>*>::const_iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::const_iterator ite= vertexSet.end();
	for (; it !=ite; it++)
		(*it)->visited=false;
	vector<T> res;
	it=v.begin();
	ite=v.end();
	for (; it !=ite; it++) {
		if( (*it)->info==t )
			return res;
		if ( (*it)->visited==false ) {
			if(dfs(*it,res,t,min))
				return res;
		}
	}

	return res;
}

template <class T>
bool Graph<T>::dfs(Vertex<T> *v,vector<T> &res, const T &t, int &min) const {
	v->visited = true;
	cout << "add " << v->info << " " << (v->adj).size() << "\n";
	res.push_back(v->info);
	if(v->info==t)
		return true;
	typename vector<Edge<T> >::iterator it= (v->adj).begin();
	typename vector<Edge<T> >::iterator ite= (v->adj).end();
	for (; it !=ite; it++)
		if ( it->dest->visited == false /*&& !it->hidden*/) {
			if(it->weight<min)
				min=it->weight;
			return dfs(it->dest, res, t,min);
		}

	return false;
}

template <class T>
vector<T> Graph<T>::bfs(Vertex<T> *v) const {
	vector<T> res;
	queue<Vertex<T> *> q;
	q.push(v);
	v->visited = true;
	while (!q.empty()) {
		Vertex<T> *v1 = q.front();
		q.pop();
		res.push_back(v1->info);
		typename vector<Edge<T> >::iterator it=v1->adj.begin();
		typename vector<Edge<T> >::iterator ite=v1->adj.end();
		for (; it!=ite; it++) {
			Vertex<T> *d = it->dest;
			if (d->visited==false) {
				d->visited=true;
				q.push(d);
			}
		}
	}
	return res;
}

template <class T>
int Graph<T>::maxNewChildren(Vertex<T> *v, T &inf) const {
	vector<T> res;
	queue<Vertex<T> *> q;
	queue<int> level;
	int maxChildren=0;
	inf =v->info;
	q.push(v);
	level.push(0);
	v->visited = true;
	while (!q.empty()) {
		Vertex<T> *v1 = q.front();
		q.pop();
		res.push_back(v1->info);
		int l=level.front();
		level.pop(); l++;
		int nChildren=0;
		typename vector<Edge<T> >::iterator it=v1->adj.begin();
		typename vector<Edge<T> >::iterator ite=v1->adj.end();
		for (; it!=ite; it++) {
			Vertex<T> *d = it->dest;
			if (d->visited==false) {
				d->visited=true;
				q.push(d);
				level.push(l);
				nChildren++;
			}
		}
		if (nChildren>maxChildren) {
			maxChildren=nChildren;
			inf = v1->info;
		}
	}
	return maxChildren;
}


template <class T>
Vertex<T>* Graph<T>::getVertex(const T &v) const {
	for(unsigned int i = 0; i < vertexSet.size(); i++)
		if (vertexSet[i]->info == v) return vertexSet[i];
	return NULL;
}

template<class T>
void Graph<T>::resetIndegrees() {
	//colocar todos os indegree em 0;
	for(unsigned int i = 0; i < vertexSet.size(); i++) vertexSet[i]->indegree = 0;

	//actualizar os indegree
	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		//percorrer o vector de Edges, e actualizar indegree
		for(unsigned int j = 0; j < vertexSet[i]->adj.size(); j++) {
			vertexSet[i]->adj[j].dest->indegree++;
		}
	}
}


template<class T>
vector<Vertex<T>*> Graph<T>::getSources() const {
	vector< Vertex<T>* > buffer;
	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		if( vertexSet[i]->indegree == 0 ) buffer.push_back( vertexSet[i] );
	}
	return buffer;
}


template <class T>
void Graph<T>::dfsVisit() {
	typename vector<Vertex<T>*>::const_iterator it= vertexSet.begin();
	typename vector<Vertex<T>*>::const_iterator ite= vertexSet.end();
	for (; it !=ite; it++)
		(*it)->visited=false;
	it=vertexSet.begin();
	for (; it !=ite; it++)
		if ( (*it)->visited==false )
			dfsVisit(*it);
}

template <class T>
void Graph<T>::dfsVisit(Vertex<T> *v) {
	v->processing = true;
	v->visited = true;
	typename vector<Edge<T> >::iterator it= (v->adj).begin();
	typename vector<Edge<T> >::iterator ite= (v->adj).end();
	for (; it !=ite; it++) {
		if ( it->dest->processing == true) numCycles++;
		if ( it->dest->visited == false ){
			dfsVisit(it->dest);
		}
	}
	v->processing = false;
}

template<class T>
vector<T> Graph<T>::topologicalOrder() {
	//vetor com o resultado da ordenacao
	vector<T> res;

	//verificar se é um DAG
	if( getNumCycles() > 0 ) {
		cout << "Ordenacao Impossivel!" << endl;
		return res;
	}

	//garantir que os "indegree" estao inicializados corretamente
	this->resetIndegrees();

	queue<Vertex<T>*> q;

	vector<Vertex<T>*> sources = getSources();
	while( !sources.empty() ) {
		q.push( sources.back() );
		sources.pop_back();
	}


	//processar fontes
	while( !q.empty() ) {
		Vertex<T>* v = q.front();
		q.pop();

		res.push_back(v->info);

		for(unsigned int i = 0; i < v->adj.size(); i++) {
			v->adj[i].dest->indegree--;
			if( v->adj[i].dest->indegree == 0) q.push( v->adj[i].dest );
		}
	}

	//testar se o procedimento foi bem sucedido
	if ( res.size() != vertexSet.size() ) {
		while( !res.empty() ) res.pop_back();
	}

	//garantir que os "indegree" ficam atualizados ao final
	this->resetIndegrees();

	return res;
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
vector<T> Graph<T>::getfloydWarshallPath(const T &origin, const T &dest){

	int originIndex = -1, destinationIndex = -1;

	for(unsigned int i = 0; i < vertexSet.size(); i++)
	{
		if(vertexSet[i]->info == origin)
			originIndex = i;
		if(vertexSet[i]->info == dest)
			destinationIndex = i;

		if(originIndex != -1 && destinationIndex != -1)
			break;
	}


	vector<T> res;

	//se nao foi encontrada solucao possivel, retorna lista vazia
	if(W[originIndex][destinationIndex] == INT_INFINITY)
		return res;

	res.push_back(vertexSet[originIndex]->info);

	//se houver pontos intermedios...
	if(P[originIndex][destinationIndex] != -1)
	{
		int intermedIndex = P[originIndex][destinationIndex];

		getfloydWarshallPathAux(originIndex, intermedIndex, res);

		res.push_back(vertexSet[intermedIndex]->info);

		getfloydWarshallPathAux(intermedIndex,destinationIndex, res);
	}

	res.push_back(vertexSet[destinationIndex]->info);


	return res;
}



template<class T>
void Graph<T>::getfloydWarshallPathAux(int index1, int index2, vector<T> & res)
{
	if(P[index1][index2] != -1)
	{
		getfloydWarshallPathAux(index1, P[index1][index2], res);

		res.push_back(vertexSet[P[index1][index2]]->info);

		getfloydWarshallPathAux(P[index1][index2],index2, res);
	}
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

template<class T>
void Graph<T>::bellmanFordShortestPath(const T &s) {

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
			Vertex<T>* w = v->adj[i].dest;
			if(v->dist + v->adj[i].weight < w->dist) {
				w->dist = v->dist + v->adj[i].weight;
				w->path = v;
				q.push(w);
			}
		}
	}
}





template<class T>
void Graph<T>::dijkstraShortestPath(const T &s) {

	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		vertexSet[i]->path = NULL;
		vertexSet[i]->dist = INT_INFINITY;
		vertexSet[i]->processing = false;
	}

	Vertex<T>* v = getVertex(s);
	v->dist = 0;

	vector< Vertex<T>* > pq;
	pq.push_back(v);

	make_heap(pq.begin(), pq.end());


	while( !pq.empty() ) {

		v = pq.front();
		pop_heap(pq.begin(), pq.end());
		pq.pop_back();

		for(unsigned int i = 0; i < v->adj.size(); i++) {
			Vertex<T>* w = v->adj[i].dest;

			if(v->dist + v->adj[i].weight < w->dist ) {

				w->dist = v->dist + v->adj[i].weight;
				w->path = v;

				//se já estiver na lista, apenas a actualiza
				if(!w->processing)
				{
					w->processing = true;
					pq.push_back(w);
				}

				make_heap (pq.begin(),pq.end(),vertex_greater_than<T>());
			}
		}
	}
}

template<class T>
int Graph<T>::edgeCost(int vOrigIndex, int vDestIndex)
{
	if(vertexSet[vOrigIndex] == vertexSet[vDestIndex])
		return 0;

	for(unsigned int i = 0; i < vertexSet[vOrigIndex]->adj.size(); i++)
	{
		if(vertexSet[vOrigIndex]->adj[i].dest == vertexSet[vDestIndex])
			return vertexSet[vOrigIndex]->adj[i].weight;
	}

	return INT_INFINITY;
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


template<class T>
void Graph<T>::floydWarshallShortestPath() {

	W = new int * [vertexSet.size()];
	P = new int * [vertexSet.size()];
	for(unsigned int i = 0; i < vertexSet.size(); i++)
	{
		W[i] = new int[vertexSet.size()];
		P[i] = new int[vertexSet.size()];
		for(unsigned int j = 0; j < vertexSet.size(); j++)
		{
			W[i][j] = edgeCost(i,j);
			P[i][j] = -1;
		}
	}

	for(unsigned int k = 0; k < vertexSet.size(); k++)
		for(unsigned int i = 0; i < vertexSet.size(); i++)
			for(unsigned int j = 0; j < vertexSet.size(); j++)
			{
				//se somarmos qualquer coisa ao valor INT_INFINITY, ocorre overflow, o que resulta num valor negativo, logo nem convém considerar essa soma
				if(W[i][k] == INT_INFINITY || W[k][j] == INT_INFINITY)
					continue;

				int val = min ( W[i][j], W[i][k]+W[k][j] );
				if(val != W[i][j])
				{
					W[i][j] = val;
					P[i][j] = k;
				}
			}

}



template <class T>
void Graph<T>::resetEdgeFlow()
{
	for (unsigned int i = 0; i < vertexSet.size(); i++)
	{
		for (unsigned int a = 0; a < vertexSet[i]->adj.size(); a++)
			vertexSet[i]->adj[a].flow = 0;
	}
}

template <class T>
Graph<T> Graph<T>::clone()
{
	Graph<T> ret;
	for (unsigned int i = 0; i < this->vertexSet.size(); i++)
		ret.addVertex(this->vertexSet[i]->info);

	for (unsigned int i = 0; i < this->vertexSet.size(); i++)
	{
		vector<Edge<T> > edges = this->vertexSet[i]->adj;
		for (unsigned int a = 0; a < edges.size(); a++)
			ret.addEdge(this->vertexSet[i]->info, edges[a].dest->info, edges[a].weight, edges[a].flow, edges[a].hidden);
	}

	return ret;
}

template <class T>
int Graph<T>::FordFulk(const T &s, const T &t) {
	int min=INT_INFINITY, cand;
	Graph<T> gr;
	vector<T> res;

	while(1) {
		gr=getGr();
		gr.unweightedShortestPath(s);
		res=gr.getPath(s,t);

		for(unsigned int i=0; i<res.size()-1; i++) {
			cand=gr.getEdgeFlow(res[i],res[i+1]);
			if(cand<min)
				min=cand;
		}
		for(unsigned int i=0; i<res.size()-1; i++)
			if(!addFlow(res[i],res[i+1],min))
				addFlow(res[i+1],res[i],-min);

		if(*res.rbegin()!=t || res.size()<2)
			break;

		min=INT_INFINITY;
	}

	int flow=0;

	vector<Edge<T> > e;
	e=getEdges();
	for(unsigned int i=0; i<e.size(); i++)
		if(e[i].dest->getInfo()==t)
			flow+=e[i].flow;

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

	return -10;
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
