#include "Graph.h"

bool Vertex::removeEdgeTo(Vertex *d) {
	d->indegree--;
	typename vector<Edge >::iterator it= adj.begin();
	typename vector<Edge >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			adj.erase(it);
			return true;
		}
		else it++;
	}
	return false;
}

bool Vertex::hideEdgeTo(Vertex *d) {
	typename vector<Edge >::iterator it= adj.begin();
	typename vector<Edge >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			it->hidden=true;
			return true;
		}
		else it++;
	}
	return false;
}

bool Vertex::addFlowTo(Vertex *d, const double flow) {
	typename vector<Edge >::iterator it= adj.begin();
	typename vector<Edge >::iterator ite= adj.end();
	while (it!=ite) {
		if (it->dest == d) {
			it->flow+=flow;
			return true;
		}
		else it++;
	}
	return false;
}

Vertex::Vertex(int in): info(in), visited(false), indegree(0), dist(0), in(0), out(0) {
	path = NULL;
}

void Vertex::addEdge(Vertex *dest, double w) {
	Edge edgeD(dest,w);
	edgeD.orig = this;
	dest->in++;
	this->out++;
	adj.push_back(edgeD);
}

void Vertex::addEdge(Vertex *dest, double w, double f)
{
	Edge edgeD(dest, w, f);
	edgeD.orig = this;
	dest->in++;
	this->out++;
	adj.push_back(edgeD);
}

void Vertex::addEdge(Vertex *dest, double w, double f, bool hide)
{
	Edge edgeD(dest, w, f);
	edgeD.orig = this;
	//edgeD.hidden=hide;
	adj.push_back(edgeD);
}

int Vertex::getInfo() const {
	return this->info;
}

int Vertex::getDist() const {
	return this->dist;
}

vector<Edge > Vertex::getAdj() const {
	return this->adj;
}

Vertex* Vertex::getPath() const {
	return this->path;
}

int Vertex::getIndegree() const {
	return this->indegree;
}

void Vertex::updateEdgeFlow(unsigned int index, float f)
{
	if (index >= adj.size())
		return;
	adj[index].flow = f;
}


Edge::Edge(Vertex *d, double w, double f): dest(d), weight(w), flow(f), hidden(false){}

Edge::Edge(Vertex *d, double w, double f, bool rev): dest(d), weight(w), flow(f),hidden(rev){}

double Edge::getFlow() const {
	return flow;
}

double Edge::getWeight() const {
	return weight;
}

Vertex * Edge::getOri() {
	return orig;
}

Vertex * Edge::getDest() {
	return dest;
}


int Graph::getNumVertex() const {
	return vertexSet.size();
}

vector<Vertex * > Graph::getVertexSet() const {
	return vertexSet;
}

bool Graph::addVertex(const int &in) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	for (; it!=ite; it++)
		if ((*it)->info == in) return false;
	Vertex *v1 = new Vertex(in);
	vertexSet.push_back(v1);
	return true;
}

bool Graph::removeVertex(const int &in) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	for (; it!=ite; it++) {
		if ((*it)->info == in) {
			Vertex * v= *it;
			vertexSet.erase(it);
			typename vector<Vertex*>::iterator it1= vertexSet.begin();
			typename vector<Vertex*>::iterator it1e= vertexSet.end();
			for (; it1!=it1e; it1++) {
				(*it1)->removeEdgeTo(v);
			}

			typename vector<Edge >::iterator itAdj= v->adj.begin();
			typename vector<Edge >::iterator itAdje= v->adj.end();
			for (; itAdj!=itAdje; itAdj++) {
				itAdj->dest->indegree--;
			}
			delete v;
			return true;
		}
	}
	return false;
}

bool Graph::addEdge(const int &sourc, const int &dest, double w, double f) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex *vS, *vD;
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

bool Graph::addEdge(const int &sourc, const int &dest, double w, double f, bool hidden) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex *vS, *vD;
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

bool Graph::removeEdge(const int &sourc, const int &dest) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex *vS, *vD;
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

bool Graph::hideEdge(const int &sourc, const int &dest) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex *vS, *vD;
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

bool Graph::addFlow(const int &sourc, const int &dest, const double flow) {
	typename vector<Vertex*>::iterator it= vertexSet.begin();
	typename vector<Vertex*>::iterator ite= vertexSet.end();
	int found=0;
	Vertex *vS, *vD;
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

Vertex* Graph::getVertex(const int &v) const {
	for(unsigned int i = 0; i < vertexSet.size(); i++)
		if (vertexSet[i]->info == v) return vertexSet[i];
	return NULL;
}


vector<Vertex*> Graph::getSources() const {
	vector< Vertex* > buffer;
	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		if( vertexSet[i]->indegree == 0 ) buffer.push_back( vertexSet[i] );
	}
	return buffer;
}

vector<int> Graph::getPath(const int &origin, const int &dest){

	list<int> buffer;
	Vertex* v = getVertex(dest);

	buffer.push_front(v->info);
	while ( v->path != NULL &&  v->path->info != origin) {
		v = v->path;
		buffer.push_front(v->info);
	}
	if( v->path != NULL )
		buffer.push_front(v->path->info);


	vector<int> res;
	while( !buffer.empty() ) {
		res.push_back( buffer.front() );
		buffer.pop_front();
	}
	return res;
}


void Graph::unweightedShortestPath(const int &s) {

	for(unsigned int i = 0; i < vertexSet.size(); i++) {
		vertexSet[i]->path = NULL;
		vertexSet[i]->dist = INT_INFINITY;
	}

	Vertex* v = getVertex(s);
	v->dist = 0;
	queue< Vertex* > q;
	q.push(v);

	while( !q.empty() ) {
		v = q.front(); q.pop();
		for(unsigned int i = 0; i < v->adj.size(); i++) {
			if(!v->adj[i].hidden) {
				Vertex* w = v->adj[i].dest;
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


int Graph::FordFulk(const int &s, const int &t) {
	int min=INT_INFINITY, cand;
	Graph gr;
	vector<int> res;

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
			if(addFlow(res[i],res[i+1],min))
				addFlow(res[i+1],res[i],-min);


		if(*res.rbegin()!=t || res.size()<2)
			break;

		min=INT_INFINITY;
	}

	int flow=0;

	vector<Edge > e;
	e=getEdges();
	for(unsigned int i=0; i<e.size(); i++)
		if(e[i].dest->getInfo()==t)
			flow+=e[i].flow;

	return flow;
}

int Graph::Dinic(const int &s, const int &t) {
	int num=getNumVertex()-1;
	int min=INT_INFINITY, cand;
	Graph gr;
	vector<int> res;
	vector<Edge> er;
	vector<Edge> e;
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


vector<Edge> Graph::getEdges() {
	vector<Vertex *> v = getVertexSet();
	vector<Edge > e;
	vector<Edge > tmp;
	for(unsigned int i=0; i<v.size(); i++) {
		tmp=v[i]->getAdj();
		e.insert(e.end(),tmp.begin(),tmp.end());
	}

	return e;
}


int Graph::getEdgeFlow(const int &s, const int &t) {
	vector<Vertex *> v = getVertexSet();
	vector<Edge > tmp;
	for(unsigned int i=0; i<v.size(); i++) {
		if(v[i]->getInfo()==s)
			tmp=v[i]->getAdj();
		for(unsigned int a=0; a<tmp.size(); a++)
			if(tmp[a].dest->getInfo()==t)
				return tmp[a].weight;
	}

	return 0;
}

Graph Graph::getGr() {
	vector<Edge> e=getEdges();
	Graph ret;

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

void Graph::printEdges(vector<Edge> e) {
	for(unsigned int i=0; i<e.size(); i++)
		cout << e[i].hidden << "  " << e[i].getFlow() << "/" << e[i].getWeight() << "  " << e[i].orig->getInfo() << " " << e[i].getDest()->getInfo() << "\n";
}

bool Graph::isNumber(string str)
{
	bool ret = false;
	for (unsigned int i = 0; i < str.size(); i++)
	{
		if (isdigit(str[i]))
			ret = true;
		else
			return false;
	}
	return ret;
}

int Graph::menuEdges() {

	string op;
	int num1, num2, capac;

	cout << "Indique a fonte da aresta: ";
	cin >> op;

	if(isNumber(op.c_str()))
		num1=atoi(op.c_str());
	else{
		cout << "Input invalido" << endl;
		return -1;
	}

	cout << "Indique o destino da aresta:";
	cin >> op;
	if(isNumber(op.c_str()))
		num2=atoi(op.c_str());
	else {
		cout << "Input invalido" << endl;
		return -1;
	}

	cout << "Indique a nova capacidade da aresta:";
	cin >> op;
	if(isNumber(op.c_str()))
		capac=atoi(op.c_str());
	else{
		cout << "Input invalido" << endl;
		return -1;
	}

	bool ret=this->removeEdge(num1, num2);;

	if(ret==false) {
		cout << endl <<  "Aresta invalida!" << endl;
		return -1;
	}

	this->addEdge(num1, num2, capac);

	return 0;
}


void Graph::plotGraph() {
	stringstream ss;
	GraphViewer *gv = new GraphViewer(600, 600, true);
	gv->createWindow(600, 600);

	vector<Edge > e=this->getEdges();
	vector<Vertex *> v=this->getVertexSet();

	for(unsigned int i=0; i<v.size(); i++) {
		gv->addNode(v[i]->getInfo());
		if((double)v[i]->out/(double)v[i]->in < 0.7 && v[i]->getInfo() != this->getNumVertex()-1)
			gv->setVertexColor(v[i]->getInfo(),"red");
		else gv->setVertexColor(v[i]->getInfo(),"green");
	}

	for(unsigned int i=0; i<e.size(); i++) {
		ss.str("");
		ss << e[i].getFlow() << "/" << e[i].getWeight();
		gv->addEdge(i,e[i].getOri()->getInfo(),e[i].getDest()->getInfo(),EdgeType::DIRECTED);
		gv->setEdgeLabel(i,ss.str());
		gv->setEdgeThickness(i,e[i].getFlow()*0.5);

		if(e[i].getFlow()==e[i].getWeight())
			gv->setEdgeColor(i,"red");
		else gv->setEdgeColor(i,"green");
	}

	gv->rearrange();

	cin.ignore(1000,'\n');
	cin.clear();
	getchar();
	gv->closeWindow();
}


void Graph::readFile(string name) {

	stringstream ss;
	ifstream f(name.c_str());
	if(!f.is_open())
		throw WrongFile(name);
	int s,t,w, edges;

	f >> s; f >> edges;
	for (int i = 0; i < edges; i++) {
		f>>s;
		f>>t;
		f>>w;
		this->addVertex(s);
		this->addVertex(t);
		this->addEdge(s, t, w, 0);
	}
}

void Graph::writeFile( string name) {
	ofstream f(name.c_str());
	vector<Edge > e=this->getEdges();

	f << this->getNumVertex() << " " << e.size() << endl;

	for(unsigned int i=0; i<e.size(); i++)
		f << e[i].getOri()->getInfo() << " " << e[i].getDest()->getInfo() << " " << e[i].getWeight() << endl;
}

void Graph::menuFord()
{
	clock_t begin, end;
	double time_spent;

	begin = clock();
	//cout << "source: " << 0 << " sink: " << this->getNumVertex()-1 << endl;
	int max = this->FordFulk(0,this->getNumVertex()-1) ;

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;

	cout << "\nFluxo maximo: "<< max << endl;
	cout << "Intervalo de tempo: " << time_spent << endl;

}

void Graph::menuDinic()
{
	clock_t begin, end;
	double time_spent;

	begin = clock();

	int max = this->Dinic(0,this->getNumVertex()-1);

	end = clock();

	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	cout << "\nFluxo maximo: "<< max << endl;
	cout << "Intervalo de tempo: " << time_spent << endl;

}

void Graph::resetEdgeFlow()
{
	for (unsigned int i = 0; i < vertexSet.size(); i++)
	{
		for (unsigned int a = 0; a < vertexSet[i]->adj.size(); a++)
			vertexSet[i]->adj[a].flow = 0;
	}
}

