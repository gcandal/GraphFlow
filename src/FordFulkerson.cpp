#include "FordFulkerson.h"


int FordFulkerson::min (int x, int y) {
	return x < y ? x : y;
}

void FordFulkerson::enqueue(int x) {
	q.push(x);
	color[x]= GRAY1;
}

int FordFulkerson::dequeue() {
	int x = q.front();
	q.pop();
	color[x] = BLACK1;
	return x;

}

bool FordFulkerson::bfs(int start, int target) {
	for(int i = 0; i < nodes; i++)
		color[i] = WHITE1;

	enqueue(start);
	pred[start] = -1;
	while(!q.empty()) {
		int u = dequeue();
		for(int i = 0; i < nodes; i++) {
			if(color[i] == WHITE1 && capacity[u][i] - flow[u][i] > 0) {
				enqueue(i);
				pred[i] = u;
			}
		}
	}
	return color[target] == BLACK1;
}

void FordFulkerson::read_input_file() {

	ifstream input;
	input.open(file_name.c_str());
	if(!input.is_open())
		throw WrongFile(file_name);
	input >> nodes >> edges;

	for (int i = 0; i < nodes;i++) {
		for(int j = 0; j < nodes; j++)
			capacity[i][j] = 0;
	}

	for (int i = 0; i < edges; i++)
	{
		int a,b,c;
		input >> a >> b >> c;
		g->addNode(a);
		gr->addNode(a);
		g->addNode(b);
		gr->addNode(b);
		capacity[a][b] = c;
	}

	input.close();
}

FordFulkerson::FordFulkerson(string file) {
	ifstream input;
	input.open(file.c_str());
	if(!input.is_open())
		throw WrongFile(file);
	input.close();

	nodes = 0;
	edges = 0;
	n_edges_g = 0;
	n_edges_gr = 0;

	g = new GraphViewer(600, 600, true);
	g->createWindow(600, 600);
	g->defineVertexColor("green");
	g->defineEdgeColor("green");

	gr = new GraphViewer(600, 600, true);
	gr->createWindow(600, 600);
	gr->defineVertexColor("blue");
	gr->defineEdgeColor("green");

	file_name = file;
	input.open(file_name.c_str());
	input >> nodes >> edges;
	input.close();

	capacity.resize(nodes);
	for(int i = 0; i < nodes; i++)
		capacity[i].resize(nodes);

	flow.resize(nodes);
	for(int i = 0; i < nodes; i++)
		flow[i].resize(nodes);

	color.resize(nodes);
	pred.resize(nodes);

	read_input_file();
}

int FordFulkerson::max_flow (int source, int sink) {

	sink = nodes-1;
	int max = 0;

	for(int i = 0; i < nodes;i++)
		for(int j = 0; j <nodes; j++)
			flow[i][j]=0;

	while(bfs(source,sink))
	{
		int increment = oo;
		for(int i = nodes-1; pred[i] >= 0; i = pred[i]) {
			increment = min (increment, capacity[pred[i]][i] - flow[pred[i]][i]);
		}

		for(int i = nodes-1; pred[i] >= 0; i = pred[i]) {
			flow[pred[i]][i] += increment;
			flow[i][pred[i]] -= increment;
		}
		max += increment;
	}
	return max;
}

void FordFulkerson::print() {
	for (int i = 0; i < n_edges_g;i++)
		g->removeEdge(i);

	for (int i = 0; i < n_edges_gr;i++)
		gr->removeEdge(i);

	n_edges_g = 0;
	n_edges_gr = 0;

	for (int i = 0; i < nodes; i++)
	{
		for (int j = 0; j < nodes; j++){
			if(capacity[i][j] != 0){

				stringstream ss;
				ss << flow[i][j] << "/" << capacity[i][j];
				g->addEdge(n_edges_g,i,j,EdgeType::DIRECTED);
				g->setEdgeLabel(n_edges_g,ss.str());
				g->setEdgeThickness(n_edges_g,0.5*flow[i][j]);
				//cout << i << " " << j << " (" << flow[i][j] << "/" << capacity[i][j] << ")" << endl;
				if(capacity[i][j] == flow[i][j])
					g->setEdgeColor(n_edges_g,"red");
				n_edges_g++;
			}
			if (flow[i][j] > 0) {
				stringstream ss;
				if(flow[i][j] < capacity[i][j])
				{
					gr->addEdge(n_edges_gr,i,j,EdgeType::DIRECTED);
					gr->setEdgeFlow(n_edges_gr,capacity[i][j] - flow[i][j]);
					n_edges_gr++;
				}

				gr->addEdge(n_edges_gr,j,i,EdgeType::DIRECTED);
				gr->setEdgeFlow(n_edges_gr,flow[i][j]);
				n_edges_gr++;
			}
		}
	}
	g->rearrange();
	gr->rearrange();
	return;
}

vector<int> FordFulkerson::getBottlenecks() const {
	vector<int> bottlenecks;

	for(int i = 0; i < nodes-1; i++)
	{
		int out_capacity = 0;
		int in_capacity = 0;

		for (int j = 0; j < nodes ; j++) {
			out_capacity += capacity[i][j];
		}

		for (int j = 0; j < nodes ; j++) {
			in_capacity += capacity[j][i];
		}

		if ((double)out_capacity / (double)in_capacity < 0.7) {
			g->setVertexColor(i,"red");
			bottlenecks.push_back(i);
		}
	}
	return bottlenecks;
}

void FordFulkerson::setEdgeCapacity(int origin, int destination, int capac)
{
	if(origin == destination || origin >= nodes || destination >= nodes || origin < 0 || destination < 0 || capacity[origin][destination] < 1)
		throw(InvalidEdge(origin,destination));
	capacity[origin][destination] = capac;
}

void FordFulkerson::savetoFile(string name) {

	ofstream out;
	out.open(name.c_str());
	out << nodes << " " << edges << endl;
	for (int i = 0; i < nodes;i++)
		for (int j = 0; j < nodes; j++)
			if(capacity[i][j] > 0)
				out << i << " " << j << " " << capacity[i][j] << endl;
	out.close();
	return;
}

void FordFulkerson::menu() {

	clock_t begin, end;
	double time_spent;
	begin = clock();//
	int max = this->max_flow(0,1);
	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	this->print();
	cout << "\nFluxo maximo: "<< max << endl;
	cout << "Intervalo de tempo: " << time_spent  << endl;

	vector<int> bottlenecks = this->getBottlenecks();
	cout << "Potenciais bottlenecks: ";
	for (unsigned int i = 0; i < bottlenecks.size() - 1; i++)
		cout << bottlenecks[i] << ", ";
	cout << bottlenecks[bottlenecks.size()-1] << "." << endl;

	Sleep (5000);

	getchar();
}

int FordFulkerson::menuSetEdge() {

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


	try {
		this->setEdgeCapacity(num1, num2, capac);
	} catch(InvalidEdge &i) {
		cout << endl <<  "Aresta invalida!" << endl;
	}

	return 0;
}

bool FordFulkerson::isNumber(string str)
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

void FordFulkerson::closeGraphs() {
	g->closeWindow();
	gr->closeWindow();
}
