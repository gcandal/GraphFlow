#include "FordFulkerson.h"
#include "Graph.h"
#include "Exceptions.h"

using namespace std;

int menu() {

	stringstream ss;
	string input;
	int choice;

	cout << endl <<
			"1. Alterar a capacidade de uma aresta" << endl <<
			"2. Visualizar o grafo" << endl <<
			"3. Sair" << endl <<
			"Escolha uma opcao: ";

	cin >> input;
	ss << input;
	ss >> choice;

	cout << endl;

	switch(choice) {
	case 1:
		return 1;
	case 2:
		return 2;
	case 3:
		return 3;
	default: {
		cout << "Opcao invalida" << endl;
		return -1;
	}

	}
}


string escolherFicheiro() {

	string file;

	cout << endl << "Indique o nome do ficheiro:";
	cin >> file;

	return file;
}

int le(int l, int h)
{
	int opcao;
	cin >> opcao;
	while (cin.fail() || opcao < l || opcao > h)
	{
		cin.clear();
		cin.ignore(1000, '\n');
		cout << "Tenta outra vez: ";
		cin >> opcao;
	}
	return opcao;
}

void generateGraph() {
	int numNodes;
	string fileName;
	int flow[3];
	ofstream o;
	int source = 0;
	int sink;
	cout << "Numero de vertices: ";
	numNodes = le(2,10000);

	cout << "Capacidade passadeiras rolantes: ";
	flow[0] = le(1,1000);

	cout << "Capacidade escadas rolantes: ";
	flow[1] = le(1,1000);

	cout << "Capacidade carruagens automatizadas: ";
	flow[2] = le(1,1000);

	cout << "Nome do Ficheiro: " << endl;
	cin.clear();
	cin.ignore(1000,'\n');
	getline(cin,fileName);

	sink = numNodes-1;

	int adj[numNodes][numNodes];

	//cout << " " << flow[0] << " " << flow[1] << " " << flow[2] << endl;
	for(int u = 0; u < numNodes; u++) {
		for(int v = 0; v < u; v++) {
			if (rand() % 2 && ((u != sink && v != source) || u == source)) {
				adj[u][v] = flow[rand()%3];
				adj[v][u] = 0;
			}
			else {
				adj[u][v] = 0;
				adj[v][u] = flow[rand()%3];
			}
		}
	}

	o.open(fileName.c_str());
	int numEdges = 0;
	for(int u = 0; u < numNodes; u++)
		for(int v = 0; v < numNodes; v++)
			if(adj[u][v] > 0 && u != v)
				numEdges++;

	o << numNodes << " " << numEdges << endl;

	for(int u = 0; u < numNodes; u++)
		for(int v = 0; v < numNodes; v++)
			if(adj[u][v] > 0 && u != v)
				o << u << " " << v << " " << adj[u][v] << endl;

	o.close();
	cout << "Done!" << endl;
}

int escolherAlgoritmo() {

	stringstream ss;
	string input, file;
	int choice, menuOp;

	cout << "Que algoritmo pretende utilizar?" << endl <<
			"1. Ford Fulkerson (otimizado em tempo)" << endl <<
			"2. Ford Fulkerson " << endl <<
			"3. Dinic" << endl <<
			"4. Gerar Grafo" << endl <<
			"5. Sair" << endl <<
			"Escolha uma opcao: ";
	cin >> input;
	ss << input;
	ss >> choice;

	bool continueMenu=true;

	switch(choice) {

	case 1: {
		file=escolherFicheiro();

		try {
			FordFulkerson f1(file);
			f1.menu();

			do {
				menuOp = menu();

				switch(menuOp) {
				case 1:
					f1.menuSetEdge();
					break;
				case 2:
					f1.menu();
					break;
				case 3:
					f1.savetoFile(file);
					f1.closeGraphs();
					continueMenu=false;
					return 1;
				}

			} while(menuOp==-1 || continueMenu);


		} catch(WrongFile & w) {
			cout << endl << "O ficheiro " << w.getName() << " não existe!" << endl << endl;
			return -1;
		}
		break;
	}
	case 2: {
		file=escolherFicheiro();
		try {
			Graph graph;
			graph.readFile(file);
			graph.menuFord();

			do {
				menuOp = menu();

				switch(menuOp) {
				case 1:
					if (graph.menuEdges() == 0){
						graph.resetEdgeFlow();
						graph.menuFord();
					}
					break;
				case 2:
					graph.plotGraph();
					break;
				case 3:
					graph.writeFile(file);
					continueMenu=false;
					return 1;
				}

			} while(menuOp==-1 || continueMenu);
		} catch(WrongFile &w){
			cout << endl << "O ficheiro " << w.getName() << " não existe!" << endl << endl;
			return -1;
		}
		break;
	}
	case 3: {
		file=escolherFicheiro();
		try {
			Graph graph;
			graph.readFile(file);
			graph.menuDinic();

			do {
				menuOp = menu();

				switch(menuOp) {
				case 1:
					if (graph.menuEdges() == 0){
						graph.resetEdgeFlow();
						graph.menuFord();
					}
					break;
				case 2:
					graph.plotGraph();
					break;
				case 3:
					graph.writeFile(file);
					continueMenu=false;
					return 1;
				}

			} while(menuOp==-1 || continueMenu);
		} catch(WrongFile &w){
			cout << endl << "O ficheiro " << w.getName() << " não existe!" << endl << endl;
			return -1;
		}
		break;
	}
	case 4:
	{
		generateGraph();
		return -1;
		break;
	}
	case 5:
		return 0;
	default:
		cout << "Opcao invalida. \n" << endl;
		return -1;
	}

	return 0;
}


int main() {
	srand(time(NULL));
	int num;

	do {
		num=escolherAlgoritmo();
	}while(num!=0);

	return 0;
}
