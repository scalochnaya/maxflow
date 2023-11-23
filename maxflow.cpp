#include <iostream>
#include <fstream>
#include <chrono>
#include <queue>
#include <vector>
#include <stack>

using namespace std;


class Graph
{
protected:
	vector<vector<int> > adjMatrix; // матрица смежности
	int vertex = 0;  // кол-во вершин


	bool bfs(vector<vector<int> >& rGraph, int source, int target, vector<int>& parent)
	{
		// BFS - поиск возможного минимального пути в сети
		queue<int> q;
		/*bool* visited = new bool[vertex];
		for (int i = 0; i < vertex; i++) visited[i] = false;*/
		vector<bool> visited(vertex, false);

		q.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v < vertex; v++)
			{
				if (!visited[v] && rGraph[u][v] > 0)
				{
					if (v == target)
					{
						parent[v] = u;
						//delete[] visited;
						return true;
					}
					q.push(v);
					parent[v] = u;
					visited[v] = true;
				}
			}
		}
		//delete[] visited;
		return false;
	}

	bool dfs(vector<vector<int> >& rGraph, int source, int target, vector<int>& parent)
	{
		// DFS; вместо очереди используется стек
		/*bool* visited = new bool[vertex];
		for (int i = 0; i < vertex; i++) visited[i] = false;*/
		vector<bool> visited(vertex, false);

		stack<int> S;

		S.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!S.empty())
		{
			int u = S.top();
			S.pop();

			for (int v = 0; v < vertex; v++)
			{
				if (!visited[v] && rGraph[u][v] > 0)
				{
					S.push(v);
					parent[v] = u;
					visited[v] = true;
				}
			}
		}
		bool res = visited[target];
		//delete[] visited;
		return res;
	}

	bool bfsForDinic(vector<vector<int> >& rGraph, int source, int target, vector<int>& level)
	{
		// BFS + создание слоистой вспомогательной сети
		queue<int> q;
		q.push(source);
		//for (int i = 0; i < vertex; i++) level[i] = -1; // это точно assign ????
		level.assign(vertex, -1);
		level[source] = 0;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v < vertex; v++)
			{
				if (level[v] < 0 && rGraph[u][v] > 0)
				{
					level[v] = level[u] + 1;
					q.push(v);
				}
			}
		}
		return level[target] >= 0;
	}

	int dfsForDinic(vector<vector<int> >& rGraph, int u, int target, int flow, vector<int>& level, vector<int>& start)
	{
		if (u == target) return flow;

		for (int& v = start[u]; v < vertex; v++)
		{
			if (rGraph[u][v] > 0 && level[v] == level[u] + 1)
			{
				if (int currFlow = dfsForDinic(rGraph, v, target, min(flow, rGraph[u][v]), level, start))
				{
					rGraph[u][v] -= currFlow;
					rGraph[v][u] += currFlow;
					return currFlow;
				}
			}
		}
		return 0;
	}


public:
	Graph()
	{
		// Базовый конструктор; выделение дин.памяти
		vertex = 5;
		/*adjMatrix = new int* [vertex];
		for (int i = 0; i < vertex; i++)
		{
			adjMatrix[i] = new int[vertex];
			for (int j = 0; j < vertex; j++)
				adjMatrix[i][j] = 0;
		}*/
		adjMatrix.resize(vertex, vector<int>(vertex, 0));
	}

	Graph(int V)
	{
		// Конструктор, принимающий количество вершин
		vertex = V;
		/*adjMatrix = new int* [vertex];
		for (int i = 0; i < vertex; i++)
		{
			adjMatrix[i] = new int[vertex];
			for (int j = 0; j < vertex; j++)
				adjMatrix[i][j] = 0;
		}*/
		adjMatrix.resize(vertex, vector<int>(vertex, 0));
	}

	Graph(const Graph& G)
	{
		// Конструктор копий
		vertex = G.vertex;
		/*adjMatrix = new int*[vertex];
		for (int i = 0; i < vertex; i++)
		{
			adjMatrix[i] = new int[vertex];
			for (int j = 0; j < vertex; j++)
				adjMatrix[i][j] = G.adjMatrix[i][j];
		}*/
		adjMatrix.resize(vertex, vector<int>(vertex, 0));
		for (int u = 0; u < vertex; u++)
			for (int v = 0; v < vertex; v++)
				adjMatrix[u][v] = G.adjMatrix[u][v];
	}

	~Graph()
	{
		// Деструктор; освобождение памяти
		/*for (int i = 0; i < vertex; i++)
			delete[] adjMatrix[i];
		delete[] adjMatrix;
		adjMatrix = NULL;*/
	}


	void addEdge(int a, int b, int weight)
	{
		adjMatrix[a][b] = weight;
	}

	int fordFulkerson(int source, int target)
	{
		// Алгоритм Форда-Фалкерсона
		/*int** rGraph = new int* [vertex];
		for (int u = 0; u < vertex; u++)
		{
			rGraph[u] = new int[vertex];
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];
		}*/
		vector<vector<int> > rGraph(vertex, vector<int>(vertex, 0));
		for (int u = 0; u < vertex; u++)
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];

		vector<int> parent(vertex, 0);

		//int* parent = new int[vertex];
		int maxFlow = 0;

		while (dfs(rGraph, source, target, parent)) {
			int pathFlow = INT_MAX;
			for (int v = target; v != source; v = parent[v]) {
				int u = parent[v];
				pathFlow = min(pathFlow, rGraph[u][v]);
			}

			for (int v = target; v != source; v = parent[v]) {
				int u = parent[v];
				rGraph[u][v] -= pathFlow;
				rGraph[v][u] += pathFlow;
			}

			maxFlow += pathFlow;
		}

		/*for (int i = 0; i < vertex; i++)
			delete[] rGraph[i];
		delete[] rGraph;

		delete[] parent;*/

		return maxFlow;
	}

	int edmondsKarp(int source, int target)
	{
		// Алгоритм Эдмондса-Карпа
		/*int** rGraph = new int* [vertex];
		for (int u = 0; u < vertex; u++)
		{
			rGraph[u] = new int[vertex];
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];
		}*/
		vector<vector<int> > rGraph(vertex, vector<int>(vertex, 0));
		for (int u = 0; u < vertex; u++)
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];

		//int* parent = new int[vertex];
		vector<int> parent(vertex, 0);
		int maxFlow = 0;

		while (bfs(rGraph, source, target, parent)) {
			int pathFlow = INT_MAX;
			for (int v = target; v != source; v = parent[v]) {
				int u = parent[v];
				pathFlow = min(pathFlow, rGraph[u][v]);
			}

			for (int v = target; v != source; v = parent[v]) {
				int u = parent[v];
				rGraph[u][v] -= pathFlow;
				rGraph[v][u] += pathFlow;
			}

			maxFlow += pathFlow;
		}

		/*for (int i = 0; i < vertex; i++)
			delete[] rGraph[i];
		delete[] rGraph;

		delete[] parent;*/

		return maxFlow;
	}

	int dinic(int source, int target)
	{
		// Алгоритм Диница
		int maxFlow = 0;
		vector<int> level(vertex, -1);
		//int* level = new int[vertex];
		//for (int i = 0; i < vertex; i++) level[i] = -1;

		/*int** rGraph = new int* [vertex];
		for (int u = 0; u < vertex; u++)
		{
			rGraph[u] = new int[vertex];
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];
		}*/

		vector<vector<int> > rGraph(vertex, vector<int>(vertex, 0));
		for (int u = 0; u < vertex; u++)
			for (int v = 0; v < vertex; v++)
				rGraph[u][v] = adjMatrix[u][v];

		while (bfsForDinic(rGraph, source, target, level))
		{
			vector<int> start(vertex, 0);
			//int* start = new int[vertex];
			for (int i = 0; i < vertex; i++) start[i] = 0;

			while (int flow = dfsForDinic(rGraph, source, target, INT_MAX, level, start))
				maxFlow += flow;

			//delete[] start;
		}

		/*for (int i = 0; i < vertex; i++)
			delete[] rGraph[i];
		delete[] rGraph;

		delete[] level;*/

		return maxFlow;
	}
	
	friend ostream& operator<<(ostream& stream, Graph& G);
	friend istream& operator>>(istream& stream, Graph& G);
};

ostream& operator<<(ostream& stream, Graph& G)
{
	stream << G.vertex << endl;
	for (int i = 0; i < G.vertex; i++)
	{
		for (int j = 0; j < G.vertex; j++)
			stream << G.adjMatrix[i][j] << " ";
		stream << endl;
	}
	
	return stream;
}

istream& operator>>(istream& stream, Graph& G)
{
	if (typeid(stream) == typeid(ifstream))
	{
		int vertex;
		stream >> vertex;

		if (vertex != G.vertex)
		{
			/*for (int i = 0; i < G.vertex; i++)
				delete[] G.adjMatrix[i];
			delete[] G.adjMatrix;

			G.vertex = vertex;
			G.adjMatrix = new int* [vertex];
			for (int j = 0; j < vertex; j++)
				G.adjMatrix[j] = new int[vertex];*/
			G.adjMatrix.resize(vertex, vector<int>(vertex, 0));
		}
	}
	for (int i = 0; i < G.vertex; i++)
		for (int j = 0; j < G.vertex; j++)
			stream >> G.adjMatrix[i][j];
	
	return stream;
}



int main()
{
	/*
	Graph G(6);
	G.addEdge(0, 1, 16);
	G.addEdge(0, 2, 13);
	G.addEdge(1, 2, 10);
	G.addEdge(1, 3, 12);
	G.addEdge(2, 1, 4);
	G.addEdge(2, 4, 14);
	G.addEdge(3, 2, 9);
	G.addEdge(3, 5, 20);
	G.addEdge(4, 3, 7);
	G.addEdge(4, 5, 4);
	*/

	/*
	Graph G(4);
	G.addEdge(0, 1, 1000);
	G.addEdge(0, 2, 1000);
	G.addEdge(1, 2, 1);
	G.addEdge(1, 3, 1000);
	G.addEdge(2, 3, 1000);
	*/

	/*
	Graph G(10);
	ifstream fin; fin.open("test0.txt");
	if (fin)
	{
		fin >> G;
		fin.close();
	}
	*/

	/*
	Graph G(7);
	ifstream fin; fin.open("test1.txt");
	if (fin)
	{
		fin >> G;
		fin.close();
	}
	*/

	/*
	Graph G(4);
	ifstream fin; fin.open("test2.txt");
	if (fin)
	{
		fin >> G;
		fin.close();
	}
	*/
	Graph G(15);

	cout << G << endl;


	long int timer_1, timer_2;
	long int ff_total = 0, ek_total = 0, d_total = 0;
	int i;
	
	// Важно указать исток и сток
	int source = 0;
	int target = 4;

	// Измерение алгоритма Форда-Фалкерсона
	int ff_flow;
	timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	ff_flow = G.fordFulkerson(source, target);
	timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	ff_total += (timer_2 - timer_1);
	for (i = 0; i < 9; i++)
	{
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		G.fordFulkerson(source, target);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ff_total += (timer_2 - timer_1);
	}
	cout << "[Ford-Fulkerson algorithm]\tFlow: " << ff_flow << endl;
	cout << "[Ford-Fulkerson algorithm]\tTotal time: " << ff_total << " ns" << endl;


	// Измерение алгоритма Эдмондса-Карпа
	int ek_flow;
	timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	ek_flow = G.edmondsKarp(source, target);
	timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	ek_total += (timer_2 - timer_1);
	for (i = 0; i < 9; i++)
	{
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		G.edmondsKarp(source, target);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ek_total += (timer_2 - timer_1);
	}
	cout << "[Edmonds-Karp algorithm]\tFlow: " << ek_flow << endl;
	cout << "[Edmonds-Karp algorithm]\tTotal time: " << ek_total << " ns" << endl;


	// Измерение алгоритма Диница
	int d_flow;
	timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	d_flow = G.dinic(source, target);
	timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
	d_total += (timer_2 - timer_1);
	for (i = 0; i < 9; i++)
	{
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		G.dinic(source, target);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		d_total += (timer_2 - timer_1);
	}
	cout << "[Dinic algorithm]\t\tFlow: " << d_flow << endl;
	cout << "[Dinic algorithm]\t\tTotal time: " << d_total << " ns" << endl;

	return 0;
}