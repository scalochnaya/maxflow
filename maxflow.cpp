#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <stack>

using namespace std;

#define DEFAULT_VERTICES 10

class InvalidVertexException : public exception	// Исключение: Неверное кол-во вершин
{
protected:
	int invalidVertex;
public:
	InvalidVertexException(int i)
	{
		invalidVertex = i;
	}

	void print()
	{
		cout << "InvalidVertexException: Graph can't have " << invalidVertex << " vertices\n";
	}
};

class IndexOutOfBoundsException : public exception	// Исключение: Выход за границы
{
protected:
	int realVertices;
	int invalidVertex;
public:
	IndexOutOfBoundsException(int v, int i)
	{
		realVertices = v;
		invalidVertex = i;
	}
	void print()
	{
		cout << "IndexOutOfBoundsException: Graph has " << realVertices << ", not " << invalidVertex << endl;
	}
};

class InvalidAlgorithmException : public exception	// Исключение: Неверно указанный алгоритм
{
protected:
	int invalidAlgorithm;
public:
	InvalidAlgorithmException(int i)
	{
		invalidAlgorithm = i;
	}
	void print()
	{
		cout << "InvalidAlgorithmException: There is no algorithm with number " << invalidAlgorithm << endl;
	}
};


class Vertex					// Класс вершины
{
protected:
	vector<int> adjList;		// Список смежности для вершины
	int vertices = 0;			// Количество вершин

	void initializeVector(vector<int>& list, int v)
	{
		// Инициализация вершины
		vertices = v;
		list.resize(vertices, 0);
	}

public:
	Vertex()
	{
		// Базовый конструктор
		initializeVector(adjList, DEFAULT_VERTICES);
	}

	Vertex(int V)
	{
		// Конструктор, принимающий количество вершин
		if (V <= 0)
			throw InvalidVertexException(V);
		initializeVector(adjList, V);
	}

	Vertex(const Vertex& ver)
	{
		// Конструктор копий
		initializeVector(adjList, ver.vertices);
		for (int i = 0; i < vertices; i++)
			adjList[i] = ver.adjList[i];
	}

	void addEdge(int pos, int value)
	{
		// Ручное исправление значения в списке смежности для вершины
		if (pos >= vertices)
		{
			throw IndexOutOfBoundsException(vertices, pos);
		}
		adjList[pos] = value;
	}

	int getVertices() { return vertices; }

	int& operator[](int pos)
	{
		// Переопределение оператора []
		if (pos >= vertices)
			throw IndexOutOfBoundsException(vertices, pos);
		return adjList[pos];
	}

	void printVertex()
	{
		// Вывод списка смежности вершины в консоль
		for (int i = 0; i < vertices; i++)
			cout << adjList[i] << " ";
		cout << endl;
	}

	friend ostream& operator<<(ostream& stream, Vertex& V)
	{
		for (int i = 0; i < V.vertices; i++)
			stream << V.adjList[i] << " ";
		stream << endl;

		return stream;
	}
};


class Graph						// Класс графа, в котором реализованы вспомогательные функции
{
protected:
	vector<Vertex> adjMatrix;	// Матрица смежности (массив из вершин)
	int vertices = 0;			// кол-во вершин

	void initialize(vector<Vertex>& matrix, int Vertices)
	{
		vertices = Vertices;
		matrix.resize(vertices, Vertex(vertices));
	}

public:
	Graph()
	{
		// Базовый конструктор
		initialize(adjMatrix, DEFAULT_VERTICES);
	}

	Graph(int V)
	{
		// Конструктор, принимающий количество вершин
		if (V <= 0)
			throw InvalidVertexException(V);
		initialize(adjMatrix, V);
	}

	Vertex& operator[](int index)
	{
		// Переопределение оператора []
		if (index >= vertices)
			throw IndexOutOfBoundsException(vertices, index);
		return adjMatrix[index];
	}


	friend ostream& operator<<(ostream& stream, Graph& G)
	{
		// Вывод в файл/консоль
		stream << G.vertices << endl;
		for (int i = 0; i < G.vertices; i++)
		{
			stream << G[i];
		}

		return stream;
	}

	friend istream& operator>>(istream& stream, Graph& G)
	{
		// Вывод из файла/консоли
		if (typeid(stream) == typeid(ifstream))
		{
			int vertex;
			stream >> vertex;

			if (vertex != G.vertices)
			{
				G.initialize(G.adjMatrix, vertex);
			}
		}
		for (int i = 0; i < G.vertices; i++)
			for (int j = 0; j < G.vertices; j++)
				stream >> G.adjMatrix[i][j];

		return stream;
	}
};



class FlowGraph : public Graph		// Класс графа, предназначенного для поиска максимального потока
{
protected:
	bool bfs(vector<Vertex>& rGraph, int source, int target, vector<int>& parent)
	{
		// BFS - поиск возможного минимального пути в сети
		vector<bool> visited(vertices, false);
		queue<int> q;

		q.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v < vertices; v++)
			{
				if (!visited[v] && rGraph[u][v] > 0)
				{
					if (v == target)
					{
						parent[v] = u;
						return true;
					}
					q.push(v);
					parent[v] = u;
					visited[v] = true;
				}
			}
		}

		return false;
	}

	bool dfs(vector<Vertex>& rGraph, int source, int target, vector<int>& parent)
	{
		// DFS; вместо очереди используется стек
		vector<bool> visited(vertices, false);

		stack<int> S;

		S.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!S.empty())
		{
			int u = S.top();
			S.pop();

			for (int v = 0; v < vertices; v++)
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
		return res;
	}

	bool bfsForDinic(vector<Vertex>& rGraph, int source, int target, vector<int>& level)
	{
		// BFS + создание слоистой вспомогательной сети
		queue<int> q;
		q.push(source);
		level.assign(vertices, -1);
		level[source] = 0;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v < vertices; v++)
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

	int dfsForDinic(vector<Vertex>& rGraph, int u, int target, int flow, vector<int>& level, vector<int>& start)
	{
		if (u == target) return flow;

		for (int& v = start[u]; v < vertices; v++)
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

	int fordFulkerson(int source, int target)
	{
		// Алгоритм Форда-Фалкерсона
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		for (int u = 0; u < vertices; u++)
			for (int v = 0; v < vertices; v++)
				rGraph[u][v] = adjMatrix[u][v];

		vector<int> parent(vertices, 0);

		int maxFlow = 0;

		while (dfs(rGraph, source, target, parent))
		{
			int pathFlow = INT_MAX;
			for (int v = target; v != source; v = parent[v])
			{
				int u = parent[v];
				pathFlow = min(pathFlow, rGraph[u][v]);
			}

			for (int v = target; v != source; v = parent[v])
			{
				int u = parent[v];
				rGraph[u][v] -= pathFlow;
				rGraph[v][u] += pathFlow;
			}

			maxFlow += pathFlow;
		}

		return maxFlow;
	}

	int edmondsKarp(int source, int target)
	{
		// Алгоритм Эдмондса-Карпа
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		for (int u = 0; u < vertices; u++)
			for (int v = 0; v < vertices; v++)
				rGraph[u][v] = adjMatrix[u][v];

		vector<int> parent(vertices, 0);
		int maxFlow = 0;

		while (bfs(rGraph, source, target, parent))
		{
			int pathFlow = INT_MAX;
			for (int v = target; v != source; v = parent[v])
			{
				int u = parent[v];
				pathFlow = min(pathFlow, rGraph[u][v]);
			}

			for (int v = target; v != source; v = parent[v])
			{
				int u = parent[v];
				rGraph[u][v] -= pathFlow;
				rGraph[v][u] += pathFlow;
			}

			maxFlow += pathFlow;
		}

		return maxFlow;
	}

	int dinic(int source, int target)
	{
		// Алгоритм Диница
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		for (int u = 0; u < vertices; u++)
			for (int v = 0; v < vertices; v++)
				rGraph[u][v] = adjMatrix[u][v];

		int maxFlow = 0;
		vector<int> level(vertices, -1);

		while (bfsForDinic(rGraph, source, target, level))
		{
			vector<int> start(vertices, 0);

			while (int flow = dfsForDinic(rGraph, source, target, INT_MAX, level, start))
				maxFlow += flow;
		}

		return maxFlow;
	}
public:
	// Наследование конструкторов из класса родителя Graph
	FlowGraph() : Graph() {}
	FlowGraph(int V) : Graph(V) {}

	int findMaxFlow(int source, int target, int algo)
	{
		// Функция, которая передает управление тому или иному алгоритму
		// в зависимости от параметра algo

		if (algo == 1)
		{
			return fordFulkerson(source, target);
		}
		else if (algo == 2)
		{
			return edmondsKarp(source, target);
		}
		else if (algo == 3)
		{
			return dinic(source, target);
		}
		else
		{
			throw InvalidAlgorithmException(algo);
		}
	}
};



int main()
{
	try
	{
		/*
		FlowGraph G(10);
		ifstream fin; fin.open("test0.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}
		*/

		/*
		FlowGraph G(4);
		ifstream fin; fin.open("test1.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}
		*/

		
		/*FlowGraph G(6);
		ifstream fin; fin.open("test2.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}*/
		

		FlowGraph G(11);
		ifstream fin; fin.open("testx.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}

		cout << G << endl;



		long int timer_1, timer_2;
		long int ff_total = 0, ek_total = 0, d_total = 0;
		int i;

		//Важно указать исток и сток
		int source = 0;
		int target = 10;

		//Измерение алгоритма Форда-Фалкерсона
		int ff_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ff_flow = G.findMaxFlow(source, target, 1);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ff_total += (timer_2 - timer_1);
		for (i = 0; i < 9; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(source, target, 1);
			timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			ff_total += (timer_2 - timer_1);
		}
		cout << "[Ford-Fulkerson algorithm]\tFlow: " << ff_flow << endl;
		cout << "[Ford-Fulkerson algorithm]\tTotal time: " << ff_total << " ns" << endl;


		//Измерение алгоритма Эдмондса-Карпа
		int ek_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ek_flow = G.findMaxFlow(source, target, 2);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ek_total += (timer_2 - timer_1);
		for (i = 0; i < 9; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(source, target, 2);
			timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			ek_total += (timer_2 - timer_1);
		}
		cout << "[Edmonds-Karp algorithm]\tFlow: " << ek_flow << endl;
		cout << "[Edmonds-Karp algorithm]\tTotal time: " << ek_total << " ns" << endl;


		//Измерение алгоритма Диница
		int d_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		d_flow = G.findMaxFlow(source, target, 3);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		d_total += (timer_2 - timer_1);
		for (i = 0; i < 9; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(source, target, 3);
			timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			d_total += (timer_2 - timer_1);
		}
		cout << "[Dinic algorithm]\t\tFlow: " << d_flow << endl;
		cout << "[Dinic algorithm]\t\tTotal time: " << d_total << " ns" << endl;

	}
	catch (InvalidVertexException e) { e.print(); }
	catch (IndexOutOfBoundsException e) { e.print(); }
	catch (InvalidAlgorithmException e) { e.print(); }
	catch (...) { cout << "Something went wrong..." << endl; }

	return 0;
}
