#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <stack>

using namespace std;

#define DEFAULT_VERTICES 10

#define SOURCE 0
#define TARGET 10


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

	inline void initializeVector(vector<int>& list, int v)
	{
		vertices = v;
		list.resize(vertices);
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

		int i = 0;
		for (auto it : adjList)
		{
			it = ver.adjList[i];
			it++;  i++;
		}
	}

	void addEdge(int pos, int value)
	{
		// Ручное исправление значения в списке смежности для вершины
		if (pos - vertices >= 0)
			throw IndexOutOfBoundsException(vertices, pos);
		adjList[pos] = value;
	}

	int getVertices() { return vertices; }

	int& operator[](int pos)
	{
		// Переопределение оператора []
		if (pos - vertices >= 0)
			throw IndexOutOfBoundsException(vertices, pos);
		return adjList[pos];
	}

	void printVertex()
	{
		// Вывод списка смежности вершины в консоль
		for (auto it : adjList)
			cout << it++ << " ";
		cout << endl;
	}

	friend ostream& operator<<(ostream& stream, Vertex& V)
	{
		for (auto it : V.adjList)
			stream << it++ << " ";
		stream << endl;

		return stream;
	}
};


class Graph						// Класс графа, в котором реализованы вспомогательные функции
{
protected:
	vector<Vertex> adjMatrix;	// Матрица смежности (массив из вершин)
	int vertices = 0;			// кол-во вершин

	inline void initialize(vector<Vertex>& matrix, int Vertices)
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
		if (index - vertices >= 0)
			throw IndexOutOfBoundsException(vertices, index);
		return adjMatrix[index];
	}


	friend ostream& operator<<(ostream& stream, Graph& G)
	{
		stream << G.vertices << endl;
		
		vector<Vertex>::iterator it = G.adjMatrix.begin();
		while (it - G.adjMatrix.end() != 0)
			stream << *it++;

		return stream;
	}

	friend istream& operator>>(istream& stream, Graph& G)
	{
		if (typeid(stream) == typeid(ifstream))
		{
			int vertex;
			stream >> vertex;

			if (vertex - G.vertices != 0)
				G.initialize(G.adjMatrix, vertex);
		}
		for (int i = 0; i - G.vertices < 0; i++)
			for (int j = 0; j - G.vertices < 0; j++)
				stream >> G.adjMatrix[i][j];

		return stream;
	}
};



class FlowGraph : public Graph		// Класс графа, предназначенного для поиска максимального потока
{
protected:
	bool bfs(vector<Vertex>& rGraph, int source, int target, int parent[])
	{
		// BFS - поиск возможного минимального пути в сети
		bool* visited = new bool[vertices];
		bool* b = &visited[0]; bool* e = &visited[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = false;
			*e-- = false;
		}

		queue<int> q;

		q.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v - vertices < 0; v += 4)
			{
				if (v - vertices + 3 < 0)
				{
					if (!visited[v] && rGraph[u][v] > 0)
					{
						if (v - target == 0)
						{
							parent[v] = u;
							delete[] visited;
							return true;
						}
						q.push(v);
						parent[v] = u;
						visited[v] = true;
					}

					if (!visited[v + 1] && rGraph[u][v + 1] > 0)
					{
						if (v + 1 - target == 0)
						{
							parent[v + 1] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 1);
						parent[v + 1] = u;
						visited[v + 1] = true;
					}

					if (!visited[v + 2] && rGraph[u][v + 2] > 0)
					{
						if (v + 2 - target == 0)
						{
							parent[v + 2] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 2);
						parent[v + 2] = u;
						visited[v + 2] = true;
					}

					if (!visited[v + 3] && rGraph[u][v + 3] > 0)
					{
						if (v + 3 - target == 0)
						{
							parent[v + 3] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 3);
						parent[v + 3] = u;
						visited[v + 3] = true;
					}
				}
				else if (v - vertices + 2 < 0)
				{
					if (!visited[v] && rGraph[u][v] > 0)
					{
						if (v - target == 0)
						{
							parent[v] = u;
							delete[] visited;
							return true;
						}
						q.push(v);
						parent[v] = u;
						visited[v] = true;
					}

					if (!visited[v + 1] && rGraph[u][v + 1] > 0)
					{
						if (v + 1 - target == 0)
						{
							parent[v + 1] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 1);
						parent[v + 1] = u;
						visited[v + 1] = true;
					}

					if (!visited[v + 2] && rGraph[u][v + 2] > 0)
					{
						if (v + 2 - target == 0)
						{
							parent[v + 2] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 2);
						parent[v + 2] = u;
						visited[v + 2] = true;
					}
				}
				else if (v - vertices + 1 < 0)
				{
					if (!visited[v] && rGraph[u][v] > 0)
					{
						if (v - target == 0)
						{
							parent[v] = u;
							delete[] visited;
							return true;
						}
						q.push(v);
						parent[v] = u;
						visited[v] = true;
					}

					if (!visited[v + 1] && rGraph[u][v + 1] > 0)
					{
						if (v + 1 - target == 0)
						{
							parent[v + 1] = u;
							delete[] visited;
							return true;
						}
						q.push(v + 1);
						parent[v + 1] = u;
						visited[v + 1] = true;
					}
				}
				else
				{
					if (!visited[v] && rGraph[u][v] > 0)
					{
						if (v - target == 0)
						{
							parent[v] = u;
							delete[] visited;
							return true;
						}
						q.push(v);
						parent[v] = u;
						visited[v] = true;
					}
				}
			}
		}
		delete[] visited;
		return false;
	}

	bool dfs(vector<Vertex>& rGraph, int source, int target, int parent[])
	{
		// DFS; вместо очереди используется стек
		bool* visited = new bool[vertices];
		bool* b = &visited[0]; bool* e = &visited[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = false;
			*e-- = false;
		}

		stack<int> S;

		S.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!S.empty())
		{
			int u = S.top();
			S.pop();

			for (int v = 0; v - vertices < 0; v += 4)
			{
				if (v - vertices + 3 < 0)
				{
					visited[v] = (!visited[v]) ? ((rGraph[u][v] > 0) ? (S.push(v), parent[v] = u, true) : false) : true;
					visited[v + 1] = (!visited[v + 1]) ? ((rGraph[u][v + 1] > 0) ? (S.push(v + 1), parent[v + 1] = u, true) : false) : true;
					visited[v + 2] = (!visited[v + 2]) ? ((rGraph[u][v + 2] > 0) ? (S.push(v + 2), parent[v + 2] = u, true) : false) : true;
					visited[v + 3] = (!visited[v + 3]) ? ((rGraph[u][v + 3] > 0) ? (S.push(v + 3), parent[v + 3] = u, true) : false) : true;
				}
				else if (v - vertices + 2 < 0)
				{
					visited[v] = (!visited[v]) ? ((rGraph[u][v] > 0) ? (S.push(v), parent[v] = u, true) : false) : true;
					visited[v + 1] = (!visited[v + 1]) ? ((rGraph[u][v + 1] > 0) ? (S.push(v + 1), parent[v + 1] = u, true) : false) : true;
					visited[v + 2] = (!visited[v + 2]) ? ((rGraph[u][v + 2] > 0) ? (S.push(v + 2), parent[v + 2] = u, true) : false) : true;
				}
				else if (v - vertices + 1 < 0)
				{
					visited[v] = (!visited[v]) ? ((rGraph[u][v] > 0) ? (S.push(v), parent[v] = u, true) : false) : true;
					visited[v + 1] = (!visited[v + 1]) ? ((rGraph[u][v + 1] > 0) ? (S.push(v + 1), parent[v + 1] = u, true) : false) : true;
				}
				else
					visited[v] = (!visited[v]) ? ((rGraph[u][v] > 0) ? (S.push(v), parent[v] = u, true) : false) : true;
			}
		}
		bool res = visited[target];
		delete[] visited;
		return res;
	}

	bool bfsForDinic(vector<Vertex>& rGraph, int source, int target, int level[])
	{
		// BFS + создание слоистой вспомогательной сети level
		queue<int> q;
		q.push(source);
		int* b = &level[0]; int* e = &level[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = -1;
			*e-- = -1;
		}
		
		level[source] = 0;

		while (!q.empty())
		{
			int u = q.front();
			q.pop();

			for (int v = 0; v - vertices < 0; v += 4)
			{
				if (v - vertices + 3 < 0)
				{
					level[v] = (level[v] >= 0) ? level[v] : ((rGraph[u][v] > 0) ? (q.push(v), level[u] + 1) : level[v]);
					level[v + 1] = (level[v + 1] >= 0) ? level[v + 1] : ((rGraph[u][v + 1] > 0) ? (q.push(v + 1), level[u] + 1) : level[v + 1]);
					level[v + 2] = (level[v + 2] >= 0) ? level[v + 2] : ((rGraph[u][v + 2] > 0) ? (q.push(v + 2), level[u] + 1) : level[v + 2]);
					level[v + 3] = (level[v + 3] >= 0) ? level[v + 3] : ((rGraph[u][v + 3] > 0) ? (q.push(v + 3), level[u] + 1) : level[v + 3]);
				}
				else if (v - vertices + 2 < 0)
				{
					level[v] = (level[v] >= 0) ? level[v] : ((rGraph[u][v] > 0) ? (q.push(v), level[u] + 1) : level[v]);
					level[v + 1] = (level[v + 1] >= 0) ? level[v + 1] : ((rGraph[u][v + 1] > 0) ? (q.push(v + 1), level[u] + 1) : level[v + 1]);
					level[v + 2] = (level[v + 2] >= 0) ? level[v + 2] : ((rGraph[u][v + 2] > 0) ? (q.push(v + 2), level[u] + 1) : level[v + 2]);
				}
				else if (v - vertices + 1 < 0)
				{
					level[v] = (level[v] >= 0) ? level[v] : ((rGraph[u][v] > 0) ? (q.push(v), level[u] + 1) : level[v]);
					level[v + 1] = (level[v + 1] >= 0) ? level[v + 1] : ((rGraph[u][v + 1] > 0) ? (q.push(v + 1), level[u] + 1) : level[v + 1]);
				}
				else
					level[v] = (level[v] >= 0) ? level[v] : ((rGraph[u][v] > 0) ? (q.push(v), level[u] + 1) : level[v]);
			}
		}
		return level[target] >= 0;
	}

	int dfsForDinic(vector<Vertex>& rGraph, int u, int target, int flow, int level[], int start[])
	{
		if (u - target == 0) return flow;

		for (int& v = start[u]; v - vertices < 0; v += 4)
		{
			if (v - vertices + 3 < 0)
			{
				if (rGraph[u][v] > 0 && level[v] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v, target, ((flow < rGraph[u][v]) ? flow : rGraph[u][v]), level, start))
					{
						rGraph[u][v] -= currFlow;
						rGraph[v][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 1] > 0 && level[v + 1] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 1, target, ((flow < rGraph[u][v + 1]) ? flow : rGraph[u][v + 1]), level, start))
					{
						rGraph[u][v + 1] -= currFlow;
						rGraph[v + 1][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 2] > 0 && level[v + 2] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 2, target, ((flow < rGraph[u][v + 2]) ? flow : rGraph[u][v + 2]), level, start))
					{
						rGraph[u][v + 2] -= currFlow;
						rGraph[v + 2][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 3] > 0 && level[v + 3] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 3, target, ((flow < rGraph[u][v + 3]) ? flow : rGraph[u][v + 3]), level, start))
					{
						rGraph[u][v + 3] -= currFlow;
						rGraph[v + 3][u] += currFlow;
						return currFlow;
					}
				}
			}
			else if (v - vertices + 2 < 0)
			{
				if (rGraph[u][v] > 0 && level[v] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v, target, ((flow < rGraph[u][v]) ? flow : rGraph[u][v]), level, start))
					{
						rGraph[u][v] -= currFlow;
						rGraph[v][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 1] > 0 && level[v + 1] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 1, target, ((flow < rGraph[u][v + 1]) ? flow : rGraph[u][v + 1]), level, start))
					{
						rGraph[u][v + 1] -= currFlow;
						rGraph[v + 1][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 2] > 0 && level[v + 2] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 2, target, ((flow < rGraph[u][v + 2]) ? flow : rGraph[u][v + 2]), level, start))
					{
						rGraph[u][v + 2] -= currFlow;
						rGraph[v + 2][u] += currFlow;
						return currFlow;
					}
				}
			}
			else if (v - vertices + 1 < 0)
			{
				if (rGraph[u][v] > 0 && level[v] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v, target, ((flow < rGraph[u][v]) ? flow : rGraph[u][v]), level, start))
					{
						rGraph[u][v] -= currFlow;
						rGraph[v][u] += currFlow;
						return currFlow;
					}
				}

				if (rGraph[u][v + 1] > 0 && level[v + 1] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v + 1, target, ((flow < rGraph[u][v + 1]) ? flow : rGraph[u][v + 1]), level, start))
					{
						rGraph[u][v + 1] -= currFlow;
						rGraph[v + 1][u] += currFlow;
						return currFlow;
					}
				}
			}
			else
			{
				if (rGraph[u][v] > 0 && level[v] - level[u] - 1 == 0)
				{
					if (int currFlow = dfsForDinic(rGraph, v, target, ((flow < rGraph[u][v]) ? flow : rGraph[u][v]), level, start))
					{
						rGraph[u][v] -= currFlow;
						rGraph[v][u] += currFlow;
						return currFlow;
					}
				}
			}
		}
		return 0;
	}

	int fordFulkerson(int source, int target)
	{
		// Алгоритм Форда-Фалкерсона
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		vector<Vertex>::iterator it1, it2;
		it1 = rGraph.begin(); it2 = adjMatrix.begin();
		while (it2 - adjMatrix.end() != 0)
			*it1++ = *it2++;

		int* parent = new int[vertices];
		int* b = &parent[0]; int* e = &parent[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = 0;
			*e-- = 0;
		}

		int maxFlow = 0;

		while (dfs(rGraph, source, target, parent))
		{
			int pathFlow = INT_MAX;
			for (int v = target; v - source != 0; v = parent[v])
				pathFlow = (pathFlow < rGraph[parent[v]][v]) ? pathFlow : rGraph[parent[v]][v];

			for (int v = target; v - source != 0; v = parent[v])
			{
				rGraph[parent[v]][v] -= pathFlow;
				rGraph[v][parent[v]] += pathFlow;
			}

			maxFlow += pathFlow;
		}
		delete[] parent;
		return maxFlow;
	}

	int edmondsKarp(int source, int target)
	{
		// Алгоритм Эдмондса-Карпа
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		vector<Vertex>::iterator it1, it2;
		it1 = rGraph.begin(); it2 = adjMatrix.begin();
		while (it2 - adjMatrix.end() != 0)
			*it1++ = *it2++;

		int* parent = new int[vertices];
		int* b = &parent[0];
		int* e = &parent[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = 0;
			*e-- = 0;
		}

		int maxFlow = 0;

		while (bfs(rGraph, source, target, parent))
		{
			int pathFlow = INT_MAX;
			for (int v = target; v - source != 0; v = parent[v])
				pathFlow = (pathFlow < rGraph[parent[v]][v]) ? pathFlow : rGraph[parent[v]][v];

			for (int v = target; v - source != 0; v = parent[v])
			{
				rGraph[parent[v]][v] -= pathFlow;
				rGraph[v][parent[v]] += pathFlow;
			}

			maxFlow += pathFlow;
		}
		delete[] parent;
		return maxFlow;
	}

	int dinic(int source, int target)
	{
		// Алгоритм Диница
		vector<Vertex> rGraph(vertices, Vertex(vertices));
		vector<Vertex>::iterator it1, it2;
		it1 = rGraph.begin(); it2 = adjMatrix.begin();
		while (it2 - adjMatrix.end() != 0)
			*it1++ = *it2++;

		int maxFlow = 0;
		
		int* level = new int[vertices];
		int* b = &level[0]; int* e = &level[vertices - 1];
		while (b - e <= 0)
		{
			*b++ = -1;
			*e-- = -1;
		}

		int* start = new int[vertices];
		while (bfsForDinic(rGraph, source, target, level))
		{
			b = &start[0]; e = &start[vertices - 1];
			while (b - e <= 0)
			{
				*b++ = 0;
				*e-- = 0;
			}
			
			while (int flow = dfsForDinic(rGraph, source, target, INT_MAX, level, start))
				maxFlow += flow;
		}
		delete[] level; delete[] start;
		return maxFlow;
	}
public:
	FlowGraph() : Graph() {}
	FlowGraph(int V) : Graph(V) {}

	int findMaxFlow(int source, int target, int algo)
	{
		if (target > vertices) throw InvalidVertexException(target);

		switch (algo)
		{
		case 1:
			return fordFulkerson(source, target);
		case 2:
			return edmondsKarp(source, target);
		case 3:
			return dinic(source, target);
		default:
			throw InvalidAlgorithmException(algo);
		}
	}
};



int main()
{
	try
	{
		/*FlowGraph G(4);
		ifstream fin; fin.open("test0.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}*/

		/*FlowGraph G(10);
		ifstream fin; fin.open("test1.txt");
		if (fin)
		{
			fin >> G;
			fin.close();
		}*/

		/*FlowGraph G(12);
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

		//Измерение алгоритма Форда-Фалкерсона
		int ff_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ff_flow = G.findMaxFlow(SOURCE, TARGET, 1);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ff_total += (timer_2 - timer_1);
		for (i = 0; i - 9 < 0; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(SOURCE, TARGET, 1);
			timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			ff_total += (timer_2 - timer_1);
		}
		cout << "[Ford-Fulkerson algorithm]\tFlow: " << ff_flow << endl;
		cout << "[Ford-Fulkerson algorithm]\tTotal time: " << ff_total << " ns" << endl;


		//Измерение алгоритма Эдмондса-Карпа
		int ek_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ek_flow = G.findMaxFlow(SOURCE, TARGET, 2);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		ek_total += (timer_2 - timer_1);
		for (i = 0; i - 9 < 0; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(SOURCE, TARGET, 2);
			timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			ek_total += (timer_2 - timer_1);
		}
		cout << "[Edmonds-Karp algorithm]\tFlow: " << ek_flow << endl;
		cout << "[Edmonds-Karp algorithm]\tTotal time: " << ek_total << " ns" << endl;


		//Измерение алгоритма Диница
		int d_flow;
		timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		d_flow = G.findMaxFlow(SOURCE, TARGET, 3);
		timer_2 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
		d_total += (timer_2 - timer_1);
		for (i = 0; i - 9 < 0; i++)
		{
			timer_1 = chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();
			G.findMaxFlow(SOURCE, TARGET, 3);
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
