#include <iostream>
#include <vector>    // arreglos
#include <algorithm> // heap y sort
#include <cmath>     // Para sqrt y pow
#include <cstdint>   // Para int64_t
#include <queue>

#include <random> // para aleatoriedad en las pruebas 
#include <chrono> // para usar std::chrono

// Peso = cuadrado de la distancia euclidiana entre dos puntos reales en [0,1]
/* peso de arista = distancia_euclidiana * distancia_euclidiana, 
double es de 64 bits, es un REAL, por esto no int64_t, esto es un entero */
double Peso_arista__cuadrado_distancia_euclidiana(double x1, double x2, double y1, double y2) {
    // verifica la condicion
    if ((x1 < 0.0 || x1 > 1.0) ||
        (x2 < 0.0 || x2 > 1.0) ||
        (y1 < 0.0 || y1 > 1.0) ||
        (y2 < 0.0 || y2 > 1.0)) {
        std::cerr << "Error: Los valores deben estar entre 0 y 1 Incluyendolosss" << std::endl;
        return -1.0;
    }

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dist2 = dx * dx + dy * dy; // aqui, tenemos el cuadrado de la distancia
    return dist2;
}

// nodo representado por coordenadas en [0, 1]
struct Nodo {
    double x, y; // Coordenadas en [0, 1]
};

// arista con índices de los nodos que conecta
struct Arista {
    int u, v;    // índices de los nodos en el vector de nodos
    double peso; // peso de la arista
    // para Minimo
    bool operator<(const Arista& otra) const {
        return peso > otra.peso; 
    }
};

// estructura UnionFind con compresion_caminos
struct UnionFind {
    // arreglo, padre del nodo y tamaño
    std::vector<int> parent, size;
    // Inicializa n conjuntos (de 0 a n-1)
    UnionFind(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }
    // Encuentra la raíz (con compresión de caminos)
    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]); // COMPRESION DE CAMINOS
        return parent[x];
    }
    // Une los conjuntos de x e y (union by size)
    void Union(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        if (rx == ry) return; // ya están unidos
        if (size[rx] < size[ry]) std::swap(rx, ry);
        parent[ry] = rx;
        size[rx] += size[ry];
    }
    // VE SI ESTÁN EN EL MISMO CAMINO, ES DECIR, find(x) = find(y)
    bool connected(int x, int y) {
        return find(x) == find(y);
    }
};


// Genera N ndos aleatorios en el rango 
std::vector<Nodo> generar_nodos(int N) {
    // genera rango aleatorio 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    // inicializa nodos 
    std::vector<Nodo> nodos;
    for (int i = 0; i < N; ++i)
        // les hace push_back -> los envia al final.
        nodos.push_back({dis(gen), dis(gen)});
    return nodos;
}


// Genera todas las aristas posibles entre los nodos
std::vector<Arista> generar_aristas(const std::vector<Nodo>& nodos) {
    int n = nodos.size();
    // crea arreglo de aristas
    std::vector<Arista> aristas;
    // genera las aristas, doble ciclo for, para que, cree las combinaciones posibles. 
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // calcula el peso de la arista.
            double peso = Peso_arista__cuadrado_distancia_euclidiana(nodos[i].x, nodos[j].x, nodos[i].y, nodos[j].y);
            // push_back es como append, de python, mantiene la estructura del vector(o del arreglo) y añade el elemento al final.
            aristas.push_back({i, j, peso});
        }
    }
    return aristas;
}

// Kruskal usando arreglo ordenado
std::vector<Arista> kruskal_sort_Compresion_caminos(int n, std::vector<Arista> aristas) {
    std::sort(aristas.begin(), aristas.end(), [](const Arista& a, const Arista& b) {
        return a.peso < b.peso;
    });
    UnionFind uf(n);
    std::vector<Arista> mst;
    for (const auto& a : aristas) {
        if (!uf.connected(a.u, a.v)) {
            uf.Union(a.u, a.v);
            mst.push_back(a);
        }
    }
    return mst;
}

// Kruskal usando heap clásico
std::vector<Arista> kruskal_heap_Compresion_caminos(int n, std::vector<Arista> aristas) {
    std::priority_queue<Arista> heap;
    for (const auto& a : aristas) heap.push(a);
    UnionFind uf(n);
    std::vector<Arista> mst;
    while (!heap.empty() && mst.size() < n - 1) {
        Arista a = heap.top();
        heap.pop();
        if (!uf.connected(a.u, a.v)) {
            uf.Union(a.u, a.v);
            mst.push_back(a);
        }
    }
    return mst;
}



// estructura UnionFind con compresion_caminos
struct UnionFindSinCompresion {
    std::vector<int> parent, size;
    UnionFindSinCompresion(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }
    // Encuentra la raíz SIN COMPRESION DE CAMINOS 
    int find(int x) {
        // entonces, hace ciclo while, hasta encontrar la raíz.
        while (parent[x] != x)
            x = parent[x];
        return x;
    }
    void Union(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        if (rx == ry) return;
        if (size[rx] < size[ry]) std::swap(rx, ry);
        parent[ry] = rx;
        size[rx] += size[ry];
    }
    bool connected(int x, int y) {
        return find(x) == find(y);
    }
};

//  kruskal_sort sin compresion de caminos 
std::vector<Arista> kruskal_sort_SinCompresion_caminos(int n, std::vector<Arista> aristas) {
    std::sort(aristas.begin(), aristas.end(), [](const Arista& a, const Arista& b) {
        return a.peso < b.peso;
    });
    UnionFindSinCompresion uf(n);
    std::vector<Arista> mst;
    for (const auto& a : aristas) {
        if (!uf.connected(a.u, a.v)) {
            uf.Union(a.u, a.v);
            mst.push_back(a);
        }
    }
    return mst;
}


// kruskal_heap_SinCompresion_de_caminos
std::vector<Arista> kruskal_heap_SinCompresion_caminos(int n, std::vector<Arista> aristas) {
    std::priority_queue<Arista> heap;
    for (const auto& a : aristas) heap.push(a);
    UnionFindSinCompresion uf(n);
    std::vector<Arista> mst;
    while (!heap.empty() && mst.size() < n - 1) {
        Arista a = heap.top();
        heap.pop();
        if (!uf.connected(a.u, a.v)) {
            uf.Union(a.u, a.v);
            mst.push_back(a);
        }
    }
    return mst;
}





int main() {
    int N = 10000; // total de nodos.
    auto nodos = generar_nodos(N);
    auto aristas = generar_aristas(nodos);
    // Compresión de caminos sort
    auto start_sort = std::chrono::high_resolution_clock::now();
    auto mst_sort = kruskal_sort_Compresion_caminos(N, aristas);
    auto end_sort = std::chrono::high_resolution_clock::now();
    // Compresión de caminos heap
    auto start_heap = std::chrono::high_resolution_clock::now();
    auto mst_heap = kruskal_heap_Compresion_caminos(N, aristas);
    auto end_heap = std::chrono::high_resolution_clock::now();
    // Sin compresión de caminos sort
    auto start_sort_nc = std::chrono::high_resolution_clock::now();
    auto mst_sort_nc = kruskal_sort_SinCompresion_caminos(N, aristas);
    auto end_sort_nc = std::chrono::high_resolution_clock::now();
    // Sin compresión de caminos heap
    auto start_heap_nc = std::chrono::high_resolution_clock::now();
    auto mst_heap_nc = kruskal_heap_SinCompresion_caminos(N, aristas);
    auto end_heap_nc = std::chrono::high_resolution_clock::now();

    // calculo de chrono -> tiempo de ejecucion sort_heap con compresion caminos
    std::chrono::duration<double, std::milli> tiempo_sort = end_sort - start_sort;
    std::chrono::duration<double, std::milli> tiempo_heap = end_heap - start_heap;


    // calculo de chrono -> tiempo de ejecucion sort_heap sin compresion caminos (nc)
    std::chrono::duration<double, std::milli> tiempo_sort_nc = end_sort_nc - start_sort_nc;
    std::chrono::duration<double, std::milli> tiempo_heap_nc = end_heap_nc - start_heap_nc;

    std::cout << "[Kruskal + sort + compresion] Tiempo: " << tiempo_sort.count() << "ms, Peso total: ";
    double peso_total = 0;
    for (const auto& a : mst_sort) peso_total += a.peso;
    std::cout << peso_total << std::endl;

    std::cout << "[Kruskal + heap + compresion] Tiempo: " << tiempo_heap.count() << "ms, Peso total: ";
    peso_total = 0;
    for (const auto& a : mst_heap) peso_total += a.peso;
    std::cout << peso_total << std::endl;

    std::cout << "[Kruskal + sort + sin compresion] Tiempo: " << tiempo_sort_nc.count() << "ms, Peso total: ";
    peso_total = 0;
    for (const auto& a : mst_sort_nc) peso_total += a.peso;
    std::cout << peso_total << std::endl;

    std::cout << "[Kruskal + heap + sin compresion] Tiempo: " << tiempo_heap_nc.count() << "ms, Peso total: ";
    peso_total = 0;
    for (const auto& a : mst_heap_nc) peso_total += a.peso;
    std::cout << peso_total << std::endl;

    return 0;
}
