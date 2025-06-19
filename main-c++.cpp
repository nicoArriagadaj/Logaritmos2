#include <iostream>
#include <vector>    // arreglos
#include <algorithm> // heap y sort
#include <cmath>     // Para sqrt y pow
#include <cstdint>   // Para int64_t
#include <queue>

#include <random> // para aleatoriedad en las pruebas 
#include <chrono> // para usar std::chrono

/* 
 * Calcula el cuadrado de la distancia euclidiana entre dos puntos en [0,1]x[0,1]
 * Parámetros:
 *      x1, y1 - coordenadas del primer punto
 *      x2, y2 - coordenadas del segundo punto
 * Retorna:
 *      double - cuadrado de la distancia euclidiana o -1.0 si hay error en el rango
 */
double Peso_arista__cuadrado_distancia_euclidiana(double x1, double x2, double y1, double y2) {
    if ((x1 < 0.0 || x1 > 1.0) ||
        (x2 < 0.0 || x2 > 1.0) ||
        (y1 < 0.0 || y1 > 1.0) ||
        (y2 < 0.0 || y2 > 1.0)) {
        std::cerr << "Error: Los valores deben estar entre 0 y 1 Incluyendolosss" << std::endl;
        return -1.0;
    }
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dist2 = dx * dx + dy * dy;
    return dist2;
}

/* 
 * Estructura que representa un nodo con coordenadas (x, y) en [0,1]x[0,1]
 * Campos:
 *      double x - coordenada x
 *      double y - coordenada y
 */
struct Nodo {
    double x, y; // Coordenadas en [0, 1]
};

/* 
 * Estructura que representa una arista entre dos nodos y su peso
 * Campos:
 *      int u - índice primer nodo
 *      int v - índice segundo nodo
 *      double peso - peso de la arista
 */
struct Arista {
    int u, v;
    double peso;
    // Sobrecarga para usar en heaps de prioridad mínima (por peso)
    bool operator<(const Arista& otra) const {
        return peso > otra.peso; 
    }
};

/*
 * Estructura Union-Find con compresión de caminos y unión por tamaño.
 * Métodos:
 *   UnionFind(int n) - Inicializa n conjuntos.
 *   int find(int x) - Retorna la raíz del conjunto de x, con compresión de caminos.
 *   void Union(int x, int y) - Une los conjuntos de x e y.
 *   bool connected(int x, int y) - Indica si x e y están en el mismo conjunto.
 */
struct UnionFind {
    std::vector<int> parent, size;
    UnionFind(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }
    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
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

/*
 * Genera un vector de N nodos, cada uno con coordenadas (x, y) aleatorias en [0, 1].
 * Parámetros:
 *   N - Número de nodos a generar.
 * Retorno:
 *   std::vector<Nodo> - Vector con los nodos generados.
 */
std::vector<Nodo> generar_nodos(int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    std::vector<Nodo> nodos;
    for (int i = 0; i < N; ++i)
        nodos.push_back({dis(gen), dis(gen)});
    return nodos;
}

/*
 * Genera todas las posibles aristas entre los nodos dados y calcula el peso de cada una.
 * Parámetros:
 *   nodos - Vector de nodos.
 * Retorno:
 *   std::vector<Arista> - Vector con todas las aristas posibles.
 */
std::vector<Arista> generar_aristas(const std::vector<Nodo>& nodos) {
    int n = nodos.size();
    std::vector<Arista> aristas;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double peso = Peso_arista__cuadrado_distancia_euclidiana(nodos[i].x, nodos[j].x, nodos[i].y, nodos[j].y);
            aristas.push_back({i, j, peso});
        }
    }
    return aristas;
}

/*
 * Algoritmo de Kruskal usando sort y Union-Find con compresión de caminos.
 * Parámetros:
 *   n - Cantidad de nodos.
 *   aristas - Vector de aristas del grafo.
 * Retorno:
 *   std::vector<Arista> - Aristas que forman el MST.
 */
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

/*
 * Algoritmo de Kruskal usando heap (priority_queue) y Union-Find con compresión de caminos.
 * Parámetros:
 *   n - Cantidad de nodos.
 *   aristas - Vector de aristas del grafo.
 * Retorno:
 *   std::vector<Arista> - Aristas que forman el MST.
 */
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

/*
 * Estructura Union-Find sin compresión de caminos (solo unión por tamaño).
 * Métodos:
 *   UnionFindSinCompresion(int n) - Inicializa n conjuntos.
 *   int find(int x) - Retorna la raíz del conjunto de x (sin compresión).
 *   void Union(int x, int y) - Une los conjuntos de x e y.
 *   bool connected(int x, int y) - Indica si x e y están en el mismo conjunto.
 */
struct UnionFindSinCompresion {
    std::vector<int> parent, size;
    UnionFindSinCompresion(int n) : parent(n), size(n, 1) {
        for (int i = 0; i < n; ++i)
            parent[i] = i;
    }
    int find(int x) {
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

/*
 * Algoritmo de Kruskal usando sort y Union-Find sin compresión de caminos.
 * Parámetros:
 *   n - Cantidad de nodos.
 *   aristas - Vector de aristas del grafo.
 * Retorno:
 *   std::vector<Arista> - Aristas que forman el MST.
 */
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

/*
 * Algoritmo de Kruskal usando heap y Union-Find sin compresión de caminos.
 * Parámetros:
 *   n - Cantidad de nodos.
 *   aristas - Vector de aristas del grafo.
 * Retorno:
 *   std::vector<Arista> - Aristas que forman el MST.
 */
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


#include <fstream> // Para archivo csv
#include <iomanip> // para más precisión.

int main() {
    std::vector<int> Nval = {32, 64, 128, 256, 512, 1024, 2048, 4096}; 
    int repeticiones = 5;

    // Abrir archivo para guardar resultados
    std::ofstream fout("resultados.csv");
    fout << "N,metodo,rep,t_ms,peso_mst\n";

    // Configura la precisión para todo el archivo (opcional, para CSV)
    fout << std::fixed << std::setprecision(8);
    std::cout << std::fixed << std::setprecision(8);

    for (int N : Nval) {
        std::cout << "---- N = " << N << " ----" << std::endl;

        // Vectores para guardar resultados individuales
        std::vector<double> tiempos_sort, tiempos_heap, tiempos_sort_nc, tiempos_heap_nc;
        std::vector<double> pesos_sort, pesos_heap, pesos_sort_nc, pesos_heap_nc;

        for (int rep = 0; rep < repeticiones; ++rep) {
            auto nodos = generar_nodos(N);
            auto aristas = generar_aristas(nodos);

            // Kruskal + sort + compresión
            auto start_sort = std::chrono::high_resolution_clock::now();
            auto mst_sort = kruskal_sort_Compresion_caminos(N, aristas);
            auto end_sort = std::chrono::high_resolution_clock::now();

            double tiempo_sort = std::chrono::duration<double, std::milli>(end_sort - start_sort).count();
            double peso_sort = 0;
            for (const auto& a : mst_sort) peso_sort += a.peso;
            tiempos_sort.push_back(tiempo_sort);
            pesos_sort.push_back(peso_sort);
            fout << N << ",sort_compresion," << rep+1 << "," << tiempo_sort << "," << peso_sort << "\n";

            // Kruskal + heap + compresión
            auto start_heap = std::chrono::high_resolution_clock::now();
            auto mst_heap = kruskal_heap_Compresion_caminos(N, aristas);
            auto end_heap = std::chrono::high_resolution_clock::now();

            double tiempo_heap = std::chrono::duration<double, std::milli>(end_heap - start_heap).count();
            double peso_heap = 0;
            for (const auto& a : mst_heap) peso_heap += a.peso;
            tiempos_heap.push_back(tiempo_heap);
            pesos_heap.push_back(peso_heap);
            fout << N << ",heap_compresion," << rep+1 << "," << tiempo_heap << "," << peso_heap << "\n";

            // Kruskal + sort + SIN compresión
            auto start_sort_nc = std::chrono::high_resolution_clock::now();
            auto mst_sort_nc = kruskal_sort_SinCompresion_caminos(N, aristas);
            auto end_sort_nc = std::chrono::high_resolution_clock::now();

            double tiempo_sort_nc = std::chrono::duration<double, std::milli>(end_sort_nc - start_sort_nc).count();
            double peso_sort_nc = 0;
            for (const auto& a : mst_sort_nc) peso_sort_nc += a.peso;
            tiempos_sort_nc.push_back(tiempo_sort_nc);
            pesos_sort_nc.push_back(peso_sort_nc);
            fout << N << ",sort_sincompresion," << rep+1 << "," << tiempo_sort_nc << "," << peso_sort_nc << "\n";

            // Kruskal + heap + SIN compresión
            auto start_heap_nc = std::chrono::high_resolution_clock::now();
            auto mst_heap_nc = kruskal_heap_SinCompresion_caminos(N, aristas);
            auto end_heap_nc = std::chrono::high_resolution_clock::now();

            double tiempo_heap_nc = std::chrono::duration<double, std::milli>(end_heap_nc - start_heap_nc).count();
            double peso_heap_nc = 0;
            for (const auto& a : mst_heap_nc) peso_heap_nc += a.peso;
            tiempos_heap_nc.push_back(tiempo_heap_nc);
            pesos_heap_nc.push_back(peso_heap_nc);
            fout << N << ",heap_sincompresion," << rep+1 << "," << tiempo_heap_nc << "," << peso_heap_nc << "\n";
        }

        // Promedios
        double avg_sort = std::accumulate(tiempos_sort.begin(), tiempos_sort.end(), 0.0) / repeticiones;
        double avg_heap = std::accumulate(tiempos_heap.begin(), tiempos_heap.end(), 0.0) / repeticiones;
        double avg_sort_nc = std::accumulate(tiempos_sort_nc.begin(), tiempos_sort_nc.end(), 0.0) / repeticiones;
        double avg_heap_nc = std::accumulate(tiempos_heap_nc.begin(), tiempos_heap_nc.end(), 0.0) / repeticiones;

        double avg_peso_sort = std::accumulate(pesos_sort.begin(), pesos_sort.end(), 0.0) / repeticiones;
        double avg_peso_heap = std::accumulate(pesos_heap.begin(), pesos_heap.end(), 0.0) / repeticiones;
        double avg_peso_sort_nc = std::accumulate(pesos_sort_nc.begin(), pesos_sort_nc.end(), 0.0) / repeticiones;
        double avg_peso_heap_nc = std::accumulate(pesos_heap_nc.begin(), pesos_heap_nc.end(), 0.0) / repeticiones;

        std::cout << "[Kruskal + sort + compresion]    Promedio tiempo: " 
            << avg_sort << "ms, Promedio peso: " << avg_peso_sort << std::endl;

        std::cout << "[Kruskal + heap + compresion]    Promedio tiempo: " 
            << avg_heap << "ms, Promedio peso: " << avg_peso_heap << std::endl;

        std::cout << "[Kruskal + sort + sin compresion]    Promedio tiempo: " 
            << avg_sort_nc << "ms, Promedio peso: " << avg_peso_sort_nc << std::endl;

        std::cout << "[Kruskal + heap + sin compresion]    Promedio tiempo: " 
            << avg_heap_nc << "ms, Promedio peso: " << avg_peso_heap_nc << std::endl;
    }

    std::cout << "Resultados individuales y promedios guardados en resultados.csv" << std::endl;
    return 0;
}