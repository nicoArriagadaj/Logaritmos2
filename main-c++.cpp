#include <iostream>
#include <vector>    // arreglos
#include <algorithm> // heap y sort
#include <cmath>     // Para sqrt y pow
#include <cstdint>   // Para int64_t
#include <queue>

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

    // Para priority_queue, necesitas sobrecargar el operador <
    // Para min-heap (menor peso tiene más prioridad)
    bool operator<(const Arista& otra) const {
        return peso > otra.peso; // Cambia a < para max-heap
    }
};

// ejempl de uso MAIN;
int main() {
    // Vector de nodos de ejemplo
    std::vector<Nodo> nodos = { {0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0} };

    // Vector de aristas PARA USAR SORT
    std::vector<Arista> aristas;
    // Calcula todas las aristas posibles entre los nodos (sin repetir ni invertir (0,1), (0,2), (1,2))
    // push_back añade elemento al final. es como append en python
    // El vector mantiene el orden de inserción 
    for (int i = 0; i < nodos.size(); ++i) {
        for (int j = i + 1; j < nodos.size(); ++j) {
            double peso = Peso_arista__cuadrado_distancia_euclidiana(
                nodos[i].x, nodos[j].x, nodos[i].y, nodos[j].y);
            aristas.push_back({i, j, peso});
        }
    }
    // Ordenar el vector de menor a mayor peso CON SORT
    std::sort(aristas.begin(), aristas.end(), [](const Arista& a, const Arista& b){
        return a.peso < b.peso;
    });

    std::cout << "Aristas ordenadassss: ";
    for (auto& a : aristas) {
        std::cout << "[(" << a.u << "," << a.v << "): " << a.peso << "] ";
    }
    std::cout << std::endl;

    // Heap de aristas (menor peso al tope)
    std::priority_queue<Arista> heap;
    //  agrega el elemento al heap y automáticamente lo reubica internamente para mantener la propiedad de heap.
    for (auto& a : aristas) {
        heap.push(a);
    }

    std::cout << "Aristas desde el heap: ";
    while (!heap.empty()) {
        auto a = heap.top();
        std::cout << "[(" << a.u << "," << a.v << "): " << a.peso << "] ";
        heap.pop();
    }
    std::cout << std::endl;
    return 0;
}
