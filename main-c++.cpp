#include <iostream>
#include <vector>// arreglos
#include <algorithm> // heap

#include <cmath>    // Para sqrt y pow
#include <cstdint>  // Para int64_t

#include <queue>
#include <vector>


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


// arista MinHeap;
struct Arista {
    double peso;

    // Para priority_queue, necesitas sobrecargar el operador <
    // Para min-heap (menor peso tiene más prioridad)
    bool operator<(const Arista& otra) const {
        return peso > otra.peso; // Cambia a < para max-heap
    }
};




// ejempl de uso MAIN;
int main() {
    // Vector de aristas PARA USAR SORT
    std::vector<Arista> aristas;
    // push_back añade elemento al final. es como append en python
    // El vector mantiene el orden de inserción 
    aristas.push_back({2.0});
    aristas.push_back({0.5});
    aristas.push_back({1.1});

    // Ordenar el vector de menor a mayor peso CON SORT
    std::sort(aristas.begin(), aristas.end(), [](const Arista& a, const Arista& b){
        return a.peso < b.peso;
    });

    std::cout << "Aristas ordenadassss: ";
    for (auto& a : aristas) std::cout << a.peso << " ";
    std::cout << std::endl;

    // Heap de aristas (menor peso al tope)
    std::priority_queue<Arista> heap;
    //  agrega el elemento al heap y automáticamente lo reubica internamente para mantener la propiedad de heap.
    heap.push({2.0});
    heap.push({0.5});
    heap.push({1.1});

    std::cout << "Aristas desde el heap: ";
    while (!heap.empty()) {
        std::cout << heap.top().peso << " ";
        heap.pop();
    }
    std::cout << std::endl;
    return 0;
}