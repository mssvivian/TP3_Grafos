#include <bits/stdc++.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <queue> 
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include "graph.h" 
using namespace std;

void printFordFulkersonResult(const tuple<int, vector<tuple<int, int, int>>>& result) {
    // Extrai o fluxo máximo e o vetor de arestas do resultado
    int fluxoMaximo = get<0>(result);
    const auto& arestas = get<1>(result);


    // Itera sobre as arestas e imprime os detalhes
    cout << "Arestas no fluxo residual:" << endl;
    for (const auto& aresta : arestas) {
        int u = get<0>(aresta); // Nó de origem
        int v = get<1>(aresta); // Nó de destino
        int fluxo = get<2>(aresta); // fluxo alocado
        cout << "Origem: " << u << ", Destino: " << v << ", Fluxo Alocado: " << fluxo << endl;
    }

    // Imprime o fluxo máximo
    cout << "Fluxo Máximo: " << fluxoMaximo << endl;

}

int main(int argc, char *argv[]) {
    Grafo grafo = Grafo("grafo_rf_2.txt",false,true,true);
    tuple<int, vector<tuple<int, int, int>>> result = grafo.algoritmo_Ford_Fulkerson(1,2);
    cout<< "fluxo máximo: " << get<0>(result) << endl;
    //printFordFulkersonResult(result);

    return 0;
}

// g++ -O3 -o a.exe main.cpp graph.cpp
// ./a.exe --file grafo_0.txt --useMatrix nao
