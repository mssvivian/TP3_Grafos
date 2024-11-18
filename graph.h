#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <queue> 
#include <chrono>
#include <algorithm>
using namespace std;

class Grafo {
private:
    vector<vector<float>> matrizAdj_peso;
    vector<vector<tuple<int,float>>> listAdj_peso;
    vector<vector<bool>> matrizAdj;
    vector<vector<int>> listAdj;
    int vertices;
    int arestas;
    bool grafo_peso;
    bool usaMatriz;
    bool grafo_direcionado;
    vector<bool> visitado;
    vector<int> pai;
    vector<int> nivel;
    vector<float> dist;
    int n_CompConexa; 
    map<string, int> mapa_nomes;
    vector<string> vetor_nomes;

public:

    vector<vector<tuple<int, int, int>>> listAdj_dic;
    // Construtor
    Grafo(const string &grafo, bool isMatriz, bool peso, bool direcionado);

    // Get vértices
    int getVertices();

    // Criar txt
    void write_general_info(const string &arquivo);
    
    // Graus de um vértice
    tuple<int, int, int, float> graus();

    // Imprime a Matriz de Adjacência
    void printMatrizAdj() const;

    // Imprime a Lista de Adjacência
    void printListAdj() const;

    // Escrever em arquivo (para retornar a árvore gerada por busca)
    void Write_file_busca(vector<int> pai, vector<int> nivel, const string& outputFile, int s);

    // Escrever em arquivo (para retornar a árvore gerada por busca) -- grafo com peso
    void Write_file_busca_peso(vector<int> pai, vector<float> dist, const string& outputFile, int s);

    // Obtém os vizinhos de um vértice na matriz de adjacência
    vector<int> getVizinhosMatriz(int v) const;

    // Obtém os vizinhos de um vértice na matriz de adjacência com pesos
    vector<tuple<int,float>> getVizinhosMatriz_peso(int v) const;

    // Algoritmo Diâmetro aproximado
    int aprox(int start, const string& outputFile, bool write_tree);

    // Algoritmo diâmetro não aproximado
    int diameter();

    // Algoritmo Distância 
    int distancia(int start, const string& outputFile, bool write_tree, int end);

    // Algoritmo caminho -- sem peso (arestas) 
    vector<int> caminho_minimo(int start, const string& outputFile, bool write_tree, int end);

    // Algoritmo BFS
    void BFS(int s, const string& outputFile, bool write_tree, int e = 0);

    // Algoritmo DFS
    void DFS(int s, const string& outputFile, bool write_tree);

    // Componentes Conexas
    map<int, vector<vector<int>>, greater<int>> ComponentesConexas();

    // Imprimir Componentes Conexas
    void imprimirComponentesConexas(const map<int, vector<vector<int>>, greater<int>>& componentesConexas);

    // Algoritmo de Dijkstra que calcula a distancia de 1 vertice para todos os outros, retorna o vetor de distancia e a arvore geradora 
    vector<float> algoritmo_Dijkstra(int s, const string& outputFile, bool write_tree, bool Heap);

    // funçao distancia -- entra com booleano de heap, vertice inicial e final, retorna a distancia entre esses 2 pontos
    float distancia_peso(int start, int end, bool Heap);

    // funçao distancia -- entra com booleano de heap, vertice inicial e final, retorna a distancia entre esses 2 pontos
    vector<int> caminho_minimo_peso(int start, int end, bool Heap);

        // Função para ler o arquivo e preencher o mapa e o vetor
    void lerArquivoParaMapEVetor(const string& nomeArquivo);

    // função transforma o caminho mínimo dos índices para o caminho mínimo com o nome de cada pesquisador
    vector <string> caminho_minimo_nomes(string inicio, string fim);

    // função que calcula e retorna a distancia entre dois nomes na rede de colaboração 
    float distancia_nomes(string inicio, string fim, bool Heap);
    
    // Criar um grafo residual com os mesmos vértices e configuração (direcionado, com pesos, etc.)
    Grafo criarGrafoResidual();

    //função retorna uma tupla,, onde a primeira posição é um inteiro que representa o fluxo máximo 
    //e a segunda é um vetor de tuplas de 3 posições, vertice i, vertice j, fluxo da aresta i -> j
    tuple<int, vector<tuple<int, int, int>>> algoritmo_Ford_Fulkerson(int s, 
    int t);

    // função que calcula gargalo de um caminho
    float calcularGargalo(const vector<int>& caminho, const Grafo& grafoResidual);

    // função que dado um caminho e um gargalo, atualiza o fluxo no grafo original e 
    // as arestas (e suas capaciodades) no grafo residual
    void atualizarFluxo(Grafo& grafoResidual, const vector<int>& caminho, int gargalo);

};
#endif // GRAPH_H
