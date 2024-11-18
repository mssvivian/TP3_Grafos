#include <bits/stdc++.h>
using namespace std;
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <queue> 
#include <chrono>
#include <algorithm>
#include "graph.h"


// Definição do construtor
Grafo::Grafo(const string &grafo, bool isMatriz, bool peso, bool direcionado)
    : usaMatriz(isMatriz), n_CompConexa(0), grafo_peso(peso), grafo_direcionado(direcionado) {

    ifstream file(grafo);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo!" << endl;
        exit(EXIT_FAILURE);
    }

    // Lê o número de vértices
    file >> vertices;
    
    visitado.resize(vertices + 1, false);
    pai.resize(vertices + 1, -1);
    nivel.resize(vertices + 1, -1);
    arestas = 0;

    if (grafo_direcionado){
        listAdj_dic.resize(vertices + 1);

        // Lê as arestas e preenche a lista de adjacência - apenas grafo direcionado
        // lendo u,v constroi a aresta u -> v
        int u, v;
        int w;
        while (file >> u >> v >> w) {
            tuple<int, int, int> t1 = make_tuple(v, w, 0);
            listAdj_dic[u].push_back(t1);
            arestas++;
        }

    } else if (usaMatriz) {
        if (grafo_peso){
            matrizAdj_peso.resize(vertices + 1, vector<float>(vertices + 1, numeric_limits<float>::infinity()));

            // Lê as arestas e preenche a matriz de adjacência
            int u, v;
            float w;
            while (file >> u >> v >> w) {
                matrizAdj_peso[u][v] = w;
                matrizAdj_peso[v][u] = w; // Se o grafo for não direcionado
                arestas++;
            }
        } else{
            matrizAdj.resize(vertices + 1, vector<bool>(vertices + 1, false));

            // Lê as arestas e preenche a matriz de adjacência
            int u, v;
            while (file >> u >> v) {
                matrizAdj[u][v] = true;
                matrizAdj[v][u] = true; // Se o grafo for não direcionado
                arestas++;

            }
        }

    } else {
        if (grafo_peso){
            listAdj_peso.resize(vertices + 1);

            // Lê as arestas e preenche a lista de adjacência
            int u, v;
            float w;
            while (file >> u >> v >> w) {
                tuple<int, float> t1 = make_tuple(u, w);
                tuple<int, float> t2 = make_tuple(v, w);
                listAdj_peso[u].push_back(t2);
                listAdj_peso[v].push_back(t1); // Se o grafo for não direcionado
                arestas++;
            }
        } else{
            listAdj.resize(vertices + 1);

            // Lê as arestas e preenche a lista de adjacência
            int u, v;
            while (file >> u >> v) {
                listAdj[u].push_back(v);
                listAdj[v].push_back(u); // Se o grafo for não direcionado
                arestas++;
            }
        }
    }

    file.close();
    }

// Get vertices
int Grafo::getVertices(){
    return vertices;
}


// Graus do grafo (máx,min,mediana e média) - funciona apenas para grafo não direcionado 
tuple<int, int, int, float> Grafo::graus() {
    if (grafo_direcionado){
        cerr << "Essa função é apenas para grafos não direcionados" << endl;
        tuple<int, int, int, float> erro = make_tuple(0,0,0,0.0);
        return erro;
    }
    int max = INT_MIN;  // Valor mínimo possível para garantir a atualização
    int min = INT_MAX;  // Valor máximo possível para garantir a atualização
    int mediana = 0;
    float media = 0.0;
    vector<int> graus(vertices + 1, -1);
    float total = 0.0;
    if (usaMatriz) {
        for (int i = 0; i < vertices; i++) {
            graus[i] = accumulate(matrizAdj[i].begin(), matrizAdj[i].end(), 0); // Calcula graus
            total = total + graus[i];
        }
    } else {
        for (int i = 0; i < vertices; i++) {
            graus[i] = listAdj[i].size(); // Calcula graus na lista de adjacência
            total = total + graus[i];
        }
    }
    sort(graus.begin(),graus.end());
    // Encontrar os graus máximo e mínimo
    max = graus.back();
    min = graus.front();
    int size = graus.size();
    mediana = graus[size/2];
    media = total/vertices;

    tuple<int, int, int, float> resultado = make_tuple(max,min,mediana,media);

    return resultado;
}

// Criar txt
void Grafo::write_general_info(const string &arquivo){
    tuple<int, int, int, float> t = graus();
    ofstream outFile(arquivo);
    if (!outFile.is_open()) {
        cerr << "Erro ao abrir arquivo de saída." << endl;
        return;
    }

    outFile << "Informações gerais do grafo.\n";
    outFile << "Vértices: " << vertices << endl;
    outFile << "Arestas: " << arestas << endl;
    outFile << "Grau máximo: " << get<0>(t) << endl;
    outFile << "Grau mínimo: " << get<1>(t) << endl;
    outFile << "Mediana do grau: " << get<2>(t) << endl;
    outFile << "Média do grau: " << get<3>(t) << endl;
    
    outFile.close();

}

// Imprime a Matriz de Adjacência
void Grafo::printMatrizAdj() const {
    for (size_t i = 1; i < matrizAdj.size(); ++i) {
        for (size_t j = 1; j < matrizAdj[i].size(); ++j) {
            cout << matrizAdj[i][j] << " ";
        }
        cout << endl;
    }
}

// Imprime a Lista de Adjacência
void Grafo::printListAdj() const {
    for (size_t i = 1; i < listAdj.size(); ++i) {
        cout << i << ": ";
        for (size_t j = 0; j < listAdj[i].size(); ++j) {
            cout << listAdj[i][j] << " ";
        }
        cout << endl;
    }
}

//escrever em arquivo (para retornar a arvore gerada por busca)
void Grafo::Write_file_busca(vector<int> pai,vector<int> nivel, const string& outputFile, int s){
    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cerr << "Erro ao abrir arquivo de saída." << endl;
        return;
    }
    outFile << "Árvore gerada por busca\n";

    for (int v = 1; v <= vertices; v++){
    if (pai[v] != -1 || v == s) {
        outFile << "Vértice: " << v << ", Pai: " << pai[v] << ", Nível: " << nivel[v] << endl;
        }
    }
    outFile.close();
}

// Escrever em arquivo (para retornar a árvore gerada por busca) -- grafo com peso
void Grafo::Write_file_busca_peso(vector<int> pai, vector<float> dist, const string& outputFile, int s){
    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cerr << "Erro ao abrir arquivo de saída." << endl;
        return;
    }
    outFile << "Árvore gerada por busca\n";

    for (int v = 1; v <= vertices; v++){
    if (pai[v] != -1 || v == s) {
        outFile << "Vértice: " << v << ", Pai: " << pai[v] << ", Distância: " << dist[v] << endl;
        }
    }
    outFile.close();
}

// Obtém os vizinhos de um vértice na matriz de adjacência sem pesos
vector<int> Grafo::getVizinhosMatriz(int v) const {
    vector<int> vizinhos;
    for (int i = 1; i <= vertices; ++i) {
        if (matrizAdj[v][i] == 1) {
            vizinhos.push_back(i);
        }
    }
    return vizinhos;
}

// Obtém os vizinhos de um vértice na matriz de adjacência com pesos
vector<tuple<int,float>> Grafo::getVizinhosMatriz_peso(int v) const {
    vector<tuple<int,float>> vizinhos;
    for (int i = 1; i <= vertices; ++i) {
        if (matrizAdj_peso[v][i] != numeric_limits<float>::infinity()) {
            tuple<int, float> t1 = make_tuple(i, matrizAdj_peso[v][i]);
            vizinhos.push_back(t1);
        }
    }
    return vizinhos;
    
}

// Algoritmo Diâmetro aproximado -- acho que dá pra tirar ou teriamos que redefinir para pesos 
int Grafo::aprox(int start, const string& outputFile, bool write_tree){
    if (grafo_peso){
        cerr << "Essa função é apenas para grafos sem peso" << endl;
        return 0;
    }
    BFS (start, outputFile, write_tree); // encontrando a primeira extremidade (x)
    int maior = 0;
    int x;
    for (int i = 1; i < nivel.size(); i++){
        if (nivel[i]> maior){
            maior = nivel[i];
            x = i;
        }
    }
    BFS (x, outputFile, write_tree); // encontrando a segunda extremidade começando da primeira extremidade (x)
    for (int i = 1; i < nivel.size(); i++){
        if (nivel[i]> maior){
            maior = nivel[i];
            x = i;
        }
    }
    return maior;
}

// Algoritmo diâmetro não aproximado -- tbm considera sem peso
int Grafo::diameter() {
    if (grafo_peso){
        cerr << "Essa função é apenas para grafos sem peso" << endl;
        return 0;
    }
    int x = 0;
    int maior = 0;

    for (int s = 1; s < vertices; s++) {
        visitado.assign(vertices + 1, false);
        nivel.assign(vertices + 1, -1);

        queue<int> Q;
        
        visitado[s] = true;
        nivel[s] = 0;
        Q.push(s);

        while (!Q.empty()) {

            int v = Q.front();
            Q.pop();

            const vector<int>& vizinhos = usaMatriz ? getVizinhosMatriz(v) : listAdj[v];
            for (int w : vizinhos) {
                if (!visitado[w]) {
                    visitado[w] = true;
                    nivel[w] = nivel[v] + 1; // Define o nível de w
                    Q.push(w);
                    if (nivel[w]>maior){
                        maior = nivel[w];
                        x = w;
                    }
                }
            }
        }
    }

    return maior;
}

// Algoritmo Distância -- sem peso 
int Grafo::distancia(int start, const string& outputFile, bool write_tree, int end){
    int level;
    BFS(start, outputFile, write_tree,end);
    level = nivel[end];
    return level;
}

// Algoritmo caminho -- sem peso (arestas) 
vector<int> Grafo::caminho_minimo(int start, const string& outputFile, bool write_tree, int end){
    BFS(start, outputFile, write_tree,end);


    if (!visitado[end]) {
    return {}; // Retorna vazio se o nó final não foi alcançado
    }
    
    // Vetor para armazenar o caminho mínimo
    vector<int> caminho;

    // Reconstrói o caminho mínimo a partir do vetor de pais
    for (int v = end; v != -1; v = pai[v]) {
        caminho.push_back(v);
    }

    // O caminho está ao contrário, então invertemos
    reverse(caminho.begin(), caminho.end());

    // Verifica se o caminho encontrado é válido (se o último vértice é 'start')
    if (caminho.front() != start) {
        cout << "caminho não é válido" << endl;
        return {}; // Retorna um vetor vazio se não houver caminho
    }

    return caminho;
}

// Algoritmo BFS -- sem peso 
void Grafo::BFS(int s, const string& outputFile, bool write_tree, int e) {

    visitado.assign(vertices + 1, false);
    pai.assign(vertices + 1, -1);
    nivel.assign(vertices + 1, -1);

    queue<int> Q;
    
    visitado[s] = true;
    nivel[s] = 0;
    Q.push(s);
    vector<int> vizinhos;

    while (!Q.empty()) {
        if (e!= 0 and visitado[e]){
            break;
        }
        int v = Q.front();
        Q.pop();
        if (grafo_direcionado){
            for (const auto& tupla : listAdj_dic[v]) {
                //vizinhos.push_back(get<0>(tupla)); // Obtém o primeiro elemento da tupla
                int w = get<0>(tupla);
                int capacidade = get<1>(tupla);
                if (!visitado[w] and capacidade > 0){
                    visitado[w] = true;
                    pai[w] = v;      // Armazena o pai de w
                    nivel[w] = nivel[v] + 1; // Define o nível de w
                    Q.push(w);
                }
            }
        } else {
            vector<int> vizinhos = usaMatriz ? getVizinhosMatriz(v) : listAdj[v];
            for (int w : vizinhos) {
                if (!visitado[w]) {
                    visitado[w] = true;
                    pai[w] = v;      // Armazena o pai de w
                    nivel[w] = nivel[v] + 1; // Define o nível de w
                    Q.push(w);
                }
            }
        }
    }
    if (write_tree){
        Write_file_busca(pai, nivel, outputFile, s);
    }
}

// Algoritmo DFS -- sem peso 
void Grafo::DFS(int s, const string& outputFile, bool write_tree) {

    visitado.assign(vertices + 1, false);
    pai.assign(vertices + 1, -1);
    nivel.assign(vertices + 1, -1);

    stack<int> pilha;

    pilha.push(s);
    nivel[s] = 0;

    while (!pilha.empty()) {
        int v = pilha.top();
        pilha.pop();

        if (!visitado[v]) {
            visitado[v] = true;

            // Para cada vizinho do vértice v
            const vector<int>& vizinhos = usaMatriz ? getVizinhosMatriz(v) : listAdj[v];
            for (int w : vizinhos) {
                if (!visitado[w]) {
                    pilha.push(w);
                    pai[w] = v;           // Define o pai de w
                    nivel[w] = nivel[v] + 1; // Define o nível de w
                }
            }
        }
    }
    if (write_tree){
        Write_file_busca(pai, nivel, outputFile, s);
    }
}

// Algoritmo Componentes Conexas
map<int, vector<vector<int>>, greater<int>> Grafo::ComponentesConexas() {   
    // Mapa onde a chave é o tamanho da componente e os valores são os vértices da componente
    map<int, vector<vector<int>>, greater<int>> componentesConexas;
    vector<bool> visitado(vertices + 1, false); // Vetor para marcar vértices visitados
    queue<int> Q; // Fila para a BFS

    // Percorre todos os vértices do grafo
    for (int i = 1; i <= vertices; ++i) {
        if (!visitado[i]) {
            // Se o vértice ainda não foi visitado, achamos uma nova componente conexa
            vector<int> componenteAtual; // Armazena os vértices da componente atual
            n_CompConexa ++;
            // Executa a BFS a partir do vértice 'i'
            Q.push(i);
            visitado[i] = true;
            vector<int> vizinhos;

            while (!Q.empty()) {
                int v = Q.front();
                Q.pop();

                componenteAtual.push_back(v); // Adiciona o vértice à componente atual

                if (grafo_direcionado){
                    vector<int> vizinhos;
                    for (const auto& tupla : listAdj_dic[v]) {
                        vizinhos.push_back(get<0>(tupla)); // Obtém o primeiro elemento da tupla
                    }
                } else {
                    vector<int> vizinhos = usaMatriz ? getVizinhosMatriz(v) : listAdj[v];
                }
                for (int w : vizinhos) {
                    if (!visitado[w]) {
                        visitado[w] = true;
                        Q.push(w);
                    }
                }
            }

            // Adiciona a componente atual ao mapa, usando o tamanho da componente como chave
            int tamanho = componenteAtual.size();
            componentesConexas[tamanho].push_back(componenteAtual); 
        }
    }

    return componentesConexas; // Retorna o mapa de componentes
}

// Imprimir Componentes Conexas
void Grafo::imprimirComponentesConexas(const map<int, vector<vector<int>>, greater<int>>& componentesConexas) {
    if (componentesConexas.empty()) {
        cout << "Nenhuma componente conexa encontrada." << endl;
        return;
    }

    // Itera sobre o map onde a chave é o tamanho da componente
    for (const auto& par : componentesConexas) {
        long long tamanho = par.first; // A chave (tamanho da componente)
        const vector<vector<int>>& componentes = par.second; // O vetor de componentes conexas de mesmo tamanho

        cout << "Componentes de tamanho " << tamanho << ":" << endl;

        // Itera sobre cada componente conexa (um vetor de vértices)
        for (const auto& componente : componentes) {
            cout << "Vértices: ";
            for (int vertice : componente) {
                cout << vertice << " "; // Imprime cada vértice da componente
            }
            cout << endl;
        }

        cout << endl; // Nova linha entre as componentes de mesmo tamanho
    }
}


// Algoritmo de Dijkstra que calcula a distância de um vértice para todos os outros
// Retorna o vetor de distâncias e opcionalmente escreve a árvore geradora em arquivo
vector<float> Grafo::algoritmo_Dijkstra(int s, const string& outputFile, bool write_tree, bool Heap) {
    // Inicializando o vetor de distâncias com infinito
    dist.assign(vertices + 1, numeric_limits<float>::infinity());
    pai.assign(vertices + 1, -1);
    vector<bool> visitado(vertices + 1, false);
    if (!grafo_peso){
        cerr << "Essa função é apenas para grafos com peso" << endl;
        return dist;
    }
    
    // A distância do vértice inicial para ele mesmo é 0
    dist[s] = 0;

    if (Heap) {
        // Usando priority_queue (min-heap) para o caso Heap == true
        // A fila armazena pares (distância, vértice)
        priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> pq;
        pq.push({0, s});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();
            if (!visitado[u]){
                // Para cada vizinho v de u
                visitado[u] = true;
                const vector<tuple<int,float>>& vizinhos = usaMatriz ? getVizinhosMatriz_peso(u) : listAdj_peso[u];
                for (const auto& vizinho : vizinhos) {
                    int v = get<0>(vizinho);
                    float peso = get<1>(vizinho);
                    if (peso < 0){
                        cerr << "A biblioteca ainda nao implementa caminhos minimos com pesos negativos." << endl;
                        return dist;
                    }
                    if (dist[v] > dist[u] + peso) {
                        dist[v] = dist[u] + peso;  // Atualiza a distância de v
                        pai[v] = u;  // Armazena o pai para reconstruir o caminho
                        pq.push({dist[v], v});  // Adiciona v à fila de prioridade
                    }
                }
            }
        }
    } else {
        // Usando vetor para o caso Heap == false
        
        for (int i = 0; i < vertices; ++i) {
            // Seleciona o vértice u não visitado com a menor distância
            float menorDist = numeric_limits<float>::infinity();
            int u;

            for (int v = 1; v <= vertices; ++v) {
                if (!visitado[v] && dist[v] < menorDist) {
                    menorDist = dist[v];
                    u = v;
                }
            }

            // Marca u como visitado
            visitado[u] = true;

            // Para cada vizinho v de u
            const vector<tuple<int,float>>& vizinhos = usaMatriz ? getVizinhosMatriz_peso(u) : listAdj_peso[u];

            for (const auto& vizinho : vizinhos) {
                int v = get<0>(vizinho);
                float peso = get<1>(vizinho);
                if (peso < 0){
                    cerr << "A biblioteca ainda nao implementa caminhos minimos com pesos negativos." << endl;
                    return dist;
                }
                if (dist[v] > dist[u] + peso) {
                    dist[v] = dist[u] + peso;  // Atualiza a distância de v
                    pai[v] = u;  // Armazena o pai para reconstruir o caminho
                }
            }
        }
    }

    // Se write_tree for true, grava a árvore geradora no arquivo de saída
    if (write_tree) {
        Write_file_busca_peso(pai, dist, outputFile, s);
    }
    return dist;
}

// funçao distancia -- entra com booleano de heap, vertice inicial e final, retorna a distancia entre esses 2 pontos
float Grafo::distancia_peso(int start, int end, bool Heap) {
    // Calcula a distância de 'start' para todos os outros vértices
    vector<float> distancias = algoritmo_Dijkstra(start, "", false, Heap);

    // Retorna a distância do vértice 'start' para o vértice 'end'
    return distancias[end];
}
// funçao caminho minimos -- entra com booleano de heap, vertice inicial e final, retorna o caminho mínimo entre esses 2 pontos
vector<int> Grafo::caminho_minimo_peso(int start, int end, bool Heap) {
    // Calcula a distância de 'start' para todos os outros vértices e constrói a árvore geradora
    algoritmo_Dijkstra(start, "", false, Heap);

    // Vetor para armazenar o caminho mínimo
    vector<int> caminho;

    // Reconstrói o caminho mínimo a partir do vetor de pais
    for (int v = end; v != -1; v = pai[v]) {
        caminho.push_back(v);
    }

    // O caminho está ao contrário, então invertemos
    reverse(caminho.begin(), caminho.end());

    // Verifica se o caminho encontrado é válido (se o último vértice é 'start')
    if (caminho.front() != start) {
        return {}; // Retorna um vetor vazio se não houver caminho
    }

    return caminho;
}

// Função para ler o arquivo e preencher o mapa e o vetor
void Grafo::lerArquivoParaMapEVetor(const string& nomeArquivo) {

    ifstream arquivo(nomeArquivo);
    string linha;

    while (getline(arquivo, linha)) {
        stringstream ss(linha);
        string indice_str, nome;
        
        while (getline(ss, indice_str, ',') && getline(ss, nome)) {
            int indice = stoi(indice_str);
            mapa_nomes[nome] = indice;
            if (indice >= vetor_nomes.size()) {
                vetor_nomes.resize(indice + 1);
            }
            vetor_nomes[indice] = nome;
        }
    }

}

// função transforma o caminho mínimo dos índices para o caminho mínimo com o nome de cada pesquisador
vector <string> Grafo::caminho_minimo_nomes(string inicio, string fim){
    
    int id_inicio = mapa_nomes[inicio];
    int id_fim = mapa_nomes[fim];
    vector<int> resultado_busca = caminho_minimo_peso(id_inicio,id_fim,true);
    vector<string> output;

    for (int i = 0; i < resultado_busca.size(); ++i){
        output.push_back(vetor_nomes[resultado_busca[i]]); // adiciona o nome do pesquisador ao novo caminho

    }

    return output;
}

// função que calcula e retorna a distancia entre dois nomes na rede de colaboração 
float Grafo::distancia_nomes(string inicio, string fim, bool Heap) {
    int id_inicio = mapa_nomes[inicio];
    int id_fim = mapa_nomes[fim];
    // Calcula a distância de 'start' para todos os outros vértices
    vector<float> distancias = algoritmo_Dijkstra(id_inicio, "", false, Heap);

    // Retorna a distância do vértice 'start' para o vértice 'end'
    return distancias[id_fim];
}

Grafo Grafo::criarGrafoResidual() {
    // Criar um grafo residual com os mesmos vértices e configuração (direcionado, com pesos, etc.)
    Grafo grafoResidual("a.txt", false, true, true); // Parâmetros: arquivo vazio, lista de adjacência, com pesos, direcionado

    // Inicializar o grafo residual com base nas arestas do grafo original
    grafoResidual.vertices = vertices;
    grafoResidual.listAdj_dic.resize(vertices + 1);

    for (int u = 1; u < listAdj_dic.size(); u++) {
        if (listAdj_dic[u].size() != 0){
            for (size_t i = 0; i < listAdj_dic[u].size() ; i++) {
                tuple<int, int, int> aresta = listAdj_dic[u][i];
                int vizinho = get<0>(aresta);
                float capacidade = get<1>(aresta);
                float fluxo = get<2>(aresta);
                // Adicione aresta original no grafo residual
                grafoResidual.listAdj_dic[u].emplace_back(vizinho, capacidade, 1);
                // Adicione aresta reversa no grafo residual
                grafoResidual.listAdj_dic[vizinho].emplace_back(u, 0, 0);
                
            } 
        } 
     
    } 

    return grafoResidual;
}


//função retorna uma tupla, onde a primeira posição é um inteiro que representa o fluxo máximo 
//e a segunda é um vetor de tuplas de 3 posições, vertice i, vertice j, fluxo da aresta i -> j
tuple<int, vector<tuple<int, int, int>>> Grafo::algoritmo_Ford_Fulkerson(int s, 
int t){
    int fluxo = 0;
    vector<tuple<int, int, int>> arestas_e_fluxo;
    Grafo residual = criarGrafoResidual();

    while (true) {
        vector<int> caminho_atual = residual.caminho_minimo(s,"a.txt",false,t);
        
        if (caminho_atual.empty()){
            break;
        } 
        float gargalo_atual = calcularGargalo(caminho_atual, residual);
        
        atualizarFluxo(residual,caminho_atual,gargalo_atual);
        
    }
    
    for (size_t i = 0; i < listAdj_dic[s].size(); i++) {
        tuple<int, int, int> aresta = listAdj_dic[s][i];
        int vizinho_origin = get<0>(aresta);
        int fluxo_origin = get<2>(aresta);
        fluxo += fluxo_origin;
    }
    for (size_t i = 1; i < vertices; i++) {
        for (size_t j = 0; j < listAdj_dic[i].size(); j++){
        tuple<int, int, int> aresta = listAdj_dic[i][j];
        int vizinho_origin = get<0>(aresta);
        int fluxo_origin = get<2>(aresta);
        arestas_e_fluxo.push_back(make_tuple(i,vizinho_origin,fluxo_origin));
        }
    }
    tuple<int, vector<tuple<int, int, int>>> result = make_tuple(fluxo,arestas_e_fluxo);
    return result;
}

// função que calcula gargalo de um caminho
float Grafo::calcularGargalo(const vector<int>& caminho, const Grafo& grafoResidual) {
    float gargalo = numeric_limits<float>::infinity(); // Comece assumindo capacidade infinita

    // Percorra o caminho para encontrar o menor valor de capacidade residual
    for (size_t i = 0; i < caminho.size() - 1; ++i) {
        int u = caminho[i];     // Vértice atual
        int v = caminho[i + 1]; // Próximo vértice no caminho

        bool arestaEncontrada = false;

        // Encontre a capacidade residual da aresta (u, v) no grafo residual
        for (size_t i = 0; i < grafoResidual.listAdj_dic[u].size() ; i++) {
            tuple<int, int, int> aresta = grafoResidual.listAdj_dic[u][i];
            int vizinho = get<0>(aresta);
            float capacidade = get<1>(aresta);
            float fluxo = get<2>(aresta);
            if (vizinho == v) { // Se a aresta (u, v) for encontrada
                gargalo = min(gargalo, capacidade); // Atualize o gargalo com o menor valor encontrado
                arestaEncontrada = true;
                break;
            }
        }

        // Verifique se a aresta foi encontrada
        if (!arestaEncontrada) {
            cerr << "Erro: Aresta " << u << " -> " << v << " não encontrada no grafo residual!" << endl;
            return 0.0; // Retorne 0 para indicar que o caminho não é válido
        }
    }

    return gargalo; // Retorne o menor valor encontrado
}


// função que dado um caminho e um gargalo, atualiza o fluxo no grafo original e 
// as arestas (e suas capaciodades) no grafo residual
void Grafo::atualizarFluxo(Grafo& grafoResidual, const vector<int>& caminho, int gargalo) {
    // Percorre os pares de vértices no caminho
    for (size_t i = 0; i < caminho.size() - 1; ++i) {
        int u = caminho[i];
        int v = caminho[i + 1];

        // Atualiza a capacidade no grafo residual no sentido original (u -> v)

        for (size_t i = 0; i < grafoResidual.listAdj_dic[u].size() ; i++){

            tuple<int,int,int> aresta = grafoResidual.listAdj_dic[u][i];
            int vizinho = get<0>(aresta);
            int capacidade = get<1>(aresta);
            int fluxo = get<2>(aresta);

            if (vizinho == v) {
                grafoResidual.listAdj_dic[u][i] = make_tuple(vizinho, capacidade - gargalo, fluxo); // Atualiza a capacidade
            
                break;
            }
        }

        // Adiciona ou atualiza a aresta no sentido reverso (v -> u) no grafo residual
        for (size_t i = 0; i < grafoResidual.listAdj_dic[v].size() ; i++){

            tuple<int,int,int> aresta = grafoResidual.listAdj_dic[v][i];
            int vizinho = get<0>(aresta);
            int capacidade = get<1>(aresta);
            int fluxo = get<2>(aresta);

            if (vizinho == u) {
                grafoResidual.listAdj_dic[v][i] = make_tuple(vizinho, capacidade + gargalo, fluxo); // Aumenta a capacidade
                break;
            }
        }

        // Atualiza o fluxo no grafo original (se existir)
        
        for (size_t i = 0; i < listAdj_dic[u].size() ; i++){

            tuple<int,int,int> aresta = listAdj_dic[u][i];
            int vizinho = get<0>(aresta);
            float capacidade = get<1>(aresta);
            float fluxo = get<2>(aresta);

            if (vizinho == v) {
                listAdj_dic[u][i] = make_tuple(vizinho, capacidade, fluxo + gargalo); // Aumenta o fluxo no grafo original
                break;
            }
        }

        // Atualiza o fluxo reverso no grafo original (se existir)
        for (size_t i = 0; i < listAdj_dic[v].size() ; i++){

            tuple<int,int,int> aresta = listAdj_dic[v][i];
            int vizinho = get<0>(aresta);
            float capacidade = get<1>(aresta);
            float fluxo = get<2>(aresta);

            if (vizinho == u) {
                listAdj_dic[v][i] = make_tuple(vizinho, capacidade, fluxo - gargalo); // Reduz o fluxo reverso
                break;
            }
        }
    }
}



