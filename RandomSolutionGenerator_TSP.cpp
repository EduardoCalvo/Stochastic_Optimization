#include <bits/stdc++.h>
#include <iostream>
#include <stdio.h>
#include <omp.h>

using namespace std;

unordered_map<int, pair<float,float>> read_tsp(const string &filename){
    unordered_map<int, pair<float,float>> nodeToCoord;          // Vector para las distancias
    ifstream file(filename);                // Para leer el archivo   
    string line, keyword;                   // Strings para guardar lineas del archivo y segmentos de esas lineas
    float size;

    while(getline(file, line)){             // Para recorrer el archivo
        stringstream partLine(line);        // Para convierte el string line en una fuente de entrada
        partLine >> keyword;                // Leo un string de partline

        if(keyword == "DIMENSION") {
            char c;
            partLine >> c >> size;
        }
        else if(keyword == "NODE_COORD_SECTION"){       // Empieza seccion de la matriz
            int node;
            pair<float,float> coordinates;

            for(int i = 0; i < size; i++){
                getline(file, line);                
                partLine = stringstream(line);

                partLine >> node >> coordinates.first >> coordinates.second;
                nodeToCoord[node] = coordinates;
            }
        }
    }
    file.close();
    return nodeToCoord;
}

void randomShuffle(vector<int> &v){
    random_device rd; mt19937 gen(rd());
    shuffle(v.begin(), v.end(), gen);
    return;
}

int ceil_2D(const pair<float,float> &x, const pair<float,float> &y){
    float first = x.first - y.first, second = x.second - y.second;

    return int(ceil(sqrt(pow(first, 2) + pow(second, 2))));
}

vector<int> randSolution(unordered_map<int, pair<float,float>> &nodeToCoord, long long int &cost){
    vector<int> posCities;  
    cost = 0; 

    for(int i = 1; i <= nodeToCoord.size(); i++) posCities.push_back(i);
    randomShuffle(posCities);

    for(int i = 0; i < posCities.size(); i++) {
        int node1 = posCities[i], node2 = posCities[(i+1)%posCities.size()];
        pair<float,float> x = nodeToCoord[node1], y = nodeToCoord[node2];
        cost += ceil_2D(x, y);
    }
    return posCities;
}

int main(){
    vector<string> files = {"pla85900"};   //Nombre del archivo "pla7397", 

    for (auto filename : files){
        string line;
        unordered_map<int,pair<float,float>> nodeToCoord = read_tsp(filename + ".tsp");

        if (nodeToCoord.empty()){
            cout << "No se pudieron cargar las coordenadas" << endl;
            return 1;
        }

        vector<int> bestCostes;                        //Para guardar los mejores costos
        int i, j, ref = 100000;
        long long int smallCost = LLONG_MAX;        //Variables para correr en paralelo
        string newFile = filename + ".txt";              //Nombre del archivo de costos para la instancia filename

        FILE *archivo = fopen(newFile.c_str(), "w");        
        fprintf(archivo, "Mejores costos de la generacion aleatoria para la instancia %s\n\n", filename.c_str());
        double t_0, t_1, time;
        t_0 = clock();

    omp_set_num_threads(32);
    #pragma omp parallel for firstprivate(ref, smallCost) private(i)
        for(i = 0; i < 100; i++){                 //Buscamos 100 veces la mejor de 100000 soluciones
            int tid = omp_get_thread_num();
            printf("Este es el hilo %d\n", tid);
            vector<int> bestSol;                      //Vemos si lo ocupamos
            for (j = 0; j < ref; j++){
                long long int pivCost;
                vector<int> sol = randSolution(nodeToCoord, pivCost);
                if (pivCost < smallCost){
                    smallCost = pivCost;
                    //bestSol = sol;                //Vemos si lo ocupamos 
                }
            }
            if(archivo != NULL) fprintf(archivo,"%lld\n", smallCost);
            bestCostes.push_back(smallCost);
        }
        fclose(archivo);
        t_1 = clock();
        time = (t_1 - t_0)/CLOCKS_PER_SEC;
        printf("\nTiempo de busqueda de la mejor entre 100000 soluciones: %lf secods\n", time);
    }
    return 0;
}  
// g++ Tarea1OptiEsto.cpp -o lol -fopenmp
