#include <bits/stdc++.h>

typedef long long int lli;

std::unordered_map<int, std::pair<float,float>> read_tsp(const std::string &filename){
    std::unordered_map<int, std::pair<float,float>> nodeToCoord;          // Vector para las distancias
    std::ifstream file(filename);                // Para leer el archivo   
    std::string line, keyword;                   // Strings para guardar lineas del archivo y segmentos de esas lineas
    float size;

    while(getline(file, line)){             // Para recorrer el archivo
        std::stringstream partLine(line);        // Para convierte el string line en una fuente de entrada
        partLine >> keyword;                // Leo un string de partline

        if(keyword == "DIMENSION") {
            char c;
            partLine >> c >> size;
        }
        else if(keyword == "NODE_COORD_SECTION"){       // Empieza seccion de la matriz
            int node;
            std::pair<float,float> coordinates;

            for(int i = 0; i < size; i++){
                getline(file, line);                
                partLine = std::stringstream(line);

                partLine >> node >> coordinates.first >> coordinates.second;
                // std::cout << node-1 << std::endl;
                nodeToCoord[node-1] = coordinates;
            }
        }
    }
    file.close();
    return nodeToCoord;
}

lli ceil_2D(const std::pair<float,float> &x, const std::pair<float,float> &y){
    float first = x.first - y.first, second = x.second - y.second;

    return int(ceil(sqrt(pow(first, 2) + pow(second, 2))));
    // return int(sqrt(pow(first, 2) + pow(second, 2)));
}

void farthenestInsertion(std::vector<int> &answer, std::unordered_map<int, std::pair<float,float>> &nodeToCoord, lli &cost){
    int n = nodeToCoord.size();
    std::unordered_set<int> remainPts;
    std::vector<int> tour;
    cost = 0;

    for(int i = 0; i < n; i++) remainPts.insert(i);

    // Elejimos el nodo inicial
    std::random_device rd; std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrib(0, n-1);
    int iniNode = distrib(gen);
    tour.push_back(iniNode);
    remainPts.erase(iniNode);

    // Tomamos el segundo nodo (el más lejano al inicial)
    lli maxDstnc = -1;
    int maxNode;
    for(auto pnt:remainPts){
        lli actDist = ceil_2D(nodeToCoord[tour[0]], nodeToCoord[pnt]);
        if (actDist > maxDstnc){
            maxNode = pnt;
            maxDstnc = actDist;
        }
    }
    tour.push_back(maxNode);
    remainPts.erase(maxNode);
    cost += 2*ceil_2D(nodeToCoord[tour[0]], nodeToCoord[tour[1]]);

    //Calculamos las distancias de cada punto en remainPnts al tour
    //Y guardamos el nodo más lejano al tour
    std::unordered_map<int, lli> tourDstnc; 
    int farthestNode;
    maxDstnc = -1;
    for(auto v:remainPts){
        tourDstnc[v] = 0;

        for(int node:tour) tourDstnc[v] += ceil_2D(nodeToCoord[node], nodeToCoord[v]);
        
        if(tourDstnc[v] > maxDstnc){
            farthestNode = v;
            maxDstnc = tourDstnc[v];
        }
    }

    //Ahora si empezamos con el farthernesInsertion
    while(!remainPts.empty()){
        lli penalization = __LONG_LONG_MAX__;
        int whereInsert;                             //Esto nos dice el mejor lugar para insertar el nodo

        for(int i = 0; i < tour.size(); i++){
            //Tomamos las coordenadas de cada nodo
            std::pair<float,float> actual = nodeToCoord[tour[i]], next = nodeToCoord[tour[(i+1)%tour.size()]];
            std::pair<float,float> coordFar = nodeToCoord[farthestNode];
            //Calculamos la penalizacion de agregar el nodo
            lli pivot = ceil_2D(coordFar, actual) + ceil_2D(coordFar, next) - ceil_2D(next, actual);
            
            //Vemos donde se minimiza la penalizacion
            if(pivot < penalization){
                penalization = pivot;
                whereInsert = i;
            }
        }
        
        //Insertamos el nodo en el tour y lo eliminamos del map de distancias
        tour.insert(tour.begin() + whereInsert + 1, farthestNode);
        remainPts.erase(farthestNode);
        tourDstnc.erase(farthestNode);
        cost += penalization;
        int piv = farthestNode;
        
        //Buscamos el nuevo nodo mas lejano
        maxDstnc = -1;
        for(auto v:remainPts){
            tourDstnc[v] += ceil_2D(nodeToCoord[piv], nodeToCoord[v]);
            
            if(tourDstnc[v] > maxDstnc){
                farthestNode = v;
                maxDstnc = tourDstnc[v];
            }
        }
    }
    
    answer = tour;
    return;
}

void randomShuffle(std::vector<int> &v, const int &seed){
    std::mt19937 gen(seed); 
    shuffle(v.begin(), v.end(), gen);
    return;
}

long long int coefBinom(int n) {
    if (n < 2) return 0;  // Si n < 2, no es posible calcular C(n, 2)
    return (n * (n - 1)) / 2;
}

std::pair<int,int> neighbourGen(const int &n, std::set<std::pair<int,int>> &usedNeigh){
    std::random_device rd; std::mt19937 gen(rd());    
    std::uniform_int_distribution<int> dis(0, n-2); 
    int i = dis(gen);
    std::uniform_int_distribution<int> dis2(i+1, n-1); 
    int j = dis2(gen);
    std::pair<int,int> nuevo = std::make_pair(i,j);

    while(usedNeigh.count(nuevo)){
        i = dis(gen);
        std::uniform_int_distribution<> dis2(i+1, n-1); 
        j = dis2(gen);
        nuevo = std::make_pair(i,j);
    }

    usedNeigh.insert(nuevo);

    return nuevo;
}

void esthocastic_hillClimbing(std::vector<int> &actualAns, lli &cost, std::unordered_map<int, std::pair<float,float>> &nodeToCoord){
    lli posNewCost = cost;
    int n = actualAns.size();
    std::pair<int,int> change;

    while(1){
        std::set<std::pair<int,int>> nbhd;

        while(!(nbhd.size() == coefBinom(n))){
            std::pair<int,int> neighbour = neighbourGen(n, nbhd);
            int i = neighbour.first, j = neighbour.second;

            std::pair<float,float> xi, yi, xj, yj;
            if(i == 0 && j == n-1){
                xi = nodeToCoord[actualAns[i]], yi = nodeToCoord[actualAns[i+1]];
                xj = nodeToCoord[actualAns[j]], yj = nodeToCoord[actualAns[j-1]];
            } else {
                xi = nodeToCoord[actualAns[i]], yi = nodeToCoord[actualAns[(i-1+n)%n]];
                xj = nodeToCoord[actualAns[j]], yj = nodeToCoord[actualAns[(j+1)%n]];
            }

            posNewCost -= ceil_2D(xi, yi) + ceil_2D(xj, yj);
            posNewCost += ceil_2D(xj, yi) + ceil_2D(xi, yj);

            if(posNewCost < cost){
                change = neighbour;
                break;
            }
            posNewCost = cost;
        }

        if(posNewCost != cost){
            if(change.first == 0 && change.second == n-1) 
                std::swap(actualAns[0],actualAns[n-1]);
            else{
                int intervalSize = abs(change.second - change.first) + 1;
                for(int i = 0; i < intervalSize/2; i++){
                    std::swap(actualAns[change.first + i], actualAns[change.second - i]);
                }
            }
            cost = posNewCost;
        } else break;
    }

    return;
}

int main(){
    std::vector<std::string> files = {"myInst2"};   //Nombre del archivo

    for (auto filename : files){
        std::string line;
        std::unordered_map<int,std::pair<float,float>> nodeToCoord = read_tsp(filename + ".tsp");

        if (nodeToCoord.empty()){
            std::cout << "No se pudieron cargar las coordenadas" << std::endl;
            return 1;
        }

        lli cost;  
        std::vector<int> sol;
        std::string newFile = filename + "s.txt";              //Nombre del archivo de costos para la instancia filename

        FILE *archivo = fopen(newFile.c_str(), "w");        
        fprintf(archivo, "Costos con la primer solucion dada por la heuristica para la instancia: %s\n\n", filename.c_str());

        for(int i = 0; i < 50; i++){
            double t_0 = std::clock(), t_1, time;

            farthenestInsertion(sol, nodeToCoord, cost);
            esthocastic_hillClimbing(sol, cost, nodeToCoord);
            t_1 = std::clock();
            time = (t_1 - t_0)/CLOCKS_PER_SEC;

            fprintf(archivo, "Solución #%d\n", i+1);
            fprintf(archivo,"Costo: %lld\n", cost);
            fprintf(archivo, "Tiempo de ejecucion: %.5f\n\n", time);
        }

        fclose(archivo);
    }

    return 0;
} 