#include <bits/stdc++.h> 
#include <omp.h>

typedef long long int lli;
typedef std::vector<std::vector<lli>> lliMat;
typedef std::vector<std::vector<int>> mat;

std::pair<int,int> searchLocalImprove(const lliMat &dstMatrix, const mat &nearNodes, const std::vector<int> &actAns, lli &cost, std::unordered_map<int,int> &where, int &adyacent){
    int N = actAns.size();
    lli bestLocalcost = cost;
    std::pair<int,int> change = std::make_pair(-1,-1);

    for(int i = 0; i < N; i++){
        //valores iniciales sobre los que corro
        int node = actAns[i];
        std::pair<int,int> bestAns = std::make_pair(-1,0);
        lli minCost = bestLocalcost;
        //Calculo el mas lejos entre el nodo anterior y el nodo siguiente
        int dstLeft = dstMatrix[node][actAns[(i-1+N)%N]], dstRight = dstMatrix[node][actAns[(i+1)%N]];
        int fartestNb = (dstLeft < dstRight) ? actAns[(i+1)%N] : actAns[(i-1+N)%N];
        //Empiezo a recorrer en la lista de nodos ordenados por mas cercanos a mas lejanos
        int nearIdx = 1;
        do {
            lli cost1 = cost, cost2 = cost;
            //Obtengo el vector cercano y su posicion en la solucion
            int nb = nearNodes[node][nearIdx], jdx = where[nb]; 

            if (jdx == (i+1)%N || jdx == (i-1+N)%N){      //Si estoy checando alguno de sus adyacentes
                int nextBack = (jdx == (i+1)%N) ? 1 : -1; //Para saber si el adyacente es el de adelante o el de atrás
                cost1 -= dstMatrix[node][actAns[(i-nextBack+N)%N]] + dstMatrix[nb][actAns[(jdx+N+nextBack)%N]];
                cost1 += dstMatrix[node][actAns[(jdx+N+nextBack)%N]] + dstMatrix[nb][actAns[(i-nextBack+N)%N]];
                if(cost1 < minCost){
                    minCost = cost1;
                    bestAns = std::make_pair(jdx, nextBack); 
                }
                continue;
            } 
            //Calculo el costo de la primer solucion generada al insertar la arista de node y nb
            cost1 -= dstMatrix[node][actAns[(i+1)%N]] + dstMatrix[nb][actAns[(jdx+1)%N]];
            cost1 += dstMatrix[node][nb] + dstMatrix[actAns[(i+1)%N]][actAns[(jdx+1)%N]];
            //Calculo el costo de la segunda solucion generada al insertar la arista de node y nb
            cost2 -= dstMatrix[node][actAns[(i+N-1)%N]] + dstMatrix[nb][actAns[(jdx+N-1)%N]];
            cost2 += dstMatrix[node][nb] + dstMatrix[actAns[(i+N-1)%N]][actAns[(jdx+N-1)%N]];  
            //Pongo el cambio de mejor costo
            if (cost1 < minCost){
                minCost = cost1;
                bestAns = std::make_pair(jdx, 0);
            } 
            if (cost2 < minCost){
                minCost = cost2;
                bestAns = std::make_pair(jdx, 2);
            }
        } while(nearNodes[node][nearIdx++] != fartestNb);
        //Cambiamos por el mejor costo y guardamos los pares de la mejor solución
        if (minCost < bestLocalcost){
            change = std::make_pair(i, bestAns.first);
            bestLocalcost = minCost;
            adyacent = bestAns.second;
        }
    } 
    if(bestLocalcost < cost) cost = bestLocalcost;
    //Regresamos nuestro vector de mejora
    return change;
}

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
                nodeToCoord[node-1] = coordinates;
            }
        }
    }
    file.close();
    return nodeToCoord;
}

lli ceil_2D(const std::pair<float,float> &x, const std::pair<float,float> &y){
    float first = x.first - y.first, second = x.second - y.second;

    //return int(ceil(sqrt(pow(first, 2) + pow(second, 2))));
    return int(sqrt(pow(first, 2) + pow(second, 2)));
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

void updateMap(std::unordered_map<int,int> &where, const std::vector<int> &ans){
    for(int i = 0; i < ans.size(); i++) where[ans[i]] = i;

    return;
}

void make_mxs(lliMat &dstMatrix, mat &nearNodes, const int &N, std::unordered_map<int,std::pair<float,float>> &nodesToCoord){
    for(int i = 0; i < N; i++){
        std::pair<float,float> x = nodesToCoord[i];

        for(int j = 0; j < N; j++){
            nearNodes[i].push_back(j);
            if(j == i) continue;
            std::pair<float,float> y = nodesToCoord[j];
            dstMatrix[i][j] = ceil_2D(x, y);
        }

        std::sort(nearNodes[i].begin(), nearNodes[i].end(), [&](int s, int j){return dstMatrix[i][s] < dstMatrix[i][j];});
    }

    return;
}

void esthocastic_hillClimbing(std::vector<int> &actualAns, lli &cost, const lliMat &dstMatrix, const mat &nearNodes){
    lli NewCost = cost;
    int N = actualAns.size(), adyacent = 0;
    std::unordered_map<int,int> where;

    int debug = 0;
    while(1){
        updateMap(where, actualAns);
        //Aplicamos la busqueda local mejorada
        std::pair<int,int> posChange = searchLocalImprove(dstMatrix, nearNodes, actualAns, NewCost, where, adyacent);
        
        if(posChange != std::make_pair(-1,-1)){
            //Establezco los indices en los que hago el cambio
            int idx = posChange.first, jdx = posChange.second;
            jdx += (jdx < idx) ? N:0;
            int intervalSize = abs(jdx - idx);
            //Hago los ajustes en el costo y en la solucion
            switch(adyacent){
            case 0:                               //Caso 0 (primer tipo de tranformacion)
                for(int i = 0; i < intervalSize/2; i++) std::swap(actualAns[(idx + i + 1)%N], actualAns[(jdx-i)%N]);
                break;
            case 2:                               //Caso 2 (segundo tipo de tranformacion)
                for(int i = 0; i < intervalSize/2; i++) std::swap(actualAns[(idx + i)%N], actualAns[(jdx-(i+1))%N]);
                break;
            default:                              //Caso 3 (adyacentes)
                std::swap(actualAns[idx],actualAns[jdx]);
                break;
            }

            cost = NewCost;
        } else break;
    }

    return;
}


int main(){
    std::vector<std::string> files = {"d493"};   //Nombre del archivo

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
            lliMat dstMatrix(sol.size(), std::vector<lli>(sol.size(), -1));
            mat nearNodes(sol.size());
            make_mxs(dstMatrix, nearNodes, sol.size(), nodeToCoord);
            esthocastic_hillClimbing(sol, cost, dstMatrix, nearNodes);
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
