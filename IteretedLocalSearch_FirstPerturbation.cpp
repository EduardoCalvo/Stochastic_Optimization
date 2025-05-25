#include <bits/stdc++.h>
#include <sys/time.h>
#include <omp.h>

typedef long long int lli;
typedef std::vector<std::vector<lli>> lliMat;
typedef std::vector<std::vector<int>> mat;

bool validEdges(const std::pair<int,int> &edge1, const std::pair<int,int> &edge2){ 
    return (edge1.first != edge2.second) && (edge1.second != edge2.first); 
}

double CurrentExecutionTime(timeval start) {
    timeval current;
    gettimeofday(&current, NULL);
    return (current.tv_sec - start.tv_sec) + (current.tv_usec - start.tv_usec)*1e-6;
}

void valPairs(int &i, int &j, int &m, int &n, const std::vector<std::pair<int,int>> &edges){
    int toSwitch = 1;
    while(1){
        std::vector<int> rest;
        for(int i = 1; i < 4; i++) if (i != toSwitch) rest.push_back(i);

        if(validEdges(edges[0], edges[toSwitch]) && validEdges(edges[rest[0]], edges[rest[1]])){
            i = 0;
            j = toSwitch;
            m = rest[0];
            n = rest[1];

            break;
        }
        toSwitch++;
    }
    
    return;
}

void updateMap(std::unordered_map<int,int> &where, const std::vector<int> &ans){
    for(int i = 0; i < ans.size(); i++) where[ans[i]] = i;

    return;
}

lli ceil_2D(const std::pair<float,float> &x, const std::pair<float,float> &y){
    float first = x.first - y.first, second = x.second - y.second;

    return int(ceil(sqrt(pow(first, 2) + pow(second, 2))));
    // return int(sqrt(pow(first, 2) + pow(second, 2)));
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

std::vector<int> sampleMaker(const int &smp_size, const int &N) {
    std::set<int> randomSample;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, N-1);
    //Generamos los indices de los vecinos
    while (randomSample.size() < smp_size) randomSample.insert(dis(gen));
    return std::vector<int>(randomSample.begin(), randomSample.end());
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

lli computeCost(const std::vector<int> &ans, const lliMat &dstMat){
    lli cost = 0;
    for(int i = 0; i < ans.size(); i++){
        cost += dstMat[ans[i]][ans[(i+1)%ans.size()]];
    }

    return cost;
}

std::vector<int> first_pert(const float &alpha, std::vector<int> tour, const lliMat &dstMatrix, const mat &nearNodes, const lli &cost, lli &newCost){
    int N = tour.size();
    std::vector<int> perturbation = tour;
    std::vector<std::pair<int,int>> edgesToDel;
    std::unordered_map<int, int> where;
    updateMap(where, perturbation);
    newCost = cost;
    // Hacemos la eleccion de v
    int f_Node, f_inT;
    lli maxDif = -1;

    for(int i = 0; i < N; i++){
        int nodeV = tour[i];   
        // Calculamos la diferencia de distancias de su nodo más cercano con su siguiente en el tour
        lli pivot = dstMatrix[nodeV][tour[(i+1) % N]] - dstMatrix[nodeV][nearNodes[nodeV][1]];
        // Si es minima lo tomamos como nuevo vector inicial
        if(pivot > maxDif){
            f_inT = i;
            f_Node = tour[i];
            maxDif = pivot;
        }
    }
    edgesToDel.push_back(std::make_pair(tour[f_inT], tour[(f_inT+1)%N])); // Agregamos la primer arista a remover
    // Muestra aleatoria de ciudades
    int smp_size = (ceil(alpha*N) < 3) ? 3 : (int)ceil(alpha*N);
    std::vector<int> randSample = sampleMaker(smp_size, N);
    // Ordenar la muestra en base a las distancias con fNode
    std::sort(randSample.begin(), randSample.end(), [&](int i, int j){return dstMatrix[f_Node][tour[(f_inT+i)%N]] < dstMatrix[f_Node][tour[(f_inT+j)%N]];});
    // Tomar los 3 primeros para agregar los edges
    for(int i = 0; i < 3; i++){
        int aux = (f_inT + randSample[i])%N;
        edgesToDel.push_back(std::make_pair(tour[aux], tour[(aux+1)%N]));
    }

    int i = -1, j = -1, m = -1, n = -1;
    valPairs(i, j, m, n, edgesToDel);
    // Hacemos el primer cambio con las aristas i, j
    int idx = (where[edgesToDel[i].first] == (where[edgesToDel[i].second]+1)%N) ? where[edgesToDel[i].first] : where[edgesToDel[m].second];
    int jdx = (where[edgesToDel[j].first] == (where[edgesToDel[j].second]+1)%N) ? where[edgesToDel[j].second] : where[edgesToDel[j].first];
    // newCost -= dstMatrix[tour[idx]][tour[(idx+N-1)%N]] + dstMatrix[tour[jdx]][tour[(jdx+1)%N]];
    // newCost += dstMatrix[tour[jdx]][tour[(idx+N-1)%N]] + dstMatrix[tour[idx]][tour[(jdx+1)%N]];

    jdx += (jdx < idx) ? N:0;
    for(int i = 0; i < (jdx-idx+1)/2; i++) std::swap(perturbation[(idx + i)%N], perturbation[(jdx-i)%N]);
    updateMap(where, perturbation);
    // Segundo cambio con las aristas m,n
    idx = (where[edgesToDel[m].first] == (where[edgesToDel[m].second]+1)%N) ? where[edgesToDel[m].first] : where[edgesToDel[m].second];
    jdx = (where[edgesToDel[n].first] == (where[edgesToDel[n].second]+1)%N) ? where[edgesToDel[n].second] : where[edgesToDel[n].first];
    // newCost -= dstMatrix[tour[idx]][tour[(idx+N-1)%N]] + dstMatrix[tour[jdx]][tour[(jdx+1)%N]];
    // newCost += dstMatrix[tour[jdx]][tour[(idx+N-1)%N]] + dstMatrix[tour[idx]][tour[(jdx+1)%N]];

    jdx += (jdx < idx) ? N:0;
    for(int i = 0; i < (jdx-idx+1)/2; i++) std::swap(perturbation[(idx + i)%N], perturbation[(jdx-i)%N]);
    updateMap(where, perturbation);

    newCost = computeCost(perturbation, dstMatrix);
    return perturbation;
}

std::pair<int,int> localSearchImprove(const lliMat &dstMatrix, const mat &nearNodes, const std::vector<int> &tour, lli& cost, std::unordered_map<int,int> &where, int &adyacent){
    int N = tour.size();
    lli costOfChange = cost;
    std::pair<int,int> change = std::make_pair(-1,-1);

    for(int i = 0; i < N; i++){
        //std::cout << "it: " << i << std::endl; 
        std::pair<int,int> posChange = std::make_pair(-1,-1);
        int node = tour[i], solType;              // Nodo del tour sobre el que busco mejora
        lli bestNodeCost = costOfChange;          // Mejor costo de los cambios basados en el nodo actual 

        int dstNext = dstMatrix[node][tour[(i+1)%N]], dstPrev = dstMatrix[node][tour[(i-1+N)%N]]; // Distancias a sus vecinos
        int fartestNb = (dstNext < dstPrev) ? tour[(i-1+N)%N] : tour[(i+1)%N];                    // Vecino mas lejano
        
        int nearIdx = 1;
        while(nearIdx < nearNodes[node].size() && nearNodes[node][nearIdx] != fartestNb){
            int nb = nearNodes[node][nearIdx]; int nb_idx = where[nb];    // Elemento en la lista de nodos cercanos y posicion en tour
            lli cost1 = cost, cost2 = cost;                               // Costo de ambos cambios

            // Caso donde el nodo de la lista es adyacente
            if(nb_idx == (i+1)%N || nb_idx == (i-1+N)%N){  
                int nxtPrv = (nb_idx == (i+1)%N) ? 1 : -1;                // Checo si nb es el previo o el siguiente a node;
                
                // Checo la unica opcion cuando son adyacentes
                cost1 -= dstMatrix[node][tour[(i-nxtPrv+N)%N]] + dstMatrix[nb][tour[(nb_idx+nxtPrv+N)%N]];
                cost1 += dstMatrix[node][tour[(nb_idx+nxtPrv+N)%N]] + dstMatrix[nb][tour[(i-nxtPrv+N)%N]];

                // Si hay mejora en el costo entonces hacemos el cambio
                if (cost1 < bestNodeCost){
                    bestNodeCost = cost1;
                    posChange = std::make_pair(i, nb_idx);
                    solType = nxtPrv;                                     // Nos dice que son adyacentes y la relacion entre ellos
                }  

                // Saltamos a la siguiente iteración
                nearIdx++;
                continue;
            }

            // Checamo que indices es menor para aplicar el cambio
            int lower = (i < nb_idx) ? i : nb_idx;
            int upper = (i < nb_idx) ? nb_idx : i;

            // Primera opcion (swapear upper hacia atras)
            cost1 -= dstMatrix[tour[lower]][tour[(lower+1)%N]] + dstMatrix[tour[upper]][tour[(upper+1)%N]];
            cost1 += dstMatrix[tour[lower]][tour[upper]] + dstMatrix[tour[(lower+1)%N]][tour[(upper+1)%N]];
            // Segunda opcion (swapear lower hacia adelante)
            cost2 -= dstMatrix[tour[lower]][tour[(lower-1+N)%N]] + dstMatrix[tour[upper]][tour[(upper-1+N)%N]];
            cost2 += dstMatrix[tour[lower]][tour[upper]] + dstMatrix[tour[(lower-1+N)%N]][tour[(upper-1+N)%N]];

            // Si la op 1 mejora el mejor costo encontrado para ese nodo lo tomamos como el mejor
            if (cost1 < bestNodeCost) {
                bestNodeCost = cost1;
                posChange = std::make_pair(lower, upper);
                solType = 0;                                        // Para saber que op tomamos
            }
            // Si la op 2 mejora el mejor costo encontrado para ese nodo lo tomamos como el mejor
            if (cost2 < bestNodeCost) {
                bestNodeCost = cost2;  
                posChange = std::make_pair(lower, upper);
                solType = 2;                                        // Para saber que op tomamos
            }

            nearIdx++;
        }

        // Si el mejor costo para node supera al mejor costo encontrado sobre todos los nodos anteriores lo convertimos en el mejor
        if(bestNodeCost < costOfChange){                            
            costOfChange = bestNodeCost;
            change = posChange;
            adyacent = solType;
        }
    }

    // Si hubo mejora esto claramente bajo, entonces actualizamos
    if(costOfChange < cost) cost = costOfChange;

    return change;
}

void localSearch(std::vector<int> &tour, lli &cost, const lliMat &dstMtx, const mat &nearNodes){
    lli newCost = cost;
    int N = tour.size(), changeTyp;
    std::unordered_map<int,int> where;
    
    while(1){
        updateMap(where, tour);
        std::pair<int,int> change = localSearchImprove(dstMtx, nearNodes, tour, newCost, where, changeTyp);

        if (change != std::make_pair(-1,-1)){
            int low = change.first, up = change.second;
            int interSize = (up - low);

            //Hago el cambio segun los casos de localSearchImprove
            switch (changeTyp){
            case 0:
                for(int i = 0; i < interSize/2; i++) std::swap(tour[low + 1 + i], tour[up - i]);
                break;

            case 2:
                for(int i = 0; i < interSize/2; i++) std::swap(tour[low + i], tour[up - 1 - i]);
                break;

            default:
                std::swap(tour[low], tour[up]);
                break;
            }

            cost = newCost;
        } else break;
    }

    return;
}

void iterative_loc_search(std::vector<int> &act_ans, lli &cost, const lliMat &dstMat, const mat &nearNodes,const float &alpha){
    int N = act_ans.size();
    std::set<std::vector<int>> pert_history, sol_history;
    timeval start;
    gettimeofday(&start, NULL);

    localSearch(act_ans, cost, dstMat, nearNodes); // Busqueda local sobre actAns y cost

    while(CurrentExecutionTime(start) <= 180){
        lli pertCost = cost;
        std::vector<int> perturbation(N, 0); 
        perturbation = first_pert(alpha, act_ans, dstMat, nearNodes, cost, pertCost); //Perturbacion

        //Checamos el historial de perturbaciones
        if (pert_history.find(perturbation) == pert_history.end())
            pert_history.insert(perturbation);
        else {
            perturbation = first_pert(alpha, perturbation, dstMat, nearNodes, pertCost, pertCost); 
            pert_history.insert(perturbation);
        }

        //std::unordered_set<int> proof;
        //for(auto ele : perturbation) std::cout << ele << " ";
        //std::cout << std::endl;
        
        // for(auto ele : perturbation) proof.insert(ele);
        // if(proof.size() != N) {  
        //     for(auto ele : proof) std::cout << ele << " " ; 
        //     std::cout << proof.size() << std::endl;
        //     std::cout << "vales verga" << std::endl;
        //     break;
        // }

        //std::cout << "No es local search" << std::endl;

        localSearch(perturbation, pertCost, dstMat, nearNodes); //Busqueda local sobre perturbation y pertCost

        std::cout << "Termino local search: " << cost << std::endl;
        if(pertCost < 0) {
            std::cout << pertCost << " " << CurrentExecutionTime(start) << std::endl;
            break;
        }

        if(pertCost < cost){
        //if(pertCost < cost && sol_history.find(perturbation) == sol_history.end()){
        //  sol_history.insert(perturbation);
            act_ans = perturbation;
            cost = pertCost;
        }
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
        //Nombre del archivo de costos para la instancia filename
        std::string newFile = filename + ".txt";             
        FILE *archivo = fopen(newFile.c_str(), "w");        
        fprintf(archivo, "Costos con la primer solucion dada por la heuristica para la instancia: %s\n\n", filename.c_str());
        //Creamos matrices de distancia
        lliMat dstMatrix(nodeToCoord.size(), std::vector<lli>(nodeToCoord.size(), -1));
        mat nearNodes(nodeToCoord.size());
        make_mxs(dstMatrix, nearNodes, nodeToCoord.size(), nodeToCoord);

        for(int i = 0; i < 1; i++){
            lli newCost;
            timeval start;
            gettimeofday(&start, NULL);
            double time;
            //Heuristica
            farthenestInsertion(sol, nodeToCoord, cost);
            //first_pert(0.2, sol, dstMatrix, nearNodes, cost, newCost);
            //BusquedaLocalIterada
            iterative_loc_search(sol, cost, dstMatrix, nearNodes, 0.02);

            time = CurrentExecutionTime(start);

            fprintf(archivo, "Solución #%d\n", i+1);
            fprintf(archivo,"Costo: %lld\n", cost);
            fprintf(archivo, "Tiempo de ejecucion: %.5f\n\n", time);
        }

        fclose(archivo);
    }

    return 0;
}
