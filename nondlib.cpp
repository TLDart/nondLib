#include "nondlib.hpp"

template<class T>
std::vector<std::vector<T>> generatePointset(int sz, int dim){
    std::vector<std::vector<T>> pointset;
    for(int i = 0; i < sz; i++){
        std::vector<T> row;
        for(int j = 0; j < dim; j++){
            row.push_back(((T) rand() / (RAND_MAX)));
        }
       pointset.push_back(row);

    }
    return pointset;
}

int main(){
    std::vector<std::vector<double>> real =  generatePointset<double>(100,3);;
    std::vector<std::vector<double>> real2(real);
    std::vector<double> mx = {1,1,1};

    filterDimSweep3D<std::vector<double>, std::vector<double>>(real,mx,0);
    //computeMaxima2D<std::vector<double>, std::vector<double>>(real, mx);
    //updateMaxima3D(test,mx,0, newpoint);
    for(int i = 0; i < real.size(); i++){
        for(int j = 0; j < real[0].size();  j++){
            std::cout << real[i][j] << " " ;
        } 
            std::cout << std::endl;
   }
    filterQuadD<std::vector<double>, std::vector<double>>(real2, mx);
}