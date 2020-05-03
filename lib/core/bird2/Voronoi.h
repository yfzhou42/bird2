#ifndef PSIM_CORE_BIRD2_VORONOI_H
#define PSIM_CORE_BIRD2_VORONOI_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include "Spring.h"

using namespace Eigen;
using namespace std;

class VoronoiPoint
{
private:
    /* data */
public:
    Eigen::MatrixX4i T;
    Eigen::Vector3d center;
    vector<Spring*> springs;
    Eigen::Vector3d Fc;

    void addForce(Eigen::Vector3d Fc){
        this->Fc += Fc;
    }

    VoronoiPoint(){
        T.resize(0, 4);
        center = Eigen::Vector3d::Zero();
        Fc = Eigen::Vector3d::Zero();
    }
    VoronoiPoint(Eigen::Vector3d center, Eigen::MatrixX4i T){
        this->T = T;
        this->center = center;
    }

    void mergeT(const VoronoiPoint& other){
        int tSize = T.rows();
        T.resize(tSize + other.T.rows(), 4);
        for(int i = 0; i < other.T.rows(); i++){
            T.row(tSize + i) = other.T.row(i);
        }
    }

    Eigen::MatrixX3d remappedVerts;
    Eigen::MatrixX4i remappedTets;

    void remapVerts(const Eigen::MatrixX3d& V) {
        vector<Vector3d> verts;
        vector<Vector4i> tets;
        map<int, int> indexList;
        for(int i = 0; i < T.rows(); i++){
            Vector4i tet;
            for(int j = 0; j < 4; j++) {
                if (indexList.find(T.row(i)(j)) != indexList.end()){
                    tet(j) = indexList[T.row(i)(j)];
                }
                else {
                    indexList[T.row(i)(j)] = verts.size();
                    verts.push_back(V.row(T.row(i)(j)));
                    tet(j) = verts.size()-1;
                }
            }
            tets.push_back(tet);
        }
        
        //Insert into remapped Matrices
        remappedVerts.resize(verts.size(), 3);
        for(int i = 0; i < verts.size(); i++){
            remappedVerts.row(i) = verts[i];
        }

        remappedTets.resize(tets.size(), 4);
        for(int i = 0; i < tets.size(); i++) {
            remappedTets.row(i) = tets[i];
        }
    }
};

#endif