#ifndef PSIM_CORE_BIRD2_VORONOI_H
#define PSIM_CORE_BIRD2_VORONOI_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include "Spring.h"
#include <iostream>

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
    Eigen::Vector3d Ftheta;

    void addForce(Eigen::Vector3d Fc){
        this->Fc += Fc;
    }

    void addThetaForce(Eigen::Vector3d Ftheta){
        this->Ftheta = Ftheta;
    }

    VoronoiPoint(){
        T.resize(0, 4);
        center = Eigen::Vector3d::Zero();
        Fc = Eigen::Vector3d::Zero();
        Ftheta = Eigen::Vector3d::Zero();
    }
    VoronoiPoint(Eigen::Vector3d center, Eigen::MatrixX4i T){
        this->T = T;
        this->center = center;
        Fc = Eigen::Vector3d::Zero();
        Ftheta = Eigen::Vector3d::Zero();
    }

    void merge(const VoronoiPoint& other){
        int tSize = T.rows();
        T.conservativeResize(tSize + other.T.rows(), 4);
        for(int i = 0; i < other.T.rows(); i++){
            T.row(tSize + i) = other.T.row(i);
        }
        Fc += other.Fc;
        Ftheta += other.Ftheta;
    }

    std::pair<MatrixX4i, MatrixX3d> remapVerts(const Eigen::MatrixX3d& V) {
        vector<Vector3d> verts;
        vector<Vector4i> tets;
        map<int, int> indexList;

        Eigen::MatrixX3d remappedVerts;
        Eigen::MatrixX4i remappedTets;
        
        for(int i = 0; i < T.rows(); i++){
            Vector4i tet;
            for(int j = 0; j < 4; j++) {
                if (indexList.find(T.row(i)(j)) != indexList.end()){
                    //cout << __LINE__ << endl;
                    tet(j) = indexList[T.row(i)(j)];
                }
                else {
                    //cout << __LINE__ << endl;
                    indexList[T.row(i)(j)] = verts.size();
                    /*cout << __LINE__ << endl;
                    cout << "i: " << i << " j: " << j << endl;
                    cout << "T.row(i)(j) " << T.row(i)(j) << endl;
                    cout << "V.row...: " << V.row(T.row(i)(j)) << endl;*/
                    verts.push_back(V.row(T.row(i)(j)));
                    //cout << __LINE__ << endl;
                    tet(j) = verts.size()-1;
                    //cout << __LINE__ << endl;
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

        return std::make_pair(remappedTets, remappedVerts);
    }
};

#endif