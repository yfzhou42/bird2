#ifndef PSIM_CORE_BIRD2_VORONOI_H
#define PSIM_CORE_BIRD2_VORONOI_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include "Spring.h"
#include <iostream>
#include "RigidBodyData.h"

using namespace Eigen;
using namespace std;

namespace bird2{

typedef struct FalseStep{
    Vector3d c;
    Vector3d v;
}FalseStep;

class VoronoiPoint
{
private:
    //If dirty then we need to remake rigidBodyData
    bool dirty;
    RigidBodyData* rigidBodyData;
    void updateRigidBodyData(){
        if(dirty){
            throw "RIGIDBODYDATA HASN'T BEEN UPDATED, PLEASE CALL remapVerts(V) BEFORE USING";
        }
    }
public:
    Eigen::MatrixX4i T;
    Eigen::Vector3d center;
    vector<Spring*> springs;
    Eigen::Vector3d Fc;
    Eigen::Vector3d Ftheta;
    FalseStep falseStep;

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
        dirty = true;
        rigidBodyData = NULL;
    }
    VoronoiPoint(Eigen::Vector3d center, Eigen::MatrixX4i T){
        this->T = T;
        this->center = center;
        Fc = Eigen::Vector3d::Zero();
        Ftheta = Eigen::Vector3d::Zero();
        dirty = true;
        rigidBodyData = NULL;
    }

    void merge(const VoronoiPoint& other){
        int tSize = T.rows();
        T.conservativeResize(tSize + other.T.rows(), 4);
        for(int i = 0; i < other.T.rows(); i++){
            T.row(tSize + i) = other.T.row(i);
        }
        Fc += other.Fc;
        Ftheta += other.Ftheta;
        dirty = true;
    }

    Eigen::Vector3d getCOM(){
        updateRigidBodyData();
        return rigidBodyData->getCenterOfMass();
    }

    Eigen::Matrix3d getInertiaTensor(){
        updateRigidBodyData();
        return rigidBodyData->getInertiaTensor();
    }

    double getVolume(){
        updateRigidBodyData();
        return rigidBodyData->getVolume();
    }

    std::pair<MatrixX4i, MatrixX3d> remapVerts(const Eigen::MatrixX3d& V) {
        if(!dirty){
            return std::make_pair(rigidBodyData->getTets(), rigidBodyData->getVerts());
        }
        vector<Vector3d> verts;
        vector<Vector4i> tets;
        map<int, int> indexList;

        Eigen::MatrixX3d remappedVerts;
        Eigen::MatrixX4i remappedTets;
        
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
        if(!rigidBodyData) delete rigidBodyData;
        rigidBodyData = new RigidBodyData(remappedVerts, remappedTets);
        dirty = false;
        return std::make_pair(remappedTets, remappedVerts);
    }
};

}

#endif