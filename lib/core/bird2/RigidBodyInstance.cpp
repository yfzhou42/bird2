#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "RigidBodyTemplate.h"
#include <Eigen/Geometry>
#include <iostream>
#include "CollisionDetection.h"

using namespace Eigen;
using namespace std;

namespace bird2 {

RigidBodyInstance::RigidBodyInstance(const RigidBodyTemplate &rbtemplate,
    const Eigen::Vector3d &c, const Eigen::Vector3d &theta,
    const Eigen::Vector3d &cvel, const Eigen::Vector3d &w,
    double density, double maxStrain)
    : c(c), theta(theta), cvel(cvel), w(w), density(density), rbtemplate_(rbtemplate), maxStrain(maxStrain)
{
    voronois = rbtemplate.voronois;
    generateSprings();
    AABB = buildAABB(this);
}

RigidBodyInstance::~RigidBodyInstance()
{    
    freeAABB(AABB);
    AABB = nullptr;
}

int RigidBodyInstance::lookupVoronoiFromTet(int tet){ return rbtemplate_.lookupVoronoiFromTet(tet); }

int RigidBodyInstance::lookupVoronoiFromVert(int vert) { return rbtemplate_.lookupVoronoiFromVert(vert); }

void makeFaces(Vector4i tet, set<set<int>>& faces) {
    set<int> face1;
    face1.insert({tet[0], tet[1], tet[2]});
    faces.insert(face1);
    set<int> face2;
    face2.insert({tet[1], tet[2], tet[3]});
    faces.insert(face2);
    set<int> face3;
    face3.insert({tet[0], tet[1], tet[3]});
    faces.insert(face3);
    set<int> face4;
    face4.insert({tet[0], tet[2], tet[3]});
    faces.insert(face4);
}

void RigidBodyInstance::generateSprings() {
    /*
    for all v in voronoi
        for all v2 voronoi from v + 1
            if v and v2 share 3 vertices in a tet
                add a spring from v to v2
    */
    for (int i = 0; i < voronois.size(); i ++) {
        MatrixX4i T1 = voronois[i].T;
        //Set of all faces from tet
        set<set<int>> tetSet1;
        for (int t1 = 0; t1 < T1.rows(); t1++) {
            //for each tet which contains 4 faces
            set<set<int>> faces;
            makeFaces(T1.row(t1), faces);
            tetSet1.insert(faces.begin(), faces.end());
        }
        for (int j = i + 1; j < voronois.size(); j++) {
            MatrixX4i T2 = voronois[j].T;
            bool connected = false;
            for (int t2 = 0; t2 < T2.rows() && !connected; t2++) {
                set<set<int>> faces;
                makeFaces(T2.row(t2), faces);
                for(set<int> f : faces){
                    if(tetSet1.find(f) != tetSet1.end()) {
                        springs.push_back(new Spring(i, j, (voronois[i].center - voronois[j].center).squaredNorm()));
                        voronois[i].springs.push_back(springs.back());
                        voronois[j].springs.push_back(springs.back());
                        connected = true;
                        break;
                    }
                }
            }
        }
    }
}

}
