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
    double density)
    : c(c), theta(theta), cvel(cvel), w(w), density(density), rbtemplate_(rbtemplate)
{
    AABB = buildAABB(this);
    voronois = rbtemplate.voronois;
    springs = rbtemplate.springs;
}

RigidBodyInstance::~RigidBodyInstance()
{    
    freeAABB(AABB);
    AABB = nullptr;
}

int RigidBodyInstance::lookupVoronoiFromTet(int tet){ return rbtemplate_.lookupVoronoiFromTet(tet); }

int RigidBodyInstance::lookupVoronoiFromVert(int vert) { return rbtemplate_.lookupVoronoiFromVert(vert); }

}
