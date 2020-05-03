#ifndef PSIM_CORE_BIRD2_RIGIDBODYINSTANCE_H
#define PSIM_CORE_BIRD2_RIGIDBODYINSTANCE_H

#include <Eigen/Core>
#include <list>
#include <vector>
#include "Voronoi.h"
#include "Spring.h"
#include "RigidBodyTemplate.h"

namespace bird2 {

class RigidBodyTemplate;
struct AABBNode;

class RigidBodyInstance
{
public:
    RigidBodyInstance(const RigidBodyTemplate &rbtemplate, const Eigen::Vector3d &c, const Eigen::Vector3d &theta, const Eigen::Vector3d &cvel, const Eigen::Vector3d &w, double density);
    ~RigidBodyInstance();

    Eigen::Vector3d c;
    Eigen::Vector3d theta;

    Eigen::Vector3d cvel;
    Eigen::Vector3d w;

    std::vector<VoronoiPoint> voronois;
    std::vector<Spring> springs;

    double density;

    AABBNode *AABB;
    
    const RigidBodyTemplate &getTemplate() const {return rbtemplate_;}

    void zeroVornoiForces(){
        for (int i = 0; i < voronois.size(); i++)
            voronois[i].Fc = Vector3d::Zero();
    }
    
    int lookupVoronoiFromTet(int tet);
    int lookupVoronoiFromVert(int vert);

    int32_t bid;
private:
    const RigidBodyTemplate &rbtemplate_;
};

}

#endif // RIGIDBODYINSTANCE_H
