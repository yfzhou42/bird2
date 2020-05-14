#ifndef PSIM_CORE_BIRD2_RIGIDBODYTEMPLATE_H
#define PSIM_CORE_BIRD2_RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <set>
#include <vector>
#include "Voronoi.h"
#include "Spring.h"

#define VORONOI_POINTS 10

namespace bird2 {

class SignedDistanceField;

class RigidBodyTemplate
{
public:
    RigidBodyTemplate(Eigen::Ref<Eigen::MatrixX3d> V,
                      Eigen::Ref<Eigen::MatrixX3i> F,
                      double scale);
    RigidBodyTemplate(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets);
    ~RigidBodyTemplate();

    double getVolume() const {return volume_;}
    Eigen::Vector3d getCenterOfMass() { return com_; }
    const Eigen::Matrix3d getInertiaTensor() const { return inertiaTensor_; }
    Eigen::VectorXd getDistances() const { return distances_; }
    
    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}      
    const Eigen::MatrixX4i &getTets() const { return T; }

    double distance(Eigen::Vector3d p, int tet) const;
    Eigen::Vector3d Ddistance(int tet) const;

    std::vector<VoronoiPoint> voronois;
    std::vector<Spring*> springs;
    
    int lookupVoronoiFromTet(int tet) const { return tetToVoronoi[tet]; }
    int lookupVoronoiFromVert(int vert) const { return vertToVoronoi[vert]; }
private:
    RigidBodyTemplate(const RigidBodyTemplate &other) = delete;
    RigidBodyTemplate &operator=(const RigidBodyTemplate &other) = delete;

    void initialize();
    
    void computeFaces();
    double computeVolume();
    Eigen::Vector3d computeCenterOfMass();
    Eigen::Matrix3d computeInertiaTensor();
    Eigen::VectorXd computeDistances();
    void generateVoronoiPoints();
    void generateSprings();
    void generateVoronoiReferences();

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::MatrixX4i T;

    Eigen::VectorXd distances_;
    
    //mapping the tetrahedrons to the index of the voronois point in the voronois list
    std::vector<int> tetToVoronoi;
    //mapping the vertices to the index of the voronois point in the voronois list
    std::vector<int> vertToVoronoi;

    
    double volume_;
    Eigen::Vector3d com_;  // Only used once, but kept for testing
    Eigen::Matrix3d inertiaTensor_;    

};

}

#endif // PSIM_CORE_BIRD2_RIGIDBODYTEMPLATE_H
