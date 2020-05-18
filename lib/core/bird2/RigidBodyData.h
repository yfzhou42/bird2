#ifndef PSIM_CORE_BIRD2_RIGIDBODYDATA_H
#define PSIM_CORE_BIRD2_RIGIDBODYDATA_H

#include <string>
#include <Eigen/Core>
#include <set>
#include <vector>



namespace bird2 {

class RigidBodyData
{
public:
    RigidBodyData() =delete;
    RigidBodyData(const Eigen::MatrixX3d &verts, const Eigen::MatrixX4i &tets);
    ~RigidBodyData();

    double getVolume() const {return volume_;}
    Eigen::Vector3d getCenterOfMass() { return com_; }
    const Eigen::Matrix3d getInertiaTensor() const { return inertiaTensor_; }
    
    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}      
    const Eigen::MatrixX4i &getTets() const { return T; }
protected:

    void initialize();
    
    void computeFaces();
    double computeVolume();
    Eigen::Vector3d computeCenterOfMass();
    Eigen::Matrix3d computeInertiaTensor();

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    Eigen::MatrixX4i T;
    
    double volume_;
    Eigen::Vector3d com_;  // Only used once, but kept for testing
    Eigen::Matrix3d inertiaTensor_;    

};

}


#endif