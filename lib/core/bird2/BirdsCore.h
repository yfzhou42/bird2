#ifndef PSIM_CORE_BIRD2_BIRDSCORE_H
#define PSIM_CORE_BIRD2_BIRDSCORE_H

#include "../PhysicsCore.h"
#include "SimParameters.h"
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <memory>
#include <set>

namespace bird2 {

class RigidBodyTemplate;
class RigidBodyInstance;
struct Collision;

class BirdsCore : public PhysicsCore {
public:
    BirdsCore();
    ~BirdsCore();

    virtual void initSimulation() override;

    virtual std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
    getCurrentMesh() const override;

    virtual bool simulateOneStep() override;

    std::shared_ptr<SimParameters> getPointerToSimParameters()
    {
        return params_;
    }

    void clearScene();

    /*
     * The facade of adding RigidBodyTemplate and its instances.
     */
    Eigen::VectorXi
    addMesh(const std::string& file_name,
            double scale,
            const Eigen::MatrixXd& Qs);

    /*
     * add an RigidBodyInstance directly
     */
    int32_t
    addSingleInstance(std::shared_ptr<RigidBodyTemplate> rbt,
                      double density,
                      const Eigen::Vector3d& c,
                      const Eigen::Vector3d& theta,
                      const Eigen::Vector3d& cvel,
                      const Eigen::Vector3d& w);

    int32_t
    addSingleInstanceStrain(std::shared_ptr<RigidBodyTemplate> rbt,
                      double density,
                      const Eigen::Vector3d& c,
                      const Eigen::Vector3d& theta,
                      const Eigen::Vector3d& cvel,
                      const Eigen::Vector3d& w,
                      const double maxStrain);

    /*
     * addInstances:
     *    add rigid body instance from Qs.
     *    Qs should be a (N, 13) matrix and each row
     */
    Eigen::VectorXi
    addInstances(std::shared_ptr<RigidBodyTemplate> rbt,
                 const Eigen::MatrixXd& Qs);

    //
    // Similar to Goo::query*, but we are returning a modifiable object.
    //
    // This interface change allows testing code to modify the phase directly.
    // Therefore the accumulation of numerical error during testing can be
    // eliminated by constantly synchronizing with the phase data from the reference code.
    //
    std::shared_ptr<RigidBodyInstance>
    queryRigidBodyInstance(int32_t bid);

    auto listRigidBodyInstances() const { return bodies_; }

    void computePenaltyCollisionForces(const std::set<Collision>& collisions,
                                       Eigen::Ref<Eigen::VectorXd> Fc,
                                       Eigen::Ref<Eigen::VectorXd> Ftheta);
    void applyCollisionImpulses(const std::set<Collision>& collisions);

private:
    int32_t rigid_body_id_;
    void computeForces(Eigen::VectorXd& Fc, Eigen::VectorXd& Ftheta);
    void shatter(int bodyIndex, std::set<int> brokenVoronoi);
    std::set<int> toShatter(int index);
    void breakVoronois(Eigen::Ref<Eigen::VectorXd> Fc, Eigen::Ref<Eigen::VectorXd> Ftheta);

    double time_;
    std::shared_ptr<SimParameters> params_;

    std::vector<std::shared_ptr<RigidBodyTemplate>> templates_;
    std::vector<std::shared_ptr<RigidBodyInstance>> bodies_;

    std::vector<Eigen::MatrixXd> init_configurations_;
};

}

#endif
