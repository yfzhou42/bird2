#include "BirdsCore.h"
#include "helper.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "CollisionDetection.h"
#include <iostream>
#include <Eigen/LU> // Required for .inverse()
#include <Eigen/Geometry> // Required for .cross()
#include <utility>

using namespace Eigen;

namespace bird2 {

BirdsCore::BirdsCore()
{
    params_ = std::make_shared<SimParameters>();
    initSimulation();
}

BirdsCore::~BirdsCore()
{
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd>
BirdsCore::getCurrentMesh() const
{
    int totverts = 0;
    int totfaces = 0;

    // floor
    totverts += 5;
    totfaces += 4;

    for (const auto& rbi : bodies_)
    {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;

    // floor
    double floory = -1.0;
    renderQ.row(0) << 0, floory, 0;
    renderQ.row(1) << 1e6, floory, 1e6;
    renderQ.row(2) << -1e6, floory, 1e6;
    renderQ.row(3) << -1e6, floory, -1e6;
    renderQ.row(4) << 1e6, floory, -1e6;
    voffset += 5;

    renderF.row(0) << 0, 2, 1;
    renderF.row(1) << 0, 3, 2;
    renderF.row(2) << 0, 4, 3;
    renderF.row(3) << 0, 1, 4;
    foffset += 4;

    for (const auto& rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
    // std::cerr << __func__ << " Nbodies " << bodies_.size() << std::endl;
    return std::make_tuple(renderQ, renderF, renderC);
}


void BirdsCore::initSimulation()
{
    rigid_body_id_ = 0;
    time_ = 0;
}

void BirdsCore::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc.resize(3*bodies_.size());
    Ftheta.resize(3*bodies_.size());
    Fc.setZero();
    Ftheta.setZero();

    if (params_->gravityEnabled) {
        for (int i=0; i<bodies_.size(); i++) {
            double m = bodies_[i]->density * bodies_[i]->getTemplate().getVolume();
            Fc[3 * i + 1] -= params_->gravityG*m;
        }
    }
}

bool BirdsCore::simulateOneStep()
{
    time_ += params_->timeStep;
    int nbodies = (int)bodies_.size();

    std::vector<Vector3d> oldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c += params_->timeStep*body.cvel;
        Matrix3d Rhw = VectorMath::rotationMatrix(params_->timeStep*body.w);
        Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

        // FIXME: angular velocity >= 3*PI
        Vector3d oldtheta = body.theta;
        body.theta = VectorMath::axisAngle(Rtheta*Rhw);
        if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI/2.0)
        {
            double oldnorm = oldtheta.norm();
            oldtheta = (oldnorm - 2.0*M_PI)*oldtheta/oldnorm;
        }

        oldthetas.push_back(oldtheta);
    }

    std::set<Collision> collisions;
    collisions = collisionDetection(bodies_);

    Eigen::VectorXd cForce(3 * nbodies);
    Eigen::VectorXd thetaForce(3 * nbodies);
    computeForces(cForce, thetaForce);

    computePenaltyCollisionForces(collisions, cForce, thetaForce);
    applyCollisionImpulses(collisions);

    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_->timeStep*cForce.segment<3>(3*bodyidx)/body.density/body.getTemplate().getVolume();

        Vector3d newwguess(body.w);

        int iter = 0;
        for(iter=0; iter<params_->NewtonMaxIters; iter++)
        {
            Vector3d term1 = (-VectorMath::TMatrix(-params_->timeStep*newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
            Vector3d term2 = (VectorMath::TMatrix(params_->timeStep*body.w).inverse()*VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
            Vector3d term3 = params_->timeStep * thetaForce.segment<3>(3*bodyidx);
            Vector3d fval = term1 + term2 + term3;
            if(fval.norm() / body.density / Mi.trace() <= params_->NewtonTolerance)
                break;

            Matrix3d Df = (-VectorMath::TMatrix(-params_->timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        // std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
        body.w = newwguess;
    }

    return false;
}

void
BirdsCore::clearScene()
{
    bodies_.clear();
    templates_.clear();
    init_configurations_.clear();
}

Eigen::VectorXi
BirdsCore::addMesh(const std::string& file_name,
                   double scale,
                   const Eigen::MatrixXd& Qs)
{
    auto tup = bird1::loadOBJ(file_name);
    templates_.emplace_back(new RigidBodyTemplate(std::get<0>(tup), std::get<1>(tup), scale));

    auto rbt = templates_.back();
    init_configurations_.emplace_back(Qs);
    return addInstances(rbt, Qs);
}

std::shared_ptr<RigidBodyInstance>
BirdsCore::queryRigidBodyInstance(int32_t bid)
{
    for (auto& b : bodies_)
        if (b->bid == bid)
            return b;
    return std::shared_ptr<RigidBodyInstance>(nullptr);
}

int32_t
BirdsCore::addSingleInstance(std::shared_ptr<RigidBodyTemplate> rbt,
                             double density,
                             const Eigen::Vector3d &c,
                             const Eigen::Vector3d &theta,
                             const Eigen::Vector3d &cvel,
                             const Eigen::Vector3d &w)
{
    bodies_.emplace_back(new RigidBodyInstance(*rbt, c, theta, cvel, w, density));
    bodies_.back()->bid = rigid_body_id_++;
    return bodies_.back()->bid;
}

Eigen::VectorXi
BirdsCore::addInstances(std::shared_ptr<RigidBodyTemplate> rbt,
                        const Eigen::MatrixXd& Qs)
{
    Eigen::VectorXi ret;
    ret.resize(Qs.rows());
    for (int i = 0; i < Qs.rows(); i++) {
        double density;
        Eigen::Vector3d c, cvel, theta, w;
        density = Qs(i, 0);
        int base = 1;
        c << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        theta << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        cvel << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        base += 3;
        w << Qs(i, base + 0), Qs(i, base + 1), Qs(i, base + 2);
        ret(i) = addSingleInstance(rbt, density, c, theta, cvel, w);
    }
    return ret;
}

void BirdsCore::computePenaltyCollisionForces(const std::set<Collision>& collisions,
                                              Eigen::Ref<VectorXd> Fc,
                                              Eigen::Ref<VectorXd> Ftheta)
{
    if (!params_->penaltyEnabled)
	    return;

    // TODO compute and add penalty forces
    for(Collision c : collisions){
        double dist;
        Vector3d derivDist;
        Matrix3d rotB1 = VectorMath::rotationMatrix(bodies_[c.body1]->theta);
        Matrix3d rotB2;

        Vector3d vertHit = bodies_[c.body1]->getTemplate().getVerts().row(c.collidingVertex);

        if(c.body2 == -1){
            rotB2 = Matrix3d::Identity();
            //Dist to floor in world coord
            dist = (rotB1 * vertHit + bodies_[c.body1]->c)(1) + 1.0;
            derivDist = Vector3d(0.0, 1.0, 0.0);
        }
        else {
            rotB2 = VectorMath::rotationMatrix(bodies_[c.body2]->theta);
            Vector3d phi = rotB2.transpose() * (rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c); 

            dist = bodies_[c.body2]->getTemplate().distance(phi, c.collidingTet);
            derivDist = bodies_[c.body2]->getTemplate().Ddistance(c.collidingTet);
        }

        Vector3d term1 = params_->penaltyStiffness * dist * derivDist;
        
        if(c.body2 != -1)
            Fc.segment<3>(c.body2 * 3) -= term1.transpose() * -rotB2.transpose();
        Fc.segment<3>(c.body1 * 3) -= term1.transpose() * rotB2.transpose();
        
        if(c.body2 != -1)
            Ftheta.segment<3>(c.body2 * 3) -= term1.transpose() * rotB2.transpose() * VectorMath::crossProductMatrix(rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c)
                * VectorMath::TMatrix(-bodies_[c.body2]->theta);
        
        Ftheta.segment<3>(c.body1 * 3) -= term1.transpose() * rotB2.transpose() * (-rotB1 * VectorMath::crossProductMatrix(vertHit)
            * VectorMath::TMatrix(bodies_[c.body1]->theta));
    
    }
}

typedef struct CollisionMeta{
    int body1;
    int body2;
    int collidingVertex;
    int collidingTet;
    double dist;
    Vector3d derivDist;
    double relativeVel;
    Matrix3d rotB1;
    Matrix3d rotB2;
    CollisionMeta(){}
    CollisionMeta(Collision c, double dist, Vector3d derivDist, double relativeVel, Matrix3d rotB1, Matrix3d rotB2) : 
                dist(dist), derivDist(derivDist), relativeVel(relativeVel), rotB1(rotB1), rotB2(rotB2) {
                    this->body1 = c.body1;
                    this->body2 = c.body2;
                    this->collidingTet = c.collidingTet;
                    this->collidingVertex = c.collidingVertex;
                }
} CollisionMeta;

void BirdsCore::applyCollisionImpulses(const std::set<Collision>& collisions)
{
    if (!params_->impulsesEnabled)
        return;

    // TODO apply collision impulses
    // Map of only the largest impules per pair.
    std::map<std::set<int>, CollisionMeta*> toApply;
    for(Collision c : collisions){
        if(c.body2 == -1) continue;
        Matrix3d rotB1 = VectorMath::rotationMatrix(bodies_[c.body1]->theta);
        Matrix3d rotB2 = VectorMath::rotationMatrix(bodies_[c.body2]->theta);
        Vector3d vertHit = bodies_[c.body1]->getTemplate().getVerts().row(c.collidingVertex);
        Vector3d phi = rotB2.transpose() * (rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c); 
        Vector3d c1 = bodies_[c.body1]->c;
        Vector3d c2 = bodies_[c.body2]->c;

        Vector3d derivDist = bodies_[c.body2]->getTemplate().Ddistance(c.collidingTet);

        Vector3d dPhi = VectorMath::crossProductMatrix(bodies_[c.body2]->w).transpose() * rotB2.transpose() * (rotB1 * vertHit + c1 - c2)
                + rotB2.transpose() * (rotB1 * VectorMath::crossProductMatrix(bodies_[c.body1]->w) * vertHit 
                + (bodies_[c.body1]->cvel) - (bodies_[c.body2]->cvel));

        double relativeVel =  derivDist.transpose() * dPhi;
        if(relativeVel > 0) continue;
        double dist = bodies_[c.body2]->getTemplate().distance(phi, c.collidingTet);
        if(std::isnan(dist) || dist >= 0.0) continue;
        std::set<int> collisionSet;
        collisionSet.insert(c.body2);
        collisionSet.insert(c.body1);
        try {
            std::map<std::set<int>, CollisionMeta*>::iterator it = toApply.find(collisionSet);
            if(it == toApply.end()) throw std::out_of_range("Not existent already");
            CollisionMeta* toComp = it->second;
            if(toComp->dist > dist){
                toApply[collisionSet] = new CollisionMeta(c, dist, derivDist, relativeVel, rotB1, rotB2);
                free(toComp);
            }
        }
        catch(const std::out_of_range& e){
            toApply[collisionSet] = new CollisionMeta(c, dist, derivDist, relativeVel, rotB1, rotB2);
        }
    }

    //Apply Impulses
    for(std::map<std::set<int>, CollisionMeta*>::iterator it = toApply.begin(); it != toApply.end(); it++){
        CollisionMeta c = *it->second;
        //std::cout << "Bodies: " << *it->first.begin() << " " << *(++(it->first.begin())) << " " << "Signed Dist: " << it->second->dist << std::endl;
        //std::cout << "Coll B: " << c.body1 << " " << c.body2 << " " << "Signed Dist: " << c.dist << std::endl;
        Matrix3d Minv1 = Matrix3d::Identity() / (bodies_[c.body1]->density* bodies_[c.body1]->getTemplate().getVolume());
        Matrix3d Minv2 = Matrix3d::Identity() / (bodies_[c.body2]->density* bodies_[c.body2]->getTemplate().getVolume());
        Vector3d vertHit = bodies_[c.body1]->getTemplate().getVerts().row(c.collidingVertex);
        Vector3d c1 = bodies_[c.body1]->c;
        Vector3d c2 = bodies_[c.body2]->c;

        Vector3d dPhi = VectorMath::crossProductMatrix(bodies_[c.body2]->w).transpose() * c.rotB2.transpose() * (c.rotB1 * vertHit + c1 - c2)
                + c.rotB2.transpose() * (c.rotB1 * VectorMath::crossProductMatrix(bodies_[c.body1]->w) * vertHit 
                + (bodies_[c.body1]->cvel) - (bodies_[c.body2]->cvel));

        Matrix3d A = c.rotB2.transpose() * VectorMath::crossProductMatrix(c.rotB1 * vertHit + c1 - c2) * VectorMath::TMatrix(-bodies_[c.body2]->theta);
        Matrix3d B = c.rotB2.transpose() * (-c.rotB1 * VectorMath::crossProductMatrix(vertHit) * VectorMath::TMatrix(bodies_[c.body1]->theta));
        Matrix3d C = -c.rotB2.transpose();
        Matrix3d D = c.rotB2.transpose();

        double relativeVel = c.derivDist.transpose() * dPhi;
        if(relativeVel > 0) continue;

        //dg calcs
        Vector3d dgThetaI = c.derivDist.transpose() * A;
        Vector3d dgThetaJ = c.derivDist.transpose() * B;
        Vector3d dgCI = c.derivDist.transpose() * C;
        Vector3d dgCJ = c.derivDist.transpose() * D;

        Matrix3d IInertiaInverse = bodies_[c.body2]->getTemplate().getInertiaTensor().inverse();
        Matrix3d JInertiaInverse = bodies_[c.body1]->getTemplate().getInertiaTensor().inverse();

        double alpha = -(1.0 + params_->CoR) * relativeVel;

        alpha /= c.derivDist.transpose() * (VectorMath::crossProductMatrix(c.rotB1 * vertHit + c1 - c2) * (dgThetaI.transpose() * VectorMath::TMatrix(bodies_[c.body2]->theta).inverse() * IInertiaInverse).transpose()
                    + c.rotB2.transpose() * c.rotB1 * VectorMath::crossProductMatrix((dgThetaJ.transpose() * VectorMath::TMatrix(bodies_[c.body1]->theta).inverse() * JInertiaInverse).transpose()) * vertHit
                    + Minv1 * dgCJ - Minv2 * dgCI);

        bodies_[c.body2]->cvel += Minv2 * (alpha * dgCI);
        bodies_[c.body2]->w = (bodies_[c.body2]->w.transpose() + (alpha * dgThetaI.transpose()) * VectorMath::TMatrix(bodies_[c.body2]->theta).inverse() * IInertiaInverse).transpose();
    
        bodies_[c.body1]->cvel += Minv1 * (alpha * dgCJ);
        bodies_[c.body1]->w += ((alpha * dgThetaJ.transpose()) * VectorMath::TMatrix(bodies_[c.body1]->theta).inverse() * JInertiaInverse).transpose();

        free(it->second);
    }
}

}
