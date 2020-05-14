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

int shattered = 0;

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
    if (params_->gravityEnabled) {
        for (int i=0; i<bodies_.size(); i++) {
            double m = bodies_[i]->density * bodies_[i]->getTemplate().getVolume();
            Fc[3 * i + 1] -= params_->gravityG*m;
        }
    }
}
//decide to shatter or not 
set<int> BirdsCore::toShatter(int index) {
    set<int> goingToShatter;
    /*
    for all voronoi in body i
        does voronoi violate contraint?
            cut all connected springs.
            
    if any springs are cut, return true
    */
   std::cout << __LINE__ << std::endl;

    if(bodies_[index]->generation >= params_->maxGeneration) return goingToShatter;

    vector<Vector3d> decompForces(bodies_[index]->voronois.size());
    for (int i = 0; i < bodies_[index]->voronois.size(); i++) {
        VoronoiPoint* vp = &bodies_[index]->voronois[i];
        if(vp->Fc.squaredNorm() < 0.0001) continue;
        for(Spring* s : vp->springs){
            Vector3d decomp_Fc;
            if(s->p1 == i){
                Vector3d springVec = vp->center - bodies_[index]->voronois[s->p2].center;
                decomp_Fc = (vp->Fc.dot(springVec) / springVec.squaredNorm()) * springVec;
            }
            else{
                Vector3d springVec = vp->center - bodies_[index]->voronois[s->p1].center;
                decomp_Fc = (vp->Fc.dot(springVec) / springVec.squaredNorm()) * springVec;
            }
            vp->addForce(-decomp_Fc);
        }
        if(vp->Fc.squaredNorm() > pow(bodies_[index]->maxStrain, 2)) {
            goingToShatter.insert(i);
            //add an extra force for the adjacent voronoi points
            /*for(Spring* s : vp->springs){
                Vector3d decomp_Fc;
                if(s->p1 == i){
                    Vector3d springVec = vp->center - bodies_[index]->voronois[s->p2].center;
                    decomp_Fc = (vp->Fc.dot(springVec) / springVec.squaredNorm()) * springVec;
                    //add an extra force for the adjacent voronoi points
                    bodies_[index]->voronois[s->p2].addForce(decomp_Fc);
                }
                else{
                    Vector3d springVec = vp->center - bodies_[index]->voronois[s->p1].center;
                    decomp_Fc = (vp->Fc.dot(springVec) / springVec.squaredNorm()) * springVec;
                    //add an extra force for the adjacent voronoi points
                    bodies_[index]->voronois[s->p1].addForce(decomp_Fc);
                }
            }*/
        }
    }
   std::cout << __LINE__ << std::endl;
    return goingToShatter;
}

void BirdsCore::breakVoronois(Eigen::Ref<VectorXd> Fc, Eigen::Ref<VectorXd> Ftheta) {
    if(!params_->fractureEnabled) return;
    int sizeToCheck = bodies_.size();
   std::cout << __LINE__ << std::endl;
    //Need this to index into Force vector
    int originalIndex = 0;
    for (int i = 0 ; i < sizeToCheck; i++, originalIndex++) {
        set<int> brokenVoronoi = toShatter(i);
        if(brokenVoronoi.size() > 0){
            // Based on spring network construct new bodies.
            shatter(i, brokenVoronoi);
            i--;
            sizeToCheck--;
        }
    }
   std::cout << __LINE__ << std::endl;

    // Update force vectors to reflect the new rigid body list
    Fc.resize(3 * bodies_.size());
    Ftheta.resize(3 * bodies_.size());
    for (int i = 0; i < bodies_.size(); i++) {
        Vector3d bodyFc = Vector3d::Zero();
        Vector3d bodyFtheta = Vector3d::Zero();
        for (VoronoiPoint vp : bodies_[i]->voronois) {
            bodyFc += vp.Fc;
            bodyFtheta += vp.Ftheta;
        }
        Fc.segment<3>(i * 3) = bodyFc;
        Ftheta.segment<3>(i * 3) = bodyFtheta;
    }
   std::cout << __LINE__ << std::endl;
}

typedef struct Node{
    set<int> ids;
    //check if two Nodes share an id
    bool overlap(Node other){
        for(int i : other.ids){
            if(ids.find(i) != ids.end()) return true;
        }
        return false;
    }
} Node;

typedef struct Tree{
    vector<Node> nodes;

    //Construct list of Nodes based on the Voronoi p
    Tree(vector<VoronoiPoint> vps){
        for(int i = 0 ; i < vps.size(); i++){
            VoronoiPoint vp = vps[i];
            Node n;
            n.ids.insert(i);
            for(Spring* s : vp.springs){
                if(s->p1 == -1 || s->p2 == -1) continue;
                n.ids.insert(s->p1);
                n.ids.insert(s->p2);
            }
            nodes.push_back(n);
        }
        merge();
    }
    //Merge all nodes together to create a forest
    void merge(){
        for(int i = 0; i < nodes.size(); i++){
            for(int j = i + 1; j < nodes.size(); j++) {
                if(nodes[i].overlap(nodes[j])){
                    nodes[i].ids.insert(nodes[j].ids.begin(), nodes[j].ids.end());
                    nodes.erase(nodes.begin() + j);
                    return merge();
                }
            }
        }
    }
} Tree;

void BirdsCore::shatter(int bodyIndex, set<int> brokenVoronoi) {
    /*
    When a rigid body breaks, the subbodies need to inherit the parent's position
        pos = parentPos + parentRot * voronoiCenterOfMass;
        rot = parentRot
        cvel = parentCVel
        w = ? parentW for now
    */
   std::cout << __LINE__ << std::endl;
    shared_ptr<RigidBodyInstance> b = bodies_[bodyIndex];
    // Update spring network to remove broken springs.
    for (int i = 0; i < b->springs.size(); i++) {
        if(brokenVoronoi.find(b->springs[i]->p1) != brokenVoronoi.end() || brokenVoronoi.find(b->springs[i]->p2) != brokenVoronoi.end()){
            b->springs[i]->p1 = -1;
            b->springs[i]->p2 = -1;
            b->springs.erase(b->springs.begin() + i);
            i--;
        }
    } 
    
    //Create new VoronoiPoints through springs that are still connected
    Tree voronoiIndexes(b->voronois);
    vector<VoronoiPoint> newBodies;
    for(int i = 0; i < voronoiIndexes.nodes.size(); i++){
        VoronoiPoint vp;
        for(int vIndex : voronoiIndexes.nodes[i].ids){
            vp.merge(b->voronois[vIndex]);
        }
        newBodies.push_back(vp);
    }
    
    for (VoronoiPoint vp : newBodies) {
        //Make a template for the voronoi fragment
        vp.remapVerts(bodies_[bodyIndex]->getTemplate().getVerts());
        templates_.emplace_back(new RigidBodyTemplate(vp.remappedVerts,vp.remappedTets));
        Vector3d pos = bodies_[bodyIndex]->c + VectorMath::rotationMatrix(bodies_[bodyIndex]->theta) * templates_.back()->getCenterOfMass();
        addSingleInstanceStrain(templates_.back(), bodies_[bodyIndex]->density, pos, bodies_[bodyIndex]->theta, bodies_[bodyIndex]->cvel, bodies_[bodyIndex]->w, bodies_[bodyIndex]->maxStrain);
        bodies_.back()->generation = bodies_[bodyIndex]->generation + 1;
        bodies_.back()->voronois[0].Fc = vp.Fc;
        bodies_.back()->voronois[0].Ftheta = vp.Ftheta;
    }
    bodies_.erase(bodies_.begin() + bodyIndex);
   std::cout << __LINE__ << std::endl;
}

bool BirdsCore::simulateOneStep()
{
        std::cout<<__LINE__<<std::endl;
    time_ += params_->timeStep;
    int nbodies = (int)bodies_.size();

    std::vector<Vector3d> oldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.zeroVornoiForces();
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

    computePenaltyCollisionForces(collisions, cForce, thetaForce);
    applyCollisionImpulses(collisions);

    breakVoronois(cForce, thetaForce);

    computeForces(cForce, thetaForce);

    if(nbodies != bodies_.size()){
        nbodies = bodies_.size();
        oldthetas.resize(nbodies);
        for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
        {
            RigidBodyInstance &body = *bodies_[bodyidx];
            //body.zeroVornoiForces();
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
            oldthetas[3] = oldtheta;
        }
    }
    

    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        std::cout<<__LINE__<<std::endl;
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
        std::cout<<__LINE__<<std::endl;


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
    bodies_.emplace_back(new RigidBodyInstance(*rbt, c, theta, cvel, w, density, params_->maxStrain));
    bodies_.back()->bid = rigid_body_id_++;
    bodies_.back()->generation = 0;
    return bodies_.back()->bid;
}

int32_t
BirdsCore::addSingleInstanceStrain(std::shared_ptr<RigidBodyTemplate> rbt,
                             double density,
                             const Eigen::Vector3d &c,
                             const Eigen::Vector3d &theta,
                             const Eigen::Vector3d &cvel,
                             const Eigen::Vector3d &w,
                             const double maxStrain)
{
    bodies_.emplace_back(new RigidBodyInstance(*rbt, c, theta, cvel, w, density, maxStrain));
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
        //std::cout << c.body1 <<c.body2 << std::endl;
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
            //std::cout << "NON FLOOR COLLISION" << std::endl;
            rotB2 = VectorMath::rotationMatrix(bodies_[c.body2]->theta);
            Vector3d phi = rotB2.transpose() * (rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c); 

            dist = bodies_[c.body2]->getTemplate().distance(phi, c.collidingTet);
            derivDist = bodies_[c.body2]->getTemplate().Ddistance(c.collidingTet);
        }

        Vector3d term1 = params_->penaltyStiffness * dist * derivDist;
        /*
        if(c.body2 != -1){
            Fc.segment<3>(c.body2 * 3) -= term1.transpose() * -rotB2.transpose();
        }
        Fc.segment<3>(c.body1 * 3) -= term1.transpose() * rotB2.transpose();
        
        if(c.body2 != -1){
            Ftheta.segment<3>(c.body2 * 3) -= term1.transpose() * rotB2.transpose() * VectorMath::crossProductMatrix(rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c)
                * VectorMath::TMatrix(-bodies_[c.body2]->theta);
        }
        
        Ftheta.segment<3>(c.body1 * 3) -= term1.transpose() * rotB2.transpose() * (-rotB1 * VectorMath::crossProductMatrix(vertHit)
            * VectorMath::TMatrix(bodies_[c.body1]->theta));
        */
        // Update both voronoi points with Fc.
        if(c.body2 != -1){
            int v2 = bodies_[c.body2]->lookupVoronoiFromTet(c.collidingTet);
            bodies_[c.body2]->voronois[v2].addForce(-term1.transpose() * -rotB2.transpose());
            bodies_[c.body2]->voronois[v2].addThetaForce(-term1.transpose() * rotB2.transpose() * VectorMath::crossProductMatrix(rotB1 * vertHit + bodies_[c.body1]->c - bodies_[c.body2]->c)
                * VectorMath::TMatrix(-bodies_[c.body2]->theta));
        }
        int v1 = bodies_[c.body1]->lookupVoronoiFromVert(c.collidingVertex);
        bodies_[c.body1]->voronois[v1].addForce(-term1.transpose() * rotB2.transpose());
        bodies_[c.body1]->voronois[v1].addThetaForce(-term1.transpose() * rotB2.transpose() * (-rotB1 * VectorMath::crossProductMatrix(vertHit)
            * VectorMath::TMatrix(bodies_[c.body1]->theta)));
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
        Matrix3d rotB1 = VectorMath::rotationMatrix(bodies_[c.body1]->theta);
        Matrix3d rotB2 = Matrix3d::Identity();
        Vector3d vertHit = bodies_[c.body1]->getTemplate().getVerts().row(c.collidingVertex);
        Vector3d c1 = bodies_[c.body1]->c;
        Vector3d c2 = Vector3d(0.0,-1.0,0.0);
        Vector3d w1 = bodies_[c.body1]->w;
        Vector3d w2 = Vector3d::Zero();
        Vector3d cVel1 = bodies_[c.body1]->cvel;
        Vector3d cVel2 = Vector3d::Zero();

        
        Vector3d derivDist = Vector3d(0.0, 1.0, 0.0);

        if(c.body2 != -1){
            derivDist = bodies_[c.body2]->getTemplate().Ddistance(c.collidingTet);
            rotB2 = VectorMath::rotationMatrix(bodies_[c.body2]->theta);
            c2 = bodies_[c.body2]->c;
            cVel2 = bodies_[c.body2]->cvel;
            w2 = bodies_[c.body2]->w;
        }

        Vector3d phi = rotB2.transpose() * (rotB1 * vertHit + c1 - c2); 

        Vector3d dPhi = VectorMath::crossProductMatrix(w2).transpose() * rotB2.transpose() * (rotB1 * vertHit + c1 - c2)
                + rotB2.transpose() * (rotB1 * VectorMath::crossProductMatrix(w1) * vertHit 
                + (cVel1) - (cVel2));

        double relativeVel =  derivDist.transpose() * dPhi;

        if(relativeVel > 0) continue;
        double dist = (rotB1 * vertHit + c1)(1) + 1.0;
        if(c.body2 != -1) dist = bodies_[c.body2]->getTemplate().distance(phi, c.collidingTet);
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
        Matrix3d Minv1 = Matrix3d::Identity() / (bodies_[c.body1]->density* bodies_[c.body1]->getTemplate().getVolume());
        Matrix3d Minv2 = Matrix3d::Zero();
        Vector3d vertHit = bodies_[c.body1]->getTemplate().getVerts().row(c.collidingVertex);
        Vector3d c1 = bodies_[c.body1]->c;
        Vector3d c2 = Vector3d(0.0,-1.0,0.0);
        Vector3d w1 = bodies_[c.body1]->w;
        Vector3d w2 = Vector3d::Zero();
        Vector3d cVel1 = bodies_[c.body1]->cvel;
        Vector3d cVel2 = Vector3d::Zero();
        Vector3d theta1 = bodies_[c.body1]->theta;
        Vector3d theta2 = Vector3d::Zero();

        Matrix3d IInertiaInverse = Matrix3d::Zero();
        Matrix3d JInertiaInverse = bodies_[c.body1]->getTemplate().getInertiaTensor().inverse();

        if(c.body2 != -1){
            Minv2 = Matrix3d::Identity() / (bodies_[c.body2]->density* bodies_[c.body2]->getTemplate().getVolume());
            c2 = bodies_[c.body2]->c;
            w2 = bodies_[c.body2]->w;
            cVel2 = bodies_[c.body2]->cvel;
            theta2 = bodies_[c.body2]->theta;
            IInertiaInverse = bodies_[c.body2]->getTemplate().getInertiaTensor().inverse();
        }

        Vector3d dPhi = VectorMath::crossProductMatrix(w2).transpose() * c.rotB2.transpose() * (c.rotB1 * vertHit + c1 - c2)
                + c.rotB2.transpose() * (c.rotB1 * VectorMath::crossProductMatrix(bodies_[c.body1]->w) * vertHit 
                + (cVel1) - (cVel2));

        Matrix3d A = c.rotB2.transpose() * VectorMath::crossProductMatrix(c.rotB1 * vertHit + c1 - c2) * VectorMath::TMatrix(-theta2);
        Matrix3d B = c.rotB2.transpose() * (-c.rotB1 * VectorMath::crossProductMatrix(vertHit) * VectorMath::TMatrix(theta1));
        Matrix3d C = -c.rotB2.transpose();
        Matrix3d D = c.rotB2.transpose();

        double relativeVel = c.derivDist.transpose() * dPhi;
        if(relativeVel > 0) continue;

        //dg calcs
        Vector3d dgThetaI = c.derivDist.transpose() * A;
        Vector3d dgThetaJ = c.derivDist.transpose() * B;
        Vector3d dgCI = c.derivDist.transpose() * C;
        Vector3d dgCJ = c.derivDist.transpose() * D;

        double alpha = -(1.0 + params_->CoR) * relativeVel;

        alpha /= c.derivDist.transpose() * (VectorMath::crossProductMatrix(c.rotB1 * vertHit + c1 - c2) * (dgThetaI.transpose() * VectorMath::TMatrix(theta2).inverse() * IInertiaInverse).transpose()
                    + c.rotB2.transpose() * c.rotB1 * VectorMath::crossProductMatrix((dgThetaJ.transpose() * VectorMath::TMatrix(theta1).inverse() * JInertiaInverse).transpose()) * vertHit
                    + Minv1 * dgCJ - Minv2 * dgCI);

        if(c.body2 != -1) {
            bodies_[c.body2]->cvel += Minv2 * (alpha * dgCI);
            bodies_[c.body2]->w = (bodies_[c.body2]->w.transpose() + (alpha * dgThetaI.transpose()) * VectorMath::TMatrix(theta2).inverse() * IInertiaInverse).transpose();
        }
        bodies_[c.body1]->cvel += Minv1 * (alpha * dgCJ);
        bodies_[c.body1]->w += ((alpha * dgThetaJ.transpose()) * VectorMath::TMatrix(theta1).inverse() * JInertiaInverse).transpose();

        free(it->second);
    }
}

}
