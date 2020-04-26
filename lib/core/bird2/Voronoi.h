#ifndef PSIM_CORE_BIRD2_VORONOI_H
#define PSIM_CORE_BIRD2_VORONOI_H
#include <Eigen/Core>
class VoronoiPoint
{
private:
    /* data */
public:
    Eigen::MatrixX4i T;
    Eigen::Vector3d center;
    VoronoiPoint(Eigen::Vector3d center, Eigen::MatrixX4i T){
        this->T = T;
        this->center = center;
    }
};

#endif