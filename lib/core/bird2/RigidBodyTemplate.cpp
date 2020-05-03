#include "RigidBodyTemplate.h"
#include "Distance.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <time.h>

#define VORONOI_POINTS 3

namespace bird2 {

using namespace std;
using namespace Eigen;

RigidBodyTemplate::RigidBodyTemplate(Eigen::Ref<Eigen::MatrixX3d> V,
                                     Eigen::Ref<Eigen::MatrixX3i> F,
                                     double scale) : volume_(0)
{
    srand(time(NULL));
    inertiaTensor_.setZero();
    Eigen::MatrixXd mV = V * scale;
    Eigen::MatrixXi mF = F;

    igl::copyleft::tetgen::tetrahedralize(mV, mF, "Qpq1.414a0.01", this->V, this->T, this->F);
    computeFaces();
    initialize();
    generateVoronoiPoints();
    std::cout << "Springs Count: " << springs.size() << std::endl;
    for (Spring s : springs) {
        std::cout << "Vor Connection: " << s.p1 << "," << s.p2 << std::endl;
    }
}

RigidBodyTemplate::RigidBodyTemplate(const Eigen::MatrixX3d& verts, const Eigen::MatrixX4i& tets)
    : volume_(0)
{
    V = verts;
    T = tets;
    computeFaces();    
    Eigen::MatrixXd mV = V;
    Eigen::MatrixXi mF = F;
    igl::copyleft::tetgen::tetrahedralize(V, mF, "Qpq1.414a0.01", this->V, this->T, this->F);
    computeFaces();
    initialize();
    generateVoronoiPoints();
}

Vector3d getTetCenter(const MatrixX3d& V, const MatrixX4i& T, int tet){
    return (V.row(T.row(tet)(0)) + V.row(T.row(tet)(1)) + V.row(T.row(tet)(2)) + V.row(T.row(tet)(3)))/4.0;
}

void RigidBodyTemplate::generateVoronoiPoints(){
    vector<Eigen::Vector3d> voronoiCenters;
    set<int> includedTets;
    for (int i = 0; i < VORONOI_POINTS; i++) {
        int tet = rand() % T.rows();
        if(includedTets.find(tet) != includedTets.end()) continue;
        voronoiCenters.push_back(getTetCenter(V, T, tet));
        includedTets.insert(tet);
    }
    
    tetToVoronoi.resize(T.rows());

    vector<vector<Vector4i>> voronoiTets;
    voronoiTets.resize(voronoiCenters.size());
    for (int i = 0; i < T.rows(); i++) {
        double minDist = 100000.0;
        int minIndex = -1;
        Vector3d tetCenter = getTetCenter(V, T, i);
        for(int j = 0; j < voronoiCenters.size(); j++){
            double curDist = (tetCenter-voronoiCenters[j]).squaredNorm();
            if ((tetCenter-voronoiCenters[j]).squaredNorm() < minDist){
                minDist = curDist;
                minIndex = j;
            }
        }
        if(minIndex == -1) std::cout << "THIS SHOULD NEVER EVER HAPPEN" << std::endl;
        voronoiTets[minIndex].push_back(T.row(i));
        //Save the voronoi->tet reference for lookup later.
        tetToVoronoi[i] = minIndex;
    }

    for(int i = 0; i< voronoiCenters.size(); i++){
        MatrixX4i tetMatrix;
        tetMatrix.resize(voronoiTets[i].size(), 4);
        for(int j = 0; j < voronoiTets[i].size(); j++){
            tetMatrix.row(j) = voronoiTets[i][j];
        }
        voronois.push_back(VoronoiPoint(voronoiCenters[i], tetMatrix));
        //std::cout << "Voronoi " << i << " Tet Count: " << voronois[i].T.rows() << std::endl;
    }
    generateSprings();
    generateVoronoiReferences();
}

void RigidBodyTemplate::generateVoronoiReferences(){
    vertToVoronoi.resize(V.rows());
    for (int i = 0; i < voronois.size(); i++) {
        for (int j = 0; j < voronois[i].T.rows(); j++) {
            for (int k = 0; k < 4; k++) {
                vertToVoronoi[voronois[i].T.row(j)[k]] = i;
            }
        }
    }
    
}

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

void RigidBodyTemplate::generateSprings() {
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
        for (int t1 = 0; t1 < T1.size(); t1++) {
            //for each tet which contains 4 faces
            set<set<int>> faces;
            makeFaces(T1.row(t1), faces);
            tetSet1.insert(faces.begin(), faces.end());
        }
        for (int j = i + 1; j < voronois.size(); j++) {
            MatrixX4i T2 = voronois[i].T;
            bool connected = false;
            for (int t2 = 0; t2 < T2.size() && !connected; t2++) {
                set<set<int>> faces;
                makeFaces(T2.row(t2), faces);
                for(set<int> f : faces){
                    if(tetSet1.find(f) != tetSet1.end()) {
                        springs.push_back(Spring(i, j));
                        voronois[i].springs.push_back(&springs.back());
                        voronois[j].springs.push_back(&springs.back());
                        connected = true;
                        break;
                    }
                }
            }
        }
    }
}

RigidBodyTemplate::~RigidBodyTemplate()
{
}

void RigidBodyTemplate::initialize()
{
    volume_ = computeVolume();
    com_ = computeCenterOfMass();
    for (int i = 0; i < V.rows(); i++)
        V.row(i) -= com_;
    inertiaTensor_ = computeInertiaTensor();

    //std::cerr << "computeDistances, F rows: " << F.rows() << std::endl;
    distances_ = computeDistances();
    //std::cerr << "computeDistances Done" << std::endl;
}

void RigidBodyTemplate::computeFaces()
{
    struct triple {
        triple(int aa, int bb, int cc)
            : a(aa)
            , b(bb)
            , c(cc)
        {
            if (a < b)
                std::swap(a, b);
            if (a < c)
                std::swap(a, c);
            if (b < c)
                std::swap(b, c);
        }

        int a, b, c;
        bool operator<(const triple& other) const
        {
            if (a < other.a)
                return true;
            else if (a > other.a)
                return false;
            if (b < other.b)
                return true;
            else if (b > other.b)
                return false;
            return c < other.c;
        }
    };

    int ntets = (int)T.rows();
    MatrixX3i allfaces(4 * ntets, 3);
    Matrix<int, 4, 3> faceidx;
    faceidx << 0, 1, 3,
               3, 1, 2,
               3, 2, 0,
               0, 2, 1;

    for (int i = 0; i < ntets; i++) {
        Vector4i tet = T.row(i);
        for (int face = 0; face < 4; face++) {
            for (int k = 0; k < 3; k++)
                allfaces(4 * i + face, k) = tet[faceidx(face, k)];
        }
    }

    map<triple, vector<int>> faces;
    for (int i = 0; i < 4 * ntets; i++) {
        triple t(allfaces(i, 0), allfaces(i, 1), allfaces(i, 2));
        faces[t].push_back(i / 4);
    }

    int nfaces = 0;
    for (map<triple, vector<int>>::iterator it = faces.begin(); it != faces.end(); ++it)
        if (it->second.size() == 1)
            nfaces++;

    F.resize(nfaces, 3);
    int idx = 0;

    for (int i = 0; i < 4 * ntets; i++) {
        triple t(allfaces(i, 0), allfaces(i, 1), allfaces(i, 2));
        if (faces[t].size() == 1) {
            F.row(idx) = allfaces.row(i);
            idx++;
        }
    }
}

double RigidBodyTemplate::computeVolume()
{
    double volume = 0;
    for (int i = 0; i < F.rows(); i++) {
        Vector3d pts[3];
        Vector3d centroid(0, 0, 0);
        for (int j = 0; j < 3; j++) {
            pts[j] = V.row(F(i, j));
            centroid += pts[j];
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        centroid /= 3.0;
        volume += centroid.dot(normal) * area / 3.0;
    }
    return volume;
}

Vector3d RigidBodyTemplate::computeCenterOfMass()
{
    Vector3d cm(0, 0, 0);
    for (int i = 0; i < F.rows(); i++) {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++) {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();

        Vector3d term(0, 0, 0);
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = k; l < 3; l++)
                    term[j] += pts[k][j] * pts[l][j];
            }
            term[j] *= area * normal[j] / 12.0;
        }

        cm += term;
    }

    return cm / volume_;
}

Eigen::Matrix3d RigidBodyTemplate::computeInertiaTensor()
{
    Eigen::Matrix3d inertiaTensor;
    inertiaTensor.setZero();

    Vector3d quads(0, 0, 0);
    Vector3d mixed(0, 0, 0);
    for (int i = 0; i < F.rows(); i++)
    {
        Vector3d pts[3];
        for (int j = 0; j < 3; j++)
        {
            pts[j] = V.row(F(i, j));
        }
        Vector3d normal = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        double area = 0.5 * normal.norm();
        normal /= normal.norm();


        Vector3d term(0, 0, 0);
        Vector3d mixterm(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = k; l < 3; l++)
                {
                    for (int m = l; m < 3; m++)
                        term[j] += pts[k][j] * pts[l][j] * pts[m][j];
                }
            }
            term[j] *= area*normal[j] / 30.0;
        }
        double mix = 0;
        for (int j = 0; j < 3; j++)
        {
            mix += 6.0*pts[j][0] * pts[j][1] * pts[j][2];
            for (int k = 0; k < 3; k++)
            {
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 1) % 3][(k + 2) % 3];
                mix += 2.0*pts[j][k] * pts[j][(k + 1) % 3] * pts[(j + 2) % 3][(k + 2) % 3];
            }
            mix += pts[j][0] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][2];
            mix += pts[j][2] * pts[(j + 1) % 3][1] * pts[(j + 2) % 3][0];
        }
        for (int j = 0; j < 3; j++)
            mixterm[j] = mix*area*normal[j] / 60.0;

        quads += term;
        mixed += mixterm;
    }

    inertiaTensor << quads[1] + quads[2], -mixed[2], -mixed[1],
        -mixed[2], quads[0] + quads[2], -mixed[0],
        -mixed[1], -mixed[0], quads[0] + quads[1];

    return inertiaTensor;
}

Eigen::VectorXd RigidBodyTemplate::computeDistances()
{
    int nverts = (int)V.rows();
    int nfaces = (int)F.rows();
    Eigen::VectorXd distances;
    distances.resize(nverts);
    for (int i = 0; i < nverts; i++) {
        double dist = numeric_limits<double>::infinity();
        for (int j = 0; j < nfaces; j++) {
            double dummy;
            if (Distance::vertexPlaneDistanceLessThan(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dist)) {
                Vector3d distvec = Distance::vertexFaceDistance(V.row(i), V.row(F(j, 0)), V.row(F(j, 1)), V.row(F(j, 2)), dummy, dummy, dummy);
                dist = min(dist, distvec.norm());
            }
        }
        distances(i) = dist;
    }
    return distances;
}


double tripleProduct(Vector3d a, Vector3d b, Vector3d c)
{
    return a.dot(b.cross(c));
}

Vector4d barycentricTet(Vector3d a, Vector3d b, Vector3d c, Vector3d d, Vector3d p)
{
    Vector3d vap = p - a;
    Vector3d vbp = p - b;

    Vector3d vab = b - a;
    Vector3d vac = c - a;
    Vector3d vad = d - a;

    Vector3d vbc = c - b;
    Vector3d vbd = d - b;

    double va6 = tripleProduct(vbp, vbd, vbc);
    double vb6 = tripleProduct(vap, vac, vad);
    double vc6 = tripleProduct(vap, vad, vab);
    double vd6 = tripleProduct(vap, vab, vac);
    double v6 = 1.0 / tripleProduct(vab, vac, vad);
    return Vector4d(va6*v6, vb6*v6, vc6*v6, vd6*v6);
}

double RigidBodyTemplate::distance(Vector3d p, int tet) const
{
    // TODO: Compute distance from point to object boundary
    if (tet == -1)
        return p(1) + 1.0;
    if (tet < -1)
        return 0.0;
    Vector4i tetraIndices = T.row(tet);
    Vector3d a = V.row(tetraIndices[0]);
    Vector3d b = V.row(tetraIndices[1]);
    Vector3d c = V.row(tetraIndices[2]);
    Vector3d d = V.row(tetraIndices[3]);

    Vector4d weights = barycentricTet(a, b, c, d, p);
    //std::cout << "distance to origin " << p.norm() << "distance to boundary " << 1.0 - p.norm()<< endl;

    //std::cout << "DIST TORET: " << distances_(tetraIndices[0]) * weights[0] + distances_(tetraIndices[1]) * weights[1] + distances_(tetraIndices[2]) * weights[2] + distances_(tetraIndices[3]) * weights[3] << std::endl;
    return -(distances_(tetraIndices[0]) * weights[0] + distances_(tetraIndices[1]) * weights[1] + distances_(tetraIndices[2]) * weights[2] + distances_(tetraIndices[3]) * weights[3]);
}

Vector3d RigidBodyTemplate::Ddistance(int tet) const
{
    // TODO: Compute derivative of distance from point to boundary
    Vector4i tetraIndices = T.row(tet);
    Vector3d a = V.row(tetraIndices[0]);
    Vector3d b = V.row(tetraIndices[1]);
    Vector3d c = V.row(tetraIndices[2]);
    Vector3d d = V.row(tetraIndices[3]);
    Matrix3d M;
    M.row(0) = (a - d).transpose();
    M.row(1) = (b - d).transpose();
    M.row(2) = (c - d).transpose();
    
    //Matrix3d M((a - d).transpose(), (b - d).transpose(), (c - d).transpose());
    double D0 = -distances_(tetraIndices[0]);
    double D1 = -distances_(tetraIndices[1]);
    double D2 = -distances_(tetraIndices[2]);
    double D3 = -distances_(tetraIndices[3]);

    Vector3d D(D0 - D3, D1 - D3, D2 - D3);

    return D.transpose() * M.inverse();
}

}