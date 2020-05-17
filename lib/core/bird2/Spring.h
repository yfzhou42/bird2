#ifndef PSIM_CORE_BIRD2_SPRING_H
#define PSIM_CORE_BIRD2_SPRING_H

class Spring
{
private:
public:
    // Index into Voronoi Point array
    int p1,p2;
    double length;
    Spring(int p1, int p2, double length){
        this->p1 = p1;
        this->p2 = p2;
        this->length = length;
    }
};

#endif