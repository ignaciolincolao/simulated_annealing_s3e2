#ifndef STRUCT_DIST
#define STRUCT_DIST


struct structDist {
    double distance;
    int id;
    bool operator<(const structDist& other) const {
        return distance < other.distance;
    }
};


#endif