#ifndef KD_TREE_H
#define KD_TREE_H


#include "EnglobingCube.h"
#include <vector>

class Mesh;

using namespace std;

class Kd_tree {
public:
    Kd_tree() : left(nullptr), right(nullptr) {}
    Kd_tree( Mesh &mesh);

    void drawTree(Vec3 color) const;
    void drawEnglobingCube(Vec3 color) const;
    RayTriangleIntersection nodeIntersection(Ray const &ray) const;
private:
    Kd_tree* left;
    Kd_tree* right;
    EnglobingCube englobingCube;
    std::vector<Triangle> triangles;
    int depth;
    int maxDepth;
    int bestAxis;

    void generateKdTree();
    std::vector<Triangle> sortTriangle(std::vector<Triangle> &triangles, int axe);
    int findBestAxis(EnglobingCube &cube);
    void createNode(Kd_tree &node, std::vector<Triangle> tri, int axis);

    void setDepth(int depth);
    void setMaxDepth(int maxDepth);
    void setEnglobingCube(EnglobingCube englobingCube);
    void setBestAxis(int bestAxis);




    bool rayIntersectEnglobingCube(Ray const &ray,const EnglobingCube &cube) const;
};


#endif // KD_TREE_H