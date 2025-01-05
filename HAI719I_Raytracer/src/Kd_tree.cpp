#include "Kd_tree.h"
#include "Mesh.h"
#include <limits>
#include <cmath>
#include <algorithm>

using namespace std;

Kd_tree::Kd_tree(Mesh &mesh) {
    // Initialisation de la racine
    left = nullptr;
    right = nullptr;
    depth = 0;

    Vec3 min = mesh.vertices[0].position;
    Vec3 max = mesh.vertices[0].position;

    for (const MeshVertex &v: mesh.vertices) {
        for (int p = 0; p < 3; p++) {
            if (v.position[p] < min[p])
                min[p] = v.position[p];
            if (v.position[p] > max[p])
                max[p] = v.position[p];
        }
    }

    englobingCube = EnglobingCube(min, max);

    triangles = mesh.trianglesObject;

    bestAxis = findBestAxis(englobingCube);

    maxDepth = static_cast<int>(log2(triangles.size()));
    generateKdTree();
}

void Kd_tree::generateKdTree() {
    if (triangles.size() <= 5 || depth >= maxDepth) {
        return; // Leaf node
    }


    // Sort triangles according to the chosen axis
    std::vector<Triangle> triangleSorted = sortTriangle(triangles, bestAxis);
    float medianValue = triangleSorted[triangleSorted.size() / 2].getMedianPoint()[bestAxis];

    // Split triangles with overlap
    std::vector<Triangle> leftTriangles, rightTriangles;
    for (const auto &triangle: triangleSorted) {
        if (triangle.getMedianPoint()[bestAxis] <= medianValue) {
            leftTriangles.push_back(triangle);
        }
        if (triangle.getMedianPoint()[bestAxis] >= medianValue) {
            rightTriangles.push_back(triangle);
        }
    }

    // Create child nodes
    left = new Kd_tree();
    createNode(*left, leftTriangles, (bestAxis + 1) % 3);

    right = new Kd_tree();
    createNode(*right, rightTriangles, (bestAxis + 1) % 3);
}

std::vector<Triangle> Kd_tree::sortTriangle(std::vector<Triangle> &triangles, int axis) {
    auto compareTriangles = [axis](const Triangle &a, const Triangle &b) {
        return a.getMedianPoint()[axis] < b.getMedianPoint()[axis];
    };

    std::sort(triangles.begin(), triangles.end(), compareTriangles);
    return triangles;
}

int Kd_tree::findBestAxis(EnglobingCube &cube) {
    Vec3 diag = cube.max - cube.min;
    if (diag[0] >= diag[1] && diag[0] >= diag[2]) return 0; // Axe X
    else if (diag[1] >= diag[0] && diag[1] >= diag[2]) return 1; // Axe Y
    else return 2; // Axe Z
}

void Kd_tree::setDepth(int depth) {
    this->depth = depth;
}

void Kd_tree::setMaxDepth(int maxDepth) {
    this->maxDepth = maxDepth;
}

void Kd_tree::setEnglobingCube(EnglobingCube englobingCube) {
    this->englobingCube = englobingCube;
}

void Kd_tree::setBestAxis(int bestAxis) {
    this->bestAxis = bestAxis;
}

void Kd_tree::drawTree(Vec3 color) const {
    this->englobingCube.draw(color);
}

void Kd_tree::drawEnglobingCube(Vec3 color) const {
    this->englobingCube.draw(color);
}

void Kd_tree::createNode(Kd_tree &node, std::vector<Triangle> tri, int axis) {
    node.setEnglobingCube(EnglobingCube(tri));
    node.triangles = tri;
    node.setDepth(depth + 1);
    node.setMaxDepth(maxDepth);
    node.setBestAxis(axis);
    node.generateKdTree();
}

bool Kd_tree::rayIntersectEnglobingCube(const Ray &ray, const EnglobingCube &cube) const {
    float invDirX = (ray.direction()[0] != 0) ? 1.0f / ray.direction()[0] : std::numeric_limits<float>::max();
    float tmin = (cube.min[0] - ray.origin()[0]) * invDirX;
    float tmax = (cube.max[0] - ray.origin()[0]) * invDirX;

    if (tmin > tmax) std::swap(tmin, tmax);

    float tymin = (cube.min[1] - ray.origin()[1]) / ray.direction()[1];
    float tymax = (cube.max[1] - ray.origin()[1]) / ray.direction()[1];

    if (tymin > tymax) std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax)) return false;

    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    float tzmin = (cube.min[2] - ray.origin()[2]) / ray.direction()[2];
    float tzmax = (cube.max[2] - ray.origin()[2]) / ray.direction()[2];

    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax)) return false;

    return true;
}

RayTriangleIntersection Kd_tree::nodeIntersection(const Ray &ray) const {
    RayTriangleIntersection result;
    result.t = std::numeric_limits<float>::max();

    if (rayIntersectEnglobingCube(ray, englobingCube)) {
        if (left == nullptr && right == nullptr) { // Nœud feuille
            for (const auto &triangle: triangles) {
                RayTriangleIntersection intersection = triangle.getIntersection(ray);
                if (intersection.intersectionExists && intersection.t < result.t) {
                    result = intersection;
                }
            }
        } else {
            // Récursion sur les enfants
            RayTriangleIntersection leftResult = left ? left->nodeIntersection(ray) : result;
            RayTriangleIntersection rightResult = right ? right->nodeIntersection(ray) : result;

            if (leftResult.intersectionExists && leftResult.t < result.t) {
                result = leftResult;
            }
            if (rightResult.intersectionExists && rightResult.t < result.t) {
                result = rightResult;
            }
        }
    }

    return result;
}