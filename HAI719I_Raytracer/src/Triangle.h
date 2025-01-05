#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection {
    bool intersectionExists;
    float t;
    float w0, w1, w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
    float u, v;
};

class Triangle {
private:
    Vec3 m_normal;
    float area;
public:
    Triangle() {}
    unsigned int index;
    Triangle(Vec3 const &c0, Vec3 const &c1, Vec3 const &c2) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }

    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }

    void setC0(Vec3 const &c0) { m_c[0] = c0; }
    void setC1(Vec3 const &c1) { m_c[1] = c1; }
    void setC2(Vec3 const &c2) { m_c[2] = c2; }
    Vec3 const &normal() const { return m_normal; }


    //we need to know if the line is parallel to the triangle
    bool isParallelTo(Line const &L) const {
        // on fait le produit scalaire entre la direction du rayon et la normale du triangle
        // si ce produit scalaire est nul, alors le rayon est parallèle au triangle
        return Vec3::dot(L.direction(), m_normal) == 0; // si le produit scalaire est nul, alors les deux vecteurs sont orthogonaux donc parallèles
    }

    //on sait que le rayon n'est pas parallèle au triangle, on va maintenant vériifer si le triangle est devant le rayon
    bool isFrontFacing(Ray const &ray) const {
        // on fait le produit scalaire entre la direction du rayon et la normale du triangle
        // si ce produit scalaire est négatif, alors le triangle est devant le rayon
        return Vec3::dot(ray.direction(), m_normal) < 0; // si le produit scalaire est négatif, alors le triangle est devant le rayon
    }


    RayTriangleIntersection getIntersection(Ray const &ray) const {
        RayTriangleIntersection result;

        Vec3 a = m_c[0]; // les 3 sommets du triangle
        Vec3 b = m_c[1];
        Vec3 c = m_c[2];

        Vec3 origin = ray.origin();
        Vec3 direction = ray.direction();



        // 1) check that the ray is not parallel to the triangle:
        if (isParallelTo(Line(origin, direction))) {
            result.intersectionExists = false;
            return result;
        }

        // 2) check that the triangle is "in front of" the ray:
        if (!isFrontFacing(ray)) {
            result.intersectionExists = false;
            return result;
        }




        // 3) check that the intersection point is inside the triangle:


        float D = Vec3::dot(m_normal, a);
        float t = (D - Vec3::dot(m_normal, origin))/Vec3::dot(m_normal, direction);

        if (t < 0) {
            result.intersectionExists = false;
            return result;
        }
        Vec3 p = origin + t * direction;


        Vec3 v0 = b - a, v1 = c - a, v2 = p - a;
        float d00 = Vec3::dot(v0, v0);
        float d01 = Vec3::dot(v0, v1);
        float d11 = Vec3::dot(v1, v1);
        float d20 = Vec3::dot(v2, v0);
        float d21 = Vec3::dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;
        float v = (d11 * d20 - d01 * d21) / denom;
        float w = (d00 * d21 - d01 * d20) / denom;
        float u = 1.0f - v - w;

        if (v < 0 || w < 0 || u < 0) {
            result.intersectionExists = false;
            return result;
        }
        else{
            result.intersectionExists = true;
            result.t = t;
            result.w0 = u;
            result.w1 = v;
            result.w2 = w;
            result.intersection = p;
            result.normal = m_normal;
            result.u = u;
            result.v = v;
        }





        return result;
    }


    Vec3 m_c[3];

    Vec3 const & getFirstPoint() const {
        return m_c[0];
    }
    Vec3 const & getSecondPoint() const {
        return m_c[1];
    }
    Vec3 const & getThirdPoint()  const{
        return m_c[2];
    }

    Vec3 getMedianPoint() const {
        return (getFirstPoint() + getSecondPoint() + getThirdPoint()) / 3.0f;
    }
};

#endif
