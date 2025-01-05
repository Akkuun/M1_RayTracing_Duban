#ifndef Sphere_H
#define Sphere_H

#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection {
    bool intersectionExists;
    float t;
    float theta, phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean(Vec3 ThetaPhiR) {
    return ThetaPhiR[2] *
           Vec3(cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]), sin(ThetaPhiR[1]));
}

static
Vec3 SphericalCoordinatesToEuclidean(float theta, float phi) {
    return {cos(theta) * cos(phi), sin(theta) * cos(phi), sin(phi)};
}

static
Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz) {
    float R = xyz.length();
    float phi = asin(xyz[2] / R);
    float theta = atan2(xyz[1], xyz[0]);
    return {theta, phi, R};
}


class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}

    Sphere(Vec3 c, float r) : Mesh(), m_center(c), m_radius(r) {}

    void build_arrays() {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt) {
            float u = (float) (thetaIt) / (float) (nTheta - 1);
            float theta = u * 2 * M_PI;
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float) (phiIt) / (float) (nPhi - 1);
                float phi = -M_PI / 2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normalsArray[3 * vertexIndex + 0] = xyz[0];
                normalsArray[3 * vertexIndex + 1] = xyz[1];
                normalsArray[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt) {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }
    }

    // equation de l'interaction rayon-shpere : t^2 * d*d + 2 *t *d *(o-c) + ||o-c||^2 - r^2 =0  avec  c: centre de la sphère , r rayon de la shpere


    //calcul la formule de l'intersection entre le rayon projeté et la sphère
    static void ComputeEquation(Vec3 &origin, Vec3 &cercleOrigin, float &sphereRadius,
                         Vec3 &d_direction, RaySphereIntersection &raySphereIntersection) {

        //t représentant la distance entre l'origine du rayon et le point d'intersection
        float t = 0;

        // équation de la forme ax^2 + bx + c = 0 avec a = d * d, b = 2 * d * (o - c), c = ||o - c||^2 - r^2

        //o - c
        Vec3 oc = origin - cercleOrigin;
        // a = d * d * t^2
        float a = Vec3::dot(d_direction, d_direction);

        // b = 2 * d * (o - c)
        float b = 2 * Vec3::dot(d_direction, oc);
        // c = ||o - c||^2 - r^2

        float c = oc.squareLength() - (sphereRadius * sphereRadius);

        //delta = b^2 - 4ac

        float delta = (b * b) - (4 * (a * c));

        if (delta == 0) { // si delta = 0 alors il y a une seule solution
            t = -b / (2 * a);


        } else if (delta > 0) { // si delta > 0 alors il y a deux solutions
            double t1 = (-b - sqrt(delta)) / (2 * a);
            double t2 = (-b + sqrt(delta)) / (2 * a);
            t = std::min(t1,
                         t2); // on prend la plus petite valeur de t car le cercle va être touché 2 fois par le rayon
            // 1 fois au plus proche et 1 fois au plus loin donc on prend le plus proche pour avoir la bonne intersection

        } else if (delta < 0) { // il n'y a pas de solution
            raySphereIntersection.intersectionExists = false;
            return;
        }

        raySphereIntersection.t = t; // on stocke la valeur de t
        raySphereIntersection.intersectionExists = true; // il y a une intersection
        raySphereIntersection.intersection = origin + t * d_direction; // on calcule le point d'intersection
        raySphereIntersection.normal = (raySphereIntersection.intersection -
                                        cercleOrigin); // on calcule la normale au point d'intersection
        raySphereIntersection.normal.normalize(); // on normalise la normale

    }

    //calcul l'intersection  entre un rayon et la sphère sphere
    RaySphereIntersection intersect(const Ray &ray) const {
        //on retourne l'intersection la plus proche
        RaySphereIntersection intersection;


        Vec3 origin = ray.origin(); // origin

        Vec3 cercleOrigin = this->m_center;

        float sphereRadius = this->m_radius;

        Vec3 d_direction = ray.direction();

        //calcul de l'intersection entre la shpere et le rayon
        ComputeEquation(origin, cercleOrigin, sphereRadius, d_direction, intersection);

        //calcul de la normale du point d'intersection
        //normal donnée par : N = (P - C) / r




        return intersection;
    }


    template<typename T>
    constexpr const T &clamp(const T &value, const T &min, const T &max) {
        return (value < min) ? min : (value > max) ? max : value;
    }

    //get the uv coordinates of the intersection into the index pixel of the texture
    Vec3 getUVCoordinates(const Vec3 &point, Material &texture) {
        Vec3 normalizedPoint = point -m_center;

        auto theta = acos(-normalizedPoint[1]);
        auto phi = atan2(-normalizedPoint[2], -normalizedPoint[0]) + M_PI;
        float u = 1.0 -phi / (2 * M_PI);
       float v = 1.0 -theta / M_PI;


        // i & j coordinates in the texture
        float i = u * (float)texture.nW;
        float j = v * (float)texture.nH;

        //get the index of the pixel in the texture (nTaille *3)
        int index = (static_cast<int>(j) * texture.nW + static_cast<int>(i)) * 3; // index of the pixel in the texture based on rows * Ncols + cols

        //if out of bounds we clamp the index
        if(index>=texture.nTaille*3){
            index = texture.nTaille*3 - 3;
        }

        return Vec3(
                static_cast<float>(texture.ArrayImageForTexture[index]) / 255.0f,
                static_cast<float>(texture.ArrayImageForTexture[index + 1]) / 255.0f,
                static_cast<float>(texture.ArrayImageForTexture[index + 2]) / 255.0f
        );



    }

    //apply texture with u and v of the intersection
    void applyTexture(Material &texture, Vec3 &point) {

        //get u v coordinates on the sphere

        Vec3 textureValue = getUVCoordinates(point, texture);
        texture.diffuse_material  = textureValue;



    }

};

#endif
