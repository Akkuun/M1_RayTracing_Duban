#ifndef SQUARE_H
#define SQUARE_H

#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>
#include <algorithm>

struct RaySquareIntersection {
    bool intersectionExists;
    float t;
    float u, v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    float width;
    float height;

    Square() : Mesh() {}

    Square(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
           float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void
    setQuad(Vec3 const &bottomLeft, Vec3 const &rightVector, Vec3 const &upVector, float width = 1., float height = 1.,
            float uMin = 0.f, float uMax = 1.f, float vMin = 0.f, float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector, upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector * width;
        m_up_vector = m_up_vector * height;
        this->width = width;
        this->height = height;
        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;
        vertices[0].u = uMin;
        vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;
        vertices[1].u = uMax;
        vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;
        vertices[2].u = uMax;
        vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;
        vertices[3].u = uMin;
        vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    void scale(Vec3 const &scale) {
        Mat3 scale_matrix(scale[0], 0., 0.,
                          0., scale[1], 0.,
                          0., 0., scale[2]); //Matrice de transformation de mise à l'échelle
        apply_transformation_matrix(scale_matrix);

        width *= scale[0];
        height *= scale[1];
    }

    //calcule l'équation de l'intersection entre le rayon et le plan (le côté du carré)
    void
    ComputeEquation(Vec3 &origin, Vec3 &direction, Vec3 &normale, RaySquareIntersection &raySquareIntersection,
                    Vec3 &bottomLeft, Vec3 &rightVector, Vec3 &upVector) {
        //t => distance entre le point d'origine du rayon et le point d'intersection
        float t = 0.0;

        //equation de t : t = (D - o.n) / d.n

        // D <> distance entre le point d'origine du rayon et le point d'intersection

        // equation du plan x.n - D = 0
        // D = x.n


        float D = Vec3::dot(bottomLeft,
                            normale); //  D = x.n   point entre le coin inférieur gauche du carré et la normale du plan

        float o_n = Vec3::dot(origin, normale);
        float d_n = Vec3::dot(direction, normale);

        //IGNORE COLLISIONS DOS (parralele 0 , >0 de dos)
        if (d_n >= 0) {
            raySquareIntersection.intersectionExists = false;
            return;
        }

        t = (D - o_n) / d_n;

        //4 cas -> infini rayon parrallele au plan donc pas d'intersection
        //t < 0 -> intersection derrière le point d'origine du rayon( derrière la caméra)
        // t > intersection devant la caméra

        if (t < 0) {
            //  std::cout << "Pas d'intersection ,  t = " << t << std::endl;
            raySquareIntersection.intersectionExists = false;
            return;
        } else {
            // on a une intersection avec le plan , maintenant on doit faire en sorte de s'assurer d'être dans le carré car sinon on va juste afficher tout le plan

            // on doit vérifier si le point d'intersection est dans le carré, on va faire une projection du point d'intersection sur le plan du carré et vérifier si le point est dans le carré

            // on a le plan du carré, on a le plan d'intersection,

            Vec3 intersection = origin + t *
                                         direction; // (o + t * d) représente le vecteur  d'intersection qui va du point d'origine du rayon au plan d'intersection

            Vec3 v = intersection -
                     bottomLeft; // on fait -bottomLeft pour avoir le vecteur qui va du point d'intersection au point bottomLeft, v est le vecteur qui va du point d'intersection du plan au point bottomLeft


            float dot1 = (Vec3::dot(v,
                                    rightVector)) /
                         rightVector.length(); // projection entre le point d'intersection et le vecteur rightVector qui est le vecteur qui va vers la droite pour dessiner le carré
            float dot2 = (Vec3::dot(v,
                                    upVector)) /
                         upVector.length(); //projection entre le point d'intersection et le vecteur upVector qui est le vecteur qui va vers le haut pour dessiner le carré

            if (dot1 < 0 || dot1 > rightVector.length() || dot2 < 0 || dot2 > upVector.length()) {
                raySquareIntersection.intersectionExists = false;
                return;
            } else {
                raySquareIntersection.intersectionExists = true;
                raySquareIntersection.t = t;
                raySquareIntersection.intersection = intersection;
                raySquareIntersection.u = dot1 / rightVector.length();
                raySquareIntersection.v = dot2 / upVector.length();
                //insertion de la normale du point d'intersection


                raySquareIntersection.normal = normale;

                raySquareIntersection.normal.normalize();

                if (raySquareIntersection.intersectionExists && material.nTaille2 != 0) {
                    Vec3 tangent = rightVector;
                    Vec3 bitangent = upVector;
                    tangent.normalize();
                    bitangent.normalize();
                    applyNormalMap(raySquareIntersection.u, raySquareIntersection.v, material, raySquareIntersection.normal, tangent, bitangent);
                }


            }


        }


    }

    //fonction qui calcule l'intersection entre un rayon et un carré
    RaySquareIntersection intersect(const Ray &ray) {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;

        Vec3 m_bottom_left = vertices[0].position;
        Vec3 m_right_vector = vertices[1].position -
                              vertices[0].position; // vecteur qui va vers la droite pour dessiner le carré (rightVector) en faisant x2 - x1 (x1 étant le point en bas à gauche du carré et x2 le point en bas à droite du carré)
        Vec3 m_up_vector = vertices[3].position -
                           vertices[0].position; // vecteur qui va vers le haut pour dessiner le carré (upVector) en faisant  x4 - x1 (x1 étant le point en bas à gauche du carré et x4 le point en haut à gauche du carré)

        Vec3 m_normal = Vec3::cross(m_right_vector, m_up_vector) / width / height;


        Vec3 origin = ray.origin();
        Vec3 direction = ray.direction();
        ComputeEquation(origin, direction, m_normal, intersection, m_bottom_left, m_right_vector, m_up_vector);

        return intersection;
    }

    //clamp function
    template<typename T>
    constexpr const T &clamp(const T &value, const T &min, const T &max) {
        return (value < min) ? min : (value > max) ? max : value;
    }


    //the idea of this function is to apply the texture to the square
    //we compute the coordinate u and v of the intersection point and change the Diffuse material of the square  of the same coordinate

    void applyNormalMap(float &u, float &v, Material &material, Vec3 &normal, Vec3 &tangent, Vec3 &bitangent) {
        // Récupérer la normale de la normal map
        Vec3 normalFromMap = getNormalUV(u, v, material);

        // Transformer la normale de l'espace tangent à l'espace monde
        Vec3 worldNormal = applyTangentToWorld(normalFromMap, tangent, bitangent, normal);

        // Normaliser la normale résultante
        normal = worldNormal;
        normal.normalize();
    }

    void applyNormalMap(float &u, float &v, Material &material) {
        // Obtenir la normale de la texture UV
        Vec3 normalValue = getNormalUV(u, v, material);
        normalValue.normalize();
        // Changer la normale de la surface


    }

    Vec3 applyTangentToWorld(Vec3 tangentNormal, Vec3 tangent, Vec3 bitangent, Vec3 normal) {
        return tangent * tangentNormal[0] +
               bitangent * tangentNormal[1] +
               normal * tangentNormal[2];
    }



    //get the Vec3 texture of the UV coordinate passed in parameter (local coordinates to index of the pixel in the texture in the array)
    Vec3 getTextureUV(float &u, float &v, Material &texture) {
        int nH = texture.nH; // height of the texture
        int nW = texture.nW; // width of the texture
        // Calculate the index of the pixel in the texture
        // Clamp u and v to the range [0, 1]
        u = clamp(u, 0.f, 1.f);
        v = clamp(v, 0.f, 1.f);


        // Calculate the index of the pixel in the texture
        float i = u * (float) nW;
        float j = v * (float) nH;

        int index = (static_cast<int>(j) * nW + static_cast<int>(i)) *
                    3; // index of the pixel in the texture based on rows * Ncols + cols

        if (index >= texture.nTaille * 3) {
            index = texture.nTaille * 3 - 3;
        }



        // Return the Vec3 of the pixel
        return Vec3(
                static_cast<float>(texture.ArrayImageForTexture[index]) / 255.0f,
                static_cast<float>(texture.ArrayImageForTexture[index + 1]) / 255.0f,
                static_cast<float>(texture.ArrayImageForTexture[index + 2]) / 255.0f
        );

    }

    Vec3 getNormalUV(float &u, float &v, Material &normal_map) {
        int nH = normal_map.nH; // height of the texture
        int nW = normal_map.nW; // width of the texture
        // Calculate the index of the pixel in the texture
        // Clamp u and v to the range [0, 1]
        u = clamp(u, 0.f, 1.f);
        v = clamp(v, 0.f, 1.f);


        // Calculate the index of the pixel in the texture
        float i = u * (float) nW;
        float j = v * (float) nH;

        int index = (static_cast<int>(j) * nW + static_cast<int>(i)) *
                    3; // index of the pixel in the texture based on rows * Ncols + cols

        if (index >= normal_map.nTaille2 * 3) {
            index = normal_map.nTaille2 * 3 - 3;
        }



        // Return the Vec3 of the pixel
        return Vec3(
                static_cast<float>(normal_map.ArrayImageForNormalMap[index]) / 255.0f,
                static_cast<float>(normal_map.ArrayImageForNormalMap[index + 1]) / 255.0f,
                static_cast<float>(normal_map.ArrayImageForNormalMap[index + 2]) / 255.0f
        );

    }


};

#endif // SQUARE_H