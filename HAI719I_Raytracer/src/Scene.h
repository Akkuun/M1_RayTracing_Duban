#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "Kd_tree.h"


#include <GL/glut.h>
#include <random>


enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace{};
    LightType type;

    Vec3 pos;
    float radius{};

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection {
    bool intersectionExists;
    unsigned int typeOfIntersectedObject{};
    unsigned int objectIndex{};
    float t;
    Vec3 intersection;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;

    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX) {}
};


class Scene {
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<Square> squares;
    std::vector<Light> lights;

public:


    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for (const auto &mesh: meshes) {
            mesh.draw();
        }
        for (const auto &sphere: spheres) {
            sphere.draw();
        }
        for (const auto &square: squares) {
            square.draw();
        }
    }

    //function that returns the normal of the object intersected by the ray
    static Vec3 getNormalOfObject(RaySceneIntersection &intersect) {
        if (intersect.typeOfIntersectedObject == 0) {
            return intersect.rayMeshIntersection.normal;
        } else if (intersect.typeOfIntersectedObject == 1) {
            return intersect.raySphereIntersection.normal;
        } else if (intersect.typeOfIntersectedObject == 2) {
            return intersect.raySquareIntersection.normal;
        }
        return {1, 1, 1};
    }

    RaySceneIntersection computeIntersection(Ray const &ray) {

        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        RaySceneIntersection result;
        RaySphereIntersection rayonSphere;
        RaySquareIntersection rayonSquare;
        RayTriangleIntersection rayonMesh;

        //Parcourir les objets de la scène ici les sphères
        for (unsigned int i = 0; i < spheres.size(); i++) {

            rayonSphere = spheres[i].intersect(ray);
            // On vérifie que le rayon a une intersection.
            // Si t < 0, alors l'origine du rayon est à l'intérieur de la sphère.
            // Si t >= 0, alors le rayon intersecte la sphère.
            // On compare les intersections pour obtenir l'objet le plus proche.

            //si on a bien une intersection avec notre rayon
            if (rayonSphere.intersectionExists && rayonSphere.t >= 0 && rayonSphere.t < result.t) {

                if (spheres[i].material.nTaille != 0) {
                    spheres[i].applyTexture(spheres[i].material, rayonSphere.intersection);
                }

                result.raySphereIntersection = rayonSphere;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 1;
                result.t = rayonSphere.t;
            }
        }

        //Parcourir les objets de la scène ici les carrés
        for (unsigned int i = 0; i < squares.size(); i++) {
            rayonSquare = squares[i].intersect(ray);
            if (rayonSquare.intersectionExists && rayonSquare.t >= 0 && rayonSquare.t < result.t) {

                if (squares[i].material.nTaille != 0) {
                    squares[i].applyTexture(rayonSquare.u, rayonSquare.v, squares[i].material);
                }
                if(squares[i].material.nTaille2 !=0){
                    squares[i].applyNormalMap(rayonSquare.u, rayonSquare.v, squares[i].material);
                }

                result.raySquareIntersection = rayonSquare;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 2;
                result.t = rayonSquare.t;
            }

        }

        //Parcourir les objets de la scène ici les meshs
        for (unsigned int i = 0; i < meshes.size(); i++) {
            rayonMesh = meshes[i].intersect(ray);
            if (rayonMesh.t >= 0 && rayonMesh.t < result.t && rayonMesh.intersectionExists) {


                if (meshes[i].material.nTaille != 0) {
                    meshes[i].applyTexture(rayonMesh.u, rayonMesh.v, meshes[i].material);
                }


                result.rayMeshIntersection = rayonMesh;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 0;
                result.t = rayonMesh.t;

            }
        }


        return result;
    }

    // 0 : mesh, 1 : sphere, 2 : square
    static Vec3 getIntersectionOfObject(RaySceneIntersection &intersect) {

        if (intersect.typeOfIntersectedObject == 0) {
            return intersect.rayMeshIntersection.intersection;
        } else if (intersect.typeOfIntersectedObject == 1) {
            return intersect.raySphereIntersection.intersection;
        } else if (intersect.typeOfIntersectedObject == 2) {
            return intersect.raySquareIntersection.intersection;
        }
        return {1, 1, 1};
    }


    // Retourne la couleur de l'objet intersecté

    // 0 : couleur du mesh
    // 1 : couleur de la sphère
    // 2 : couleur du carré
    Vec3 getColorOfObject(RaySceneIntersection intersect) {
        if (intersect.typeOfIntersectedObject == 0) {
            return meshes[intersect.objectIndex].material.diffuse_material;
        } else if (intersect.typeOfIntersectedObject == 1) {
            return spheres[intersect.objectIndex].material.diffuse_material;
        } else if (intersect.typeOfIntersectedObject == 2) {
            return squares[intersect.objectIndex].material.diffuse_material;
        }
        return {1, 1, 1};
    }

    Material getMaterial(RaySceneIntersection intersect) {
        if (intersect.typeOfIntersectedObject == 0) {
            return meshes[intersect.objectIndex].material;
        } else if (intersect.typeOfIntersectedObject == 1) {
            return spheres[intersect.objectIndex].material;
        } else if (intersect.typeOfIntersectedObject == 2) {
            return squares[intersect.objectIndex].material;
        }
        return {};

    }

    //function qui retourne si un point est dans l'ombre
    // pour cela on va vérifier si le point est dans l'ombre d'un objet
    // on calcul le vecteur de lumière réfléchie, si
    bool isInShadow(const Vec3 &point, const Light &light) {
        Vec3 lightDir = light.pos - point;
        float distanceToLight = lightDir.length();
        lightDir.normalize();

        Ray shadowRay(point + lightDir * 0.001f, lightDir); // Légèrement décalé pour éviter l'auto-intersection
        RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

        return shadowIntersection.intersectionExists && shadowIntersection.t < distanceToLight;
    }


    float softShadowCoeff(const Vec3 &point, const Light &light) {
        int numSamples = 10;
        int numHits = 0;
        for (int i = 0; i < numSamples; i++) {
            Vec3 lightDir = light.pos - point;
            lightDir.normalize();
            Vec3 randomOffset = Vec3::random(-0.5, 0.5);
            Vec3 newLightPos = light.pos + randomOffset;
            Vec3 newLightDir = newLightPos - point;
            float distanceToLight = newLightDir.length();
            newLightDir.normalize();

            Ray shadowRay(point + newLightDir * 0.001f,
                          newLightDir); // Légèrement décalé pour éviter l'auto-intersection
            RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

            if (shadowIntersection.intersectionExists && shadowIntersection.t < distanceToLight) {
                numHits++;
            }
        }
        return 1.0f - (float) numHits / (float) numSamples;

    }

    Vec3 PhongModel(RaySceneIntersection &intersect, Ray &ray) {


        Vec3 ambiantReflexion, diffuseReflexion, specularReflexion;

        Material material = getMaterial(intersect); // on récupère le materiaux de l'objet intersecté

        Vec3 color = Vec3(0, 0, 0); // résultat final à envoyé composé des 3 refléxions
        //donnée d'entrée pour le calcul du phong


        // P le point d'intersection sur la surface
        // N la normale du point d'intersection
        // V le vecteur vue à partir du point d'intersection
        // L le vecteur de lumière réfléchie
        Vec3 P = getIntersectionOfObject(intersect);
        Vec3 N = getNormalOfObject(intersect); // récupération de la normale de l'objet intersecté
        Vec3 V = ray.origin() - P;
        V.normalize(); // on normalise le vecteur V car on veut la direction de la vue


        // Réflexion ambiante : Ia
        //Ia = Isa  * Ka
        //Isa intensite de la lumière ambiante (donnée par le materiaux)
        //Ka coeffecient de la reflexion ambiante de l'object

        ambiantReflexion = Vec3::compProduct(lights[0].material, material.ambient_material);
        color += ambiantReflexion;

        // calculé pour chaqu lumière

        // Réflexion diffuse : Id
        //Id = Isd * Kd * max(0, N.L)
        //Isd intensite de la lumière diffuse (donnée par le materiaux)
        //Kd coeffecient de la reflexion diffuse de l'object

        // Réflexion spéculaire : Is
        //Is = Iss * Ks * max(0, R.V)^n
        //Iss intensite de la lumière spéculaire (donnée par le materiaux)
        //Ks coeffecient de la reflexion spéculaire de l'object
        //R vecteur de reflexion
        //V vecteur vue
        //n brillance de l'object



        for (int i = 0; i < lights.size(); i++) {

            //float softShadow = 1;
            //calcul du facteur pour l'ombre douce
            float softShadow = softShadowCoeff(P, lights[i]);

            Vec3 L = lights[0].pos - P;
            L.normalize();


            int nombreEchantillon = 10;
            float facteurOmbre = 0; //facteur d'ombre final issu des calculs échantilloné

            // fonction qui calcule le facteur d'ombre pour un point donné

            //si le point est dans l'ombre, on laisse en noir, on ne fait pas de calcul -> UTILISE pour l'ombre DURE
//            if (isInShadow(P, lights[i])) {
//                continue;
//            }
            //sinon on fait les calculs car le point n'est pas dans l'ombre

            // Réflexion diffuse : Id
            float LN = Vec3::dot(N, L);
            if (LN > 0) { // si le produit scalaire est positif, on fait le calcul
                diffuseReflexion = Vec3::compProduct(lights[i].material, material.diffuse_material) * LN;
                color += diffuseReflexion;
            }

            //Réflexion spéculaire : Is

            Vec3 R = 2 * Vec3::dot(N, L) * N - L;
            float RV = Vec3::dot(R, V);
            if (RV > 0) {
                specularReflexion =
                        Vec3::compProduct(lights[i].material, material.specular_material) * pow(RV, material.shininess);
                color += specularReflexion;
            }
            color *= softShadow;


        }


        return color;
    }

    static double reflectance(double cosine, double refraction_index) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - refraction_index) / (1 + refraction_index);
        r0 = r0 * r0;
        return r0 + (1 - r0) * std::pow((1 - cosine), 5);
    }

    static float random_float() {
        static std::uniform_real_distribution<float> distribution(0.0, 1.0);
        static std::mt19937 generator(static_cast<unsigned int>(time(nullptr)));
        return distribution(generator);
    }


    Vec3 rayTraceRecursive(Ray ray, int NRemainingBounces) {

        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color;

        if (!raySceneIntersection.intersectionExists) {
            //si notre rayon envoyé ne touche pas un objet, on affiche le ciel
            float a = 0.5 * (ray.direction()[1] + 1.0);
            return 0.1 * Vec3(1.0, 1.0, 1.0) + a * Vec3(0.5, 0.7, 1.0);
        }

        Material m = getMaterial(raySceneIntersection);

        // Traitement pour les matériaux miroir
        if (m.type == Material_Mirror && NRemainingBounces > 0) {
            Vec3 P = getIntersectionOfObject(raySceneIntersection);
            Vec3 N = getNormalOfObject(raySceneIntersection);
            Vec3 V = ray.origin() - P;
            V.normalize();

            Vec3 R = 2 * Vec3::dot(N, V) * N - V;
            Ray reflectedRay(P + R * 0.001f, R); // Légèrement décalé pour éviter l'auto-intersection
            color = rayTraceRecursive(reflectedRay, NRemainingBounces - 1); // Calculer la couleur réfléchie

            // Traitement pour les matériaux en verre
        }
            //refraction case for the glass material
        else if (m.type == Material_Glass && NRemainingBounces > 0) {
            Vec3 P = getIntersectionOfObject(raySceneIntersection);
            Vec3 N = getNormalOfObject(raySceneIntersection);
            Vec3 V = ray.direction();
            V.normalize();

            float cosTheta = Vec3::dot(V, N);
            float ni = 1.0f; // Index of refraction for air
            float nt = m.index_medium; // Index of refraction for the material
            Vec3 normal;

            if (cosTheta < 0) {
                // Entering the material
                normal = N;
                cosTheta = -cosTheta;
            } else {
                // Exiting the material
                ni = m.index_medium;
                nt = 1.0f;
                normal = -N;
            }

            float eta = ni / nt;
            float k = 1.0f - eta * eta * (1.0f - cosTheta * cosTheta);


            if ((eta * cosTheta) - 0.6f > 1.0 || reflectance(cosTheta, ni) > (double) random_float()) {
                // Total internal reflection
                Vec3 R = 2 * Vec3::dot(N, V) * N - V;
                Ray reflectedRay(P + R * 0.001f, R); // Slightly offset to avoid self-intersection
                color = rayTraceRecursive(reflectedRay, NRemainingBounces - 1); // Calculate the reflected color
            } else {


                Vec3 refraction = eta * V + (eta * cosTheta - sqrtf(k)) * normal;
                Ray refractedRay(P + refraction * 0.001f, refraction); // Slightly offset to avoid self-intersection
                color = rayTraceRecursive(refractedRay, NRemainingBounces - 1); // Calculate the refracted color
            }

        } else {
            // Calculer la couleur en utilisant le modèle de Phong
            color = PhongModel(raySceneIntersection, ray);
        }

        return color;
    }

// Fonction principale de l'algorithme de ray tracing

    Vec3 rayTrace(Ray const &rayStart) {
        //TODO appeler la fonction recursive
        Vec3 color;
        RaySceneIntersection intersect = computeIntersection(rayStart); // on calcul nos intersections (cercle,ect...)

        if (!intersect.intersectionExists) {
            //si notre rayon envoyé ne touche pas un objet, on affiche le ciel
            float a = 0.5f * (rayStart.direction()[1] + 1.0);
            return 0.1 * Vec3(1.0, 1.0, 1.0) + a * Vec3(0.5, 0.7, 1.0);
        } else {
            color = rayTraceRecursive(rayStart, 5);

            return color;
        }


    }

    void setup_cornel_box2() {

        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.quad = Square(Vec3(light.pos[0] - light.radius / 2,
                                     light.pos[1] - light.radius / 2,
                                     light.pos[2] - light.radius / 2),
                                Vec3(1, 0, 0),
                                Vec3(0, 0, 1),
                                light.radius,
                                light.radius);

            light.powerCorrection = 2.f;
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material = Material("img/planesTextures/vegetation.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);


//            s.material.diffuse_material = Vec3(1., 1., 1.);
//            s.material.specular_material = Vec3(1., 1., 1.);
//            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material = Material("img/planesTextures/paint.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);

//            s.material.diffuse_material = Vec3(1., 0., 0.);
//            s.material.specular_material = Vec3(1., 0., 0.);
//            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material = Material("img/planesTextures/galaxy.ppm", Vec3(0.0, 1.0, 0.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);

//            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
//            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
//            s.material.shininess = 16;
        }

        { //Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material = Material("img/planesTextures/bricks.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);
//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }
        { //Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();

            s.material = Material("img/planesTextures/marbre.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);
//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            //s.material.index_medium = 1.5;

            // s.material.index_medium = 0.64532232323;
            s.material.index_medium = 2.3; // marche pas mal avec 1.7
        }


        { //MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        //Mesh
        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            // bones et bear pas mal
            m.loadOFF("./src/mesh/SeaMonster.off");
            m.translate(Vec3(-2.5, 3, -2));
            m.rotate_y(0);
            m.rotate_x(0);
            m.scale(Vec3(0.4, 0.4, 0.4));
            m.build_arrays();
//            m.material = Material("img/planesTextures/galaxy.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
//                                  Vec3(0.0, 0.0, 0.0), 16);
            m.setHasTree();
            m.build_kdTree(); // Tree construction

//DECOMMENTEZ POUR AVOIR LE MESH DE BASE f
            m.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 16;
//DECOMMENTEZ POUR APPLIQUER UN MATERIAUX AU MESH
            //   m.material.type = Material_Mirror;
        }
        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            // bones et bear pas mal
            m.loadOFF("./src/mesh/dragon.off");
            m.translate(Vec3(0.15, -0.15, -0.1));
            m.rotate_y(0);
            m.rotate_x(0);
            m.scale(Vec3(6, 6, 6));
            m.build_arrays();
            m.material = Material("img/planesTextures/galaxy.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);
            m.setHasTree();
            m.build_kdTree(); // Tree construction

//DECOMMENTEZ POUR AVOIR LE MESH DE BASE f
            m.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.ambient_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 16;

        }

        { //Top Left Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            // s.m_center = Vec3(-1.5, 1.5, -1.5);
            s.m_center = Vec3(0, 0, -2);
            //s.m_radius = 0.5f;
            s.m_radius = 0.79f;
            s.build_arrays();

            s.material = Material("img/sphereTextures/s1.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);
//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE

//            s.material.diffuse_material = Vec3(1., 1., 1.);
//            s.material.specular_material = Vec3(1., 1., 1.);
//
//
//            s.material.shininess = 16;
//            s.material.transparency = 0.;
//            s.material.index_medium = 0.;
        }

        { //Top Right Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.5, 1.5, -1.5);
            s.m_radius = 0.5f;
            s.build_arrays();

            s.material = Material("img/sphereTextures/s6.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);

//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE
//            s.material.diffuse_material = Vec3(1., 1., 1.);
//            s.material.specular_material = Vec3(1., 1., 1.);
//            s.material.shininess = 16;
//            s.material.transparency = 0.;
//            s.material.index_medium = 0.;
        }

        { //Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material = Material("img/planesTextures/woodWall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }

        std::cout << "\n\n" << "\r" << "\033[35mPress R anytime to run the rendering process !" << "\033[0m" << std::endl;



    }

    void setup_scene_multiple_spheres() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        // Add a light source
        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 5.0, 0.0);
            light.radius = 2.5f;
            light.quad = Square(Vec3(light.pos[0] - light.radius / 2,
                                     light.pos[1] - light.radius / 2,
                                     light.pos[2] - light.radius / 2),
                                Vec3(1, 0, 0),
                                Vec3(0, 0, 1),
                                light.radius,
                                light.radius);
            light.powerCorrection = 2.f;
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        // Add the plane
        {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 5, 5.);
            s.translate(Vec3(-1, -1, -1.4));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();

//            s.material = Material("img/planesTextures/grass.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
//                                  Vec3(0.0, 0.0, 0.0), 16);


            s.material.diffuse_material = Vec3(0.8, 0.8, 0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
        }
        // Add small colored spheres
        int numSpheres = 20;
        float radius = 0.1f;
        for (int i = 0; i < numSpheres; ++i) {
            float x = static_cast<float>(rand()) / RAND_MAX * 8.0f - 4.0f;
            float z = static_cast<float>(rand()) / RAND_MAX * 8.0f - 4.0f;
            Vec3 color = Vec3(static_cast<float>(rand()) / RAND_MAX,
                              static_cast<float>(rand()) / RAND_MAX,
                              static_cast<float>(rand()) / RAND_MAX);
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(x, -1.3, z);
            s.m_radius = radius;
            s.build_arrays();
            s.material.diffuse_material = color;
            s.material.specular_material = color;
            s.material.shininess = 16;
        }
        // Add central refractive sphere
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-0.5, -0.9, 0.0);
            s.m_radius = 0.5f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 0., 0.);
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }
        // Add central reflective sphere
        {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0.5, -0.9, 0.0);
            s.m_radius = 0.5f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

    void setup_debug_refraction() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-1.0, 8., 2.0);
            light.radius = 1.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }
        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(-2., 2., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 0., 0.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }
        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(-2., -2., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(0., 1., 0.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }
        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(2., 2., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(0., 0., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }
        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(2., -2., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        }

        { //GLASS Sphere

            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(1, 1, 1);
            s.material.specular_material = Vec3(1, 1, 1);
            s.material.shininess = 16;

            s.material.transparency = 1.0;
            //s.material.index_medium = 1.5;

            // s.material.index_medium = 0.64532232323;
            s.material.index_medium = 2.3; // marche pas mal avec 1.7
        }
    }
    void setup_cornel_box_mesh() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.quad = Square(Vec3(light.pos[0] - light.radius / 2,
                                     light.pos[1] - light.radius / 2,
                                     light.pos[2] - light.radius / 2),
                                Vec3(1, 0, 0),
                                Vec3(0, 0, 1),
                                light.radius,
                                light.radius);

            light.powerCorrection = 2.f;
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material = Material("img/planesTextures/brickwall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");


//            s.material.diffuse_material = Vec3(1., 1., 1.);
//            s.material.specular_material = Vec3(1., 1., 1.);
//            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material = Material("img/planesTextures/brickwall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");

//            s.material.diffuse_material = Vec3(1., 0., 0.);
//            s.material.specular_material = Vec3(1., 0., 0.);
//            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material = Material("img/planesTextures/brickwall.ppm", Vec3(0.0, 1.0, 0.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");

//            s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
//            s.material.specular_material = Vec3(0.0, 1.0, 0.0);
//            s.material.shininess = 16;
        }

        { //Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material = Material("img/planesTextures/brickwall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");
//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }
        { //Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();

            s.material = Material("img/planesTextures/brickwall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");
//DECOMMENTEZ POUR AVOIR LES SPHERES DE BASE

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }






        //Mesh
        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            // bones et bear pas mal
            m.loadOFF("./src/mesh/pyramid.off");
            m.translate(Vec3(-2.5, 0, -2));
            m.rotate_y(0);
            m.rotate_x(0);
            m.scale(Vec3(0.8, 0.8, 0.8));
            m.build_arrays();
            m.material = Material("img/planesTextures/grass.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16,"img/normalMaps/brickwall_normal.ppm");
//            m.setHasTree();
//            m.build_kdTree(); // Tree construction

//DECOMMENTEZ POUR AVOIR LE MESH DE BASE f
//            m.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            m.material.shininess = 16;
//DECOMMENTEZ POUR APPLIQUER UN MATERIAUX AU MESH
            //   m.material.type = Material_Mirror;
        }
        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            // bones et bear pas mal
            m.loadOFF("./src/mesh/tripod.off");
            m.translate(Vec3(0, 0, -1));
            m.rotate_y(0);
            m.rotate_x(0);
            m.scale(Vec3(2, 2, 2));

            m.build_arrays();
            m.material = Material("img/planesTextures/marbre.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);
//            m.setHasTree();
//            m.build_kdTree(); // Tree construction

//DECOMMENTEZ POUR AVOIR LE MESH DE BASE f
//            m.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            m.material.shininess = 16;

        }





        { //Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material = Material("img/planesTextures/woodWall.ppm", Vec3(1.0, 1.0, 1.0), Vec3(1.0, 1.0, 1.0),
                                  Vec3(0.0, 0.0, 0.0), 16);

//            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
//            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
//            s.material.shininess = 16;
        }

        std::cout << "\n\n" << "\r" << "\033[35mPress R anytime to run the rendering process !" << "\033[0m" << std::endl;


    }


};


#endif