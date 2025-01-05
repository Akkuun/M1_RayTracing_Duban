#ifndef MATERIAL_H
#define MATERIAL_H

#include "image_utils.h"
#include <cmath>
#include <string>
#include <GL/glut.h>
#include <cstring>
#include <unordered_map>

enum MaterialType {
    Material_Diffuse_Blinn_Phong,
    Material_Glass,
    Material_Mirror,

};



struct Material {
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;
    double shininess;

    float index_medium;
    float transparency;

    // lecture .ppm

    int nH, nW, nTaille, nH2,nW2,nTaille2;


    MaterialType type;
    std::string texture_name;
    std::string normal_map_name;

    OCTET *ArrayImageForTexture;
    OCTET *ArrayImageForNormalMap;

    // Default constructor
    Material() : type(Material_Diffuse_Blinn_Phong), transparency(0.0), index_medium(1.0),
                 ambient_material(Vec3(0., 0., 0.)) {
    }

    Material(std::string texture_name, Vec3 diffuse_material, Vec3 specular_material, Vec3 ambient_material,
             double shininess,std::string normal_map_name) : texture_name(texture_name), diffuse_material(diffuse_material),
                                 specular_material(specular_material), ambient_material(ambient_material),
                                 shininess(shininess), normal_map_name(normal_map_name) {
        this->texture_name = texture_name;
        this->normal_map_name = normal_map_name;
        this->diffuse_material = diffuse_material;
        this->specular_material = specular_material;
        this->ambient_material = ambient_material;
        this->shininess = shininess;
        //conversion string en char*
        char texture_name_char[texture_name.size() + 1];
        strcpy(texture_name_char, texture_name.c_str());

        char normal_map_name_char[normal_map_name.size() + 1];
        strcpy(normal_map_name_char, normal_map_name.c_str());

        //recuperation de hauteur/largeur
        lire_nb_lignes_colonnes_image_ppm(texture_name_char, &nH, &nW);

        lire_nb_lignes_colonnes_image_ppm(normal_map_name_char, &nH2, &nW2);

       //allocation tableau
        nTaille = nH * nW;
        nTaille2 = nH2 * nW2;
        OCTET *ImgIn;
        OCTET *ImgIn2;
        allocation_tableau(ImgIn, OCTET, nTaille*3);
        allocation_tableau(ImgIn2,OCTET,nTaille2*3);
        allocation_tableau(ArrayImageForTexture, OCTET, nTaille*3);
        allocation_tableau(ArrayImageForNormalMap,OCTET,nTaille2*3);
        lire_image_ppm(texture_name_char, ImgIn, nH * nW);
        lire_image_ppm(normal_map_name_char,ImgIn2,nH2*nW2);



        //lecture
        // Print image data into ArrayImageForTexture
        //  i : [u,v]


        // pixel i {r,g,b} -> ArrayImageForTexture[i] = r, ArrayImageForTexture[i+1] = g, ArrayImageForTexture[i+2] = b
        for (int i = 0; i < nTaille * 3; i++) {
            ArrayImageForTexture[i] = ImgIn[i];
        }

        for (int i = 0; i < nTaille2 * 3; i++) {
            ArrayImageForNormalMap[i] = ImgIn2[i];
        }


        free(ImgIn);
        free(ImgIn2);

    }

    Material(std::string texture_name, Vec3 diffuse_material, Vec3 specular_material, Vec3 ambient_material,
             double shininess) : texture_name(texture_name), diffuse_material(diffuse_material),
                                                             specular_material(specular_material), ambient_material(ambient_material),
                                                             shininess(shininess) {
        this->texture_name = texture_name;
        this->diffuse_material = diffuse_material;
        this->specular_material = specular_material;
        this->ambient_material = ambient_material;
        this->shininess = shininess;
        //conversion string en char*
        char texture_name_char[texture_name.size() + 1];
        strcpy(texture_name_char, texture_name.c_str());


        //recuperation de hauteur/largeur
        lire_nb_lignes_colonnes_image_ppm(texture_name_char, &nH, &nW);



        //allocation tableau
        nTaille = nH * nW;

        OCTET *ImgIn;

        allocation_tableau(ImgIn, OCTET, nTaille*3);

        allocation_tableau(ArrayImageForTexture, OCTET, nTaille*3);

        lire_image_ppm(texture_name_char, ImgIn, nH * nW);




        //lecture
        // Print image data into ArrayImageForTexture
        //  i : [u,v]


        // pixel i {r,g,b} -> ArrayImageForTexture[i] = r, ArrayImageForTexture[i+1] = g, ArrayImageForTexture[i+2] = b
        for (int i = 0; i < nTaille * 3; i++) {
            ArrayImageForTexture[i] = ImgIn[i];
        }

        free(ImgIn);

    }
};


#endif // MATERIAL_H
