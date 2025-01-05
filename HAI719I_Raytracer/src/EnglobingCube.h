#ifndef ENGLOBINGCUBE_H
#define ENGLOBINGCUBE_H

#include <vector>
#include <string>
#include "Vec3.h"
#include "Triangle.h"
#include <algorithm>
#include <iostream>

#include <GL/glut.h>



class EnglobingCube{

public :
    Vec3 min;
    Vec3 max;
    EnglobingCube(){}

    EnglobingCube(Vec3 minVec, Vec3 maxVec) : min(minVec), max(maxVec){}

    EnglobingCube(std::vector<Triangle> triangles){

        Vec3 minB = triangles[0].getFirstPoint();
        Vec3 maxB = triangles[0].getFirstPoint();


        for(Triangle &t : triangles){


            Vec3 v0 = t.getFirstPoint();
            Vec3 v1 = t.getFirstPoint();
            Vec3 v2 = t.getFirstPoint();

            for (int i = 0; i < 3; i++){
                minB[i] = std::min(minB[i], v0[i]);
                maxB[i] = std::max(maxB[i], v0[i]);
                minB[i] = std::min(minB[i], v1[i]);
                maxB[i] = std::max(maxB[i], v1[i]);
                minB[i] = std::min(minB[i], v2[i]);
                maxB[i] = std::max(maxB[i], v2[i]);

            }


        }

        min = minB;
        max = maxB;

    }


    void draw(Vec3 color) const {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_LINES);

        // Bottom face
        glVertex3d(min[0], min[1], min[2]);
        glVertex3d(max[0], min[1], min[2]);

        glVertex3d(min[0], min[1], min[2]);
        glVertex3d(min[0], max[1], min[2]);

        glVertex3d(max[0], min[1], min[2]);
        glVertex3d(max[0], max[1], min[2]);

        glVertex3d(min[0], max[1], min[2]);
        glVertex3d(max[0], max[1], min[2]);

        // Top face
        glVertex3d(min[0], min[1], max[2]);
        glVertex3d(max[0], min[1], max[2]);

        glVertex3d(min[0], min[1], max[2]);
        glVertex3d(min[0], max[1], max[2]);

        glVertex3d(max[0], min[1], max[2]);
        glVertex3d(max[0], max[1], max[2]);

        glVertex3d(min[0], max[1], max[2]);
        glVertex3d(max[0], max[1], max[2]);

        // Connect the corners
        glVertex3d(min[0], min[1], min[2]);
        glVertex3d(min[0], min[1], max[2]);

        glVertex3d(max[0], min[1], min[2]);
        glVertex3d(max[0], min[1], max[2]);

        glVertex3d(min[0], max[1], min[2]);
        glVertex3d(min[0], max[1], max[2]);

        glVertex3d(max[0], max[1], min[2]);
        glVertex3d(max[0], max[1], max[2]);

        glEnd();
    }


};
#endif