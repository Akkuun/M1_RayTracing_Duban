// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/Scene.h"
#include <GL/glut.h>
#include "src/matrixUtilities.h"
#include "src/Kd_tree.h"
#include "src/progressbar.cpp"

using namespace std;





// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 480;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX = 0, lastY = 0, lastZoom = 0;
static unsigned int FPS = 0;
static bool fullScreen = false;

std::vector<Scene> scenes;
unsigned int selected_scene;

std::vector<std::pair<Vec3, Vec3> > rays;

#include <ctime>
#include <iomanip>

MatrixUtilities matrixUtilities;
progressbar progressbar;


void Clear() {
#if defined _WIN32
    system("cls");
    //clrscr(); // including header file : conio.h
#elif defined (__LINUX__) || defined(__gnu_linux__) || defined(__linux__)
    system("clear");
    //std::cout<< u8"\033[2J\033[1;1H"; //Using ANSI Escape Sequences
#elif defined (__APPLE__)
    system("clear");
#endif
}


void printUsage() {
    Clear();
    cerr << endl
         << "gMini: a minimal OpenGL/GLUT application" << endl
         << "for 3D graphics." << endl
         << "Author : Tamy Boubekeur (http://www.labri.fr/~boubek)" << endl << endl
         << "Usage : ./gmini [<file.off>]" << endl
         << "Keyboard commands" << endl
         << "------------------" << endl
         << " ?: Print help" << endl
         << " w: Toggle Wireframe Mode" << endl
         << " g: Toggle Gouraud Shading Mode" << endl
         << " f: Toggle full screen mode" << endl
         << " <drag>+<left button>: rotate model" << endl
         << " <drag>+<right button>: move model" << endl
         << " <drag>+<middle button>: zoom" << endl
         << " q, <esc>: Quit" << endl << endl;

    std::cout << "\033[1;31mRESIZE THE TERMINAL TO NOT HAVE DISPLAY ERRORS\033[0m" << std::endl;

}

void usage() {
    printUsage();
    exit(EXIT_FAILURE);
}


// ------------------------------------
void initLight() {
    GLfloat light_position[4] = {0.0, 1.5, 0.0, 1.0};
    GLfloat color[4] = {1.0, 1.0, 1.0, 1.0};
    GLfloat ambient[4] = {1.0, 1.0, 1.0, 1.0};

    glLightfv(GL_LIGHT1, GL_POSITION, light_position);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT1, GL_SPECULAR, color);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
}

void init() {

    camera.resize(SCREENWIDTH, SCREENHEIGHT);
    initLight();
    //glCullFace (GL_BACK);
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
}


// ------------------------------------
// Replace the code of this
// functions for cleaning memory,
// closing sockets, etc.
// ------------------------------------

void clear() {

}

// ------------------------------------
// Replace the code of this
// functions for alternative rendering.
// ------------------------------------


void draw() {
    glEnable(GL_LIGHTING);
    scenes[selected_scene].draw();

    // draw rays : (for debug)
    //  std::cout << rays.size() << std::endl;
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(6);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (auto &ray: rays) {
        glVertex3f(ray.first[0], ray.first[1], ray.first[2]);
        glVertex3f(ray.second[0], ray.second[1], ray.second[2]);
    }
    glEnd();
}

void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

void idle() {
    static float lastTime = glutGet((GLenum) GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet((GLenum) GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle[64];
        sprintf(winTitle, "Raytracer - FPS: %d", FPS);
        glutSetWindowTitle(winTitle);
        lastTime = currentTime;
    }
    glutPostRedisplay();
}


void ray_trace_from_camera() {
    int w = glutGet(GLUT_WINDOW_WIDTH), h = glutGet(GLUT_WINDOW_HEIGHT);
    std::cout << "Ray tracing a " << w << " x " << h << " image" << std::endl;
    camera.apply();
    clock_t start = clock();
    clock_t lastPrintTime = start;
    matrixUtilities.updateBool();
    matrixUtilities.updateMatrice();
    Vec3 pos, dir;

    //unsigned int nsamples = 125; // UNCOMMENT FOR HD RENDER
    unsigned int nsamples = 10; // SIMPLE RENDER
    std::vector<Vec3> image(w * h, Vec3(0, 0, 0));

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            image[x + y * w] /= nsamples;
            for (unsigned int s = 0; s < nsamples; ++s) {
                float u = ((float) (x) + (float) (rand()) / (float) (RAND_MAX)) / w;
                float v = ((float) (y) + (float) (rand()) / (float) (RAND_MAX)) / h;
                matrixUtilities.screen_space_to_world_space_ray(u, v, pos, dir);
                Vec3 color = scenes[selected_scene].rayTrace(Ray(pos, dir));
                image[x + y * w] += color;
            }
            image[x + y * w] /= nsamples;
        }

        // ESTIMATE REMAINING TIME
        clock_t currentTime = clock();
        float elapsedTime = (float) (currentTime - start) / CLOCKS_PER_SEC;
        float estimatedTotalTime = (elapsedTime / (y + 1)) * h;
        float remainingTime = estimatedTotalTime - elapsedTime;


        // PROGRESS BAR
        float progress = (float) (y + 1) / h * 100;
        progressbar.update(progress);
        progressbar.print(remainingTime); // CHANGE


        // Print remaining time every second
        if ((currentTime - lastPrintTime) >= CLOCKS_PER_SEC) {


        }
    }

    clock_t end = clock();
    std::cout << "\n\tDone" << std::endl;
    std::cout << "\tRendering time: " << (float) (end - start) / CLOCKS_PER_SEC << "s" << std::endl;

    std::string filename = "./rendu.ppm";
    std::ofstream f(filename.c_str(), std::ios::binary);
    if (f.fail()) {
        std::cout << "Could not open file: " << filename << std::endl;
        return;
    }
    f << "P3" << std::endl << w << " " << h << std::endl << 255 << std::endl;
    for (int i = 0; i < w * h; i++)
        f << (int) (255.f * std::min<float>(1.f, image[i][0])) << " "
          << (int) (255.f * std::min<float>(1.f, image[i][1])) << " "
          << (int) (255.f * std::min<float>(1.f, image[i][2])) << " ";
    f << std::endl;
    f.close();
}


void key(unsigned char keyPressed, int x, int y) {
    Vec3 pos, dir;
    switch (keyPressed) {
        case 'f':
            if (fullScreen) {
                glutReshapeWindow(SCREENWIDTH, SCREENHEIGHT);
                fullScreen = false;
            } else {
                glutFullScreen();
                fullScreen = true;
            }
            break;
        case 'q':
        case 27:
            clear();
            exit(0);
            break;
        case 'w':
            GLint polygonMode[2];
            glGetIntegerv(GL_POLYGON_MODE, polygonMode);
            if (polygonMode[0] != GL_FILL)
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            else
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            break;

        case 'r':
            camera.apply();
            rays.clear();
            ray_trace_from_camera();
            break;
        case '+':
            selected_scene++;
            if (selected_scene >= scenes.size()) selected_scene = 0;
            break;
        default:
            printUsage();
            break;
    }
    idle();
}

void mouse(int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate(x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle();
}

void motion(int x, int y) {
    if (mouseRotatePressed) {
        camera.rotate(x, y);
    } else if (mouseMovePressed) {
        camera.move((x - lastX) / static_cast<float>(SCREENWIDTH), (lastY - y) / static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    } else if (mouseZoomPressed) {
        camera.zoom(float(y - lastZoom) / SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize(w, h);
}


int main(int argc, char **argv) {
    if (argc > 2) {
        printUsage();
        exit(EXIT_FAILURE);
    }
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow("gMini");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    key('?', 0, 0);


    camera.move(0., 0., -3.1);
    // MATRIX UTILITIES -> REDUCE USELESS COMPUTING
    matrixUtilities = MatrixUtilities();
    //SCECNES

    selected_scene = 0;
    scenes.resize(4);

    scenes[0].setup_cornel_box2();
    scenes[1].setup_scene_multiple_spheres();
    scenes[2].setup_debug_refraction();
    scenes[3].setup_cornel_box_mesh();

    //MAIN
    glutMainLoop();
    return EXIT_SUCCESS;
}