#ifndef __VIZ__
#define __VIZ__
#include <cstdlib>
#include "opencv2/opencv.hpp"
#include <cstring>
#include <vector>
//#include "highgui.h"
//#include "cv.h"

#ifdef PLY
#include "ply.hpp"
#endif

using namespace cv;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
//#include <GL/freeglut.h>
#endif
#include <stdio.h>
#include <math.h>

using namespace std;
#define forMat(i, j, mat) for (int i = 0; i < mat.rows; i++) for (int j = 0; j < mat.cols; j++)
class Viz {
  private:
    static void initOpenGL();
    static void mouseClicked(int , int , int , int );
    static void mouseMoved(int , int );
    static void reshape(int, int);
    static void keyboardUp(unsigned char key, int x, int y);
    static void keyboard(unsigned char key, int x, int y);
    static float oldAng[2];
    static int mx, my;
    static float oldCameraDistance;
    static int mouseButton;
    static float pointSize;
    static void (*displayFunc)();
    static void (*modelviewFunc)();
    static void (*initFunc)();
    static bool isShift;
    static bool isCtrl;
  public:
    static void findNearFarPlane(vector<GLfloat> &v, GLdouble *mv, float &near, float &far);
    static float ox, oy, oz; // Object central position.
    static float oox, ooy, ooz; // Object central position.
    static void saveImage(std::string fn);
    static void bindTexture(Mat &image, GLuint *text); 
    static void bindTexture(Mat &image); 
    static void setTexture(Mat &img); 
    static void updateLight();
    static Mat getImage(); 
    static void saveMask(const char *fn);
    static void startWindow(int w_, int h_, int x_ = 100, int y_ = 100);
    static void setCenter(float, float, float);
    static void centerMedian();
    static void display();
    static void setDisplayCallback(void (*func)());
    static void setModelviewCallback(void (*func)());
    static void setInitCallback(void (*func)());
    static void setKeyboardCallback(void (*func)(unsigned char key, int x, int y));
    static void (*keyboardFunc)(unsigned char key, int x, int y);
    static void setViewportWidthHeight(int, int);
    static void issueNormal(double mult = 1);
    static void issueNormal(std::vector<GLfloat> &normal, std::vector<GLfloat> &vertices, std::vector<GLuint> &indices, double mult = 1);
    static void updateProjection();
    static void setLightPosition(double a, double b);
    static void toObj(std::string filename, std::vector<GLfloat> &vertices, std::vector<GLuint> &indices, std::vector<GLfloat> &normal);
    static void toObj(std::string filename);
#ifdef PLY
    static void toPly(std::string filename, std::vector<GLfloat> &vertices, std::vector<GLuint> &indices, std::vector<GLfloat> &normal, std::vector<GLubyte> &colors);
    static void toPly(std::string filename);
#endif
    static char keyStates[256];
    static bool drawMesh, useColor, useTexture, shade;
    static int cameraMode;
    static int x, y, w, h;
    static int vw, vh;
    static float cx, cy, cz; // Camera position
    static float angDiv;
    static float cameraDistance;
    static float ang[2];
    static float lang[2];
    static float fov, oldFov;
    static float lx, ly, lz;
    static void loadPly(std::string filename);
    static void loadPly(std::string filename, std::vector<GLfloat> &vertices, std::vector<GLuint> &indices, std::vector<GLubyte> &colors);
    static void loadModelXYZ(Mat &model);
    static void loadModelXYZ(Mat &model, Mat &id);
    static void loadModelXYZTexture(Mat &model, Mat &texture);
    static void loadModelXYZ(Mat &model, Mat _texture, std::vector<GLfloat> &vertices, std::vector<GLuint> &indices, Mat &id);
    static void normalizeSurface(std::vector<GLfloat> &v); 
    static void normalizeSurface(); 
    static void setMaterial(GLenum pname, float r, float g, float b, float a);
    static void setMaterial(GLenum pname, float *c, float a);
    static void setMaterial(GLenum pname, float a); 

    static std::vector<GLfloat> vertices;
    static std::vector<GLuint> indices;
    static std::vector<GLfloat> normals;
    static std::vector<GLubyte> colors;
    static std::vector<GLfloat> textureCoords;

    static float normScale, normTranslate[3];
    static GLuint texture;
    static Mat textureMat;
};


#endif
