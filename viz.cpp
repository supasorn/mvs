#include "viz.h"
using namespace std;
vector<GLfloat> Viz::vertices;
vector<GLuint> Viz::indices;
vector<GLfloat> Viz::normals;
vector<GLubyte> Viz::colors;
vector<GLfloat> Viz::textureCoords;
GLuint Viz::texture;

char Viz::keyStates[256];
bool Viz::drawMesh = true;
bool Viz::useColor = false;
bool Viz::useTexture = false;
bool Viz::shade = false;
bool Viz::isShift;
bool Viz::isCtrl;
int Viz::vw = -1, Viz::vh = -1;
int Viz::w, Viz::h, Viz::x, Viz::y;
float Viz::ang[2] = {M_PI/2,0}, Viz::oldAng[2];
float Viz::lang[2] = {0, 2}; // For pose_results2
float Viz::angDiv = 100;
int Viz::mx, Viz::my;
float Viz::cameraDistance = 0.6, Viz::oldCameraDistance;
int Viz::mouseButton;
float Viz::oox = 0, Viz::ooy = 0, Viz::ooz = 0;
float Viz::ox = 0, Viz::oy = 0, Viz::oz = 0;
float Viz::cx = 0, Viz::cy = 0, Viz::cz = 0;
float Viz::lx, Viz::ly, Viz::lz;
void (*Viz::displayFunc)();
void (*Viz::modelviewFunc)();
void (*Viz::initFunc)();
void (*Viz::keyboardFunc)(unsigned char key, int x, int y);
int Viz::cameraMode;
float Viz::fov = 5, Viz::oldFov;
float Viz::normScale;
float Viz::normTranslate[3];
Mat Viz::textureMat;
void Viz::setLightPosition(double a, double b) {
  lang[0] = a;
  lang[1] = b;
  lx = cos(lang[0]) * cos(lang[1]); 
  ly = sin(lang[1]);
  lz = sin(lang[0]) * cos(lang[1]);
}
void Viz::toObj(string filename, vector<GLfloat> &vertices, vector<GLuint> &indices, vector<GLfloat> &normal) {
  FILE *fo = fopen(filename.c_str(), "w");
  for (int i = 0; i < vertices.size(); i+= 3) {
    fprintf(fo, "v %lf %lf %lf\n", vertices[i], vertices[i+1], vertices[i+2]);
  }
  for (int i = 0; i < indices.size(); i+= 3) {
    fprintf(fo, "f %d %d %d\n", indices[i] + 1, indices[i+1] + 1, indices[i+2] + 1);
  }
  fclose(fo);
}
void Viz::toObj(string filename) {
  Viz::toObj(filename, Viz::vertices, Viz::indices, Viz::normals);
}

void Viz::toPly(string filename, vector<GLfloat> &vertices, vector<GLuint> &indices, vector<GLfloat> &normal, vector<GLubyte> &colors) {
  FILE *fo = fopen(filename.c_str(), "w");
  fprintf(fo, "ply\nformat ascii 1.0\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face %d\nproperty list uchar int vertex_indices\nend_header\n", vertices.size() / 3, indices.size() / 3);
  for (int i = 0; i < vertices.size(); i+= 3) {
    fprintf(fo, "%lf %lf %lf %u %u %u\n", vertices[i], vertices[i+1], vertices[i+2], colors[i], colors[i+1], colors[i+2]);
  }
  for (int i = 0; i < indices.size(); i+= 3) {
    fprintf(fo, "3 %d %d %d\n", indices[i] , indices[i+1], indices[i+2]);
  }
  fclose(fo);
}
void Viz::toPly(string filename) {
  Viz::toPly(filename, Viz::vertices, Viz::indices, Viz::normals, Viz::colors);
}

void Viz::issueNormal(vector<GLfloat> &normal, vector<GLfloat> &vertices, vector<GLuint> &indices, double mult) {
  normal.clear();
  int mxId = 0;
  for (int i = 0; i < indices.size(); i++) 
    if (indices[i] > mxId)
      mxId = indices[i];

  vector<Vec3d> norms(mxId + 1);
  vector<int> count(mxId + 1);
  for (int i = 0; i < indices.size(); i+= 3) {
    unsigned int *p = &indices[i];
    Vec3d a(vertices[p[0] * 3], vertices[p[0] * 3 + 1], vertices[p[0] * 3 + 2]);
    Vec3d b(vertices[p[1] * 3], vertices[p[1] * 3 + 1], vertices[p[1] * 3 + 2]);
    Vec3d c(vertices[p[2] * 3], vertices[p[2] * 3 + 1], vertices[p[2] * 3 + 2]);
    Vec3d ab = b - a;
    Vec3d ac = c - a;
    Vec3d cp = ab.cross(ac);
    cp /= norm(cp);

    for (int j = 0; j < 3; j++) {
      //count[p[j]]++;
      //slerp(norms[p[j]], cp, norms[p[j]], 1.0 / count[p[j]]);
      norms[p[j]] += cp;
      count[p[j]] ++;
    }
  }
  for (int i = 0; i < norms.size(); i++) {
    norms[i] /= norm(norms[i]);
    norms[i] *= mult;
    normal.push_back(norms[i][0]);
    normal.push_back(norms[i][1]);
    normal.push_back(norms[i][2]);
  }
}
void Viz::issueNormal(double mult) {
  Viz::issueNormal(Viz::normals, Viz::vertices, Viz::indices, mult);
}
void Viz::setViewportWidthHeight(int w, int h) {
  vw = w;
  vh = h;
}

Mat Viz::getImage() {
  display();
  glReadBuffer(GL_FRONT);
  Mat m(h, w, CV_32FC3);
  glReadPixels(0, 0, w, h, GL_BGR, GL_FLOAT, &m.at<Vec3f>(0, 0));
  flip(m, m, 0);
  m.convertTo(m, CV_64FC3);
  return m;
}
void Viz::saveImage(string fn) {
  display();
  glReadBuffer(GL_FRONT);
  Mat m(h, w, CV_32FC3);
  glReadPixels(0, 0, w, h, GL_BGR, GL_FLOAT, &m.at<Vec3f>(0, 0));
  flip(m, m, 0);
  convertScaleAbs(m, m, 255, 0);
  imwrite(fn, m);
}
void Viz::saveMask(const char *fn) {
  GLfloat global_ambient[] = {1.0,1.0,1.0,1.0};
  
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
  glColor3ub(255,255,255);
  glDisable(GL_COLOR_MATERIAL);
  
  saveImage(fn);
  
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  GLfloat global_ambient2[] = {0.2, 0.2, 0.2, 1.0};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient2);
  glEnable(GL_COLOR_MATERIAL);
  
}
void Viz::setInitCallback(void (*func)()) {
  initFunc = func;
}
void Viz::setKeyboardCallback(void (*func)(unsigned char key, int x, int y)) {
  keyboardFunc = func;
}
void Viz::setModelviewCallback(void (*func)()) {
  modelviewFunc = func;
}
void Viz::setDisplayCallback(void (*func)()) {
  displayFunc = func;
}
void Viz::reshape(int w_, int h_) {
  w = w_;
  h = h_;
  glViewport (0, 0, (GLsizei)vw, (GLsizei)vh);
  updateProjection();
}
void Viz::setCenter(float x, float y, float z) {
  ox = x;
  oy = y;
  oz = z;
}

void Viz::centerMedian() {
  vector<double> xs, ys, zs;
  for (int i = 0; i < Viz::vertices.size(); i+= 3) {
    xs.push_back(Viz::vertices[i]);
    ys.push_back(Viz::vertices[i+1]);
    zs.push_back(Viz::vertices[i+2]);
  }
  sort(xs.begin(), xs.end());
  sort(ys.begin(), ys.end());
  sort(zs.begin(), zs.end());
  setCenter(xs[xs.size() / 2], ys[ys.size() / 2], zs[zs.size() / 2]);
}

void Viz::updateLight() {
  GLfloat light_position[] = { lx, ly, lz, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}
void Viz::setMaterial(GLenum pname, float r, float g, float b, float a) {
  float mat[4] = {r, g, b, a};
  glMaterialfv(GL_FRONT, pname, mat);
}

void Viz::setMaterial(GLenum pname, float *c, float a) {
  float mat[4] = {c[0] * a, c[1] * a, c[2] * a, a};
  glMaterialfv(GL_FRONT, pname, mat);
}
void Viz::setMaterial(GLenum pname, float a) {
  float mat[4] = {a, a, a, a};
  glMaterialfv(GL_FRONT, pname, mat);
}
void Viz::updateProjection() {
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  float mult = 1.0/tan((float)fov/2*M_PI/180) * cameraDistance;
  float minDst = mult - 5 < 0.1 ? 0.1 : mult - 5;
  float maxDst = mult + 5;
  gluPerspective (fov, (GLfloat)vw / (GLfloat)vh, minDst, maxDst );
  //glOrtho(-1, 1, -1, 1, -20, 20);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  cx = ox + mult * cos(ang[0]) * cos(ang[1]); 
  cy = oy + mult * sin(ang[1]);
  cz = oz + mult * sin(ang[0]) * cos(ang[1]);
  gluLookAt(cx, cy, cz, ox, oy, oz, 0, cos(ang[1])>0?1:-1, 0);

  updateLight();
  if(modelviewFunc)
    modelviewFunc();
}


void Viz::bindTexture(Mat &image, GLuint *text) {
  if (image.type() != CV_8UC3) {
    printf("Viz::bindTexture image not CV_8UC3");
    exit(0);
  }
  glGenTextures(1, text);

  glBindTexture( GL_TEXTURE_2D, *text ); //bind the texture to it's array
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
 
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
 
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image.cols, image.rows, 0,  GL_BGR, GL_UNSIGNED_BYTE, &image.at<Vec3b>(0, 0)[0]);
}

void Viz::bindTexture(Mat &image) {
  Viz::bindTexture(image, &Viz::texture);
}

void Viz::setTexture(Mat &img) {
  textureMat = img;
}

void Viz::initOpenGL() {

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();
  setLightPosition(lang[0], lang[1]);

  const float ambient = 0.2;
  GLfloat lmodel_ambient[] = { ambient, ambient, ambient, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  //setMaterial(GL_AMBIENT, 0.2); // default light ambient is 0, so this doesn't do anything unless light is specified
  setMaterial(GL_DIFFUSE, 0.9);
  setMaterial(GL_SPECULAR, 0.0);
  GLfloat mat_shininess[] = { 50.0 };
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  glDisable(GL_COLOR_MATERIAL);


  updateProjection();

}

void Viz::display() {
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  if (drawMesh) {
    if (normals.size()) {
      glEnableClientState(GL_NORMAL_ARRAY);
      glNormalPointer(GL_FLOAT, 0, &Viz::normals[0]);
    }
    if (useColor) {
      glEnableClientState(GL_COLOR_ARRAY);
      glColorPointer(3, GL_UNSIGNED_BYTE, 0, &Viz::colors[0]);
      glEnable(GL_COLOR_MATERIAL);
    } else {
      glDisable(GL_COLOR_MATERIAL);
    }

    if (useTexture) {
      glEnable(GL_TEXTURE_2D);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      glTexCoordPointer(2, GL_FLOAT, 0, &Viz::textureCoords[0]);
    } else {
      glDisable(GL_TEXTURE_2D);
    }

    if (shade) {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
    } else {
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, &Viz::vertices[0]);
    glDrawElements(GL_TRIANGLES, Viz::indices.size(), GL_UNSIGNED_INT, &Viz::indices[0]);
  }

  if(displayFunc) {
    displayFunc();
  }
  
  glutSwapBuffers();
}

void Viz::keyboardUp(unsigned char key, int x, int y) {
  keyStates[key] = 0;
  isCtrl = glutGetModifiers() & GLUT_ACTIVE_CTRL;
  isShift = glutGetModifiers() & GLUT_ACTIVE_SHIFT;
}
void Viz::keyboard(unsigned char key, int x, int y) {
  keyStates[key] = 1;
  isCtrl = glutGetModifiers() & GLUT_ACTIVE_CTRL;
  isShift = glutGetModifiers() & GLUT_ACTIVE_SHIFT;
  if(keyboardFunc) {
    keyboardFunc(key, x, y);
  }
  printf("%d\n", key);
  if (key == 20) { // Ctrl + t
    Viz::useColor ^= 1;
    printf("useColor = %d\n", Viz::useColor);
  } else if (key == 25) { // Ctrl + y
    Viz::shade ^= 1;
    printf("useShade = %d\n", Viz::shade);
  } else if (key == 21) { // Ctlr + u
    Viz::useTexture ^= 1;
    printf("useTexture = %d\n", Viz::useTexture);
  }

}

void Viz::mouseMoved(int x, int y) {
  if (mouseButton == GLUT_LEFT_BUTTON) {
    ang[0] = oldAng[0] + (x-mx) / angDiv;
    ang[1] = oldAng[1] + (y-my) / angDiv;
    
    while(ang[1] > M_PI) ang[1] -= 2*M_PI;
    while(ang[1] < -M_PI) ang[1] += 2*M_PI;
    
  } else if (mouseButton == GLUT_RIGHT_BUTTON) {
    if (isShift) {
      setLightPosition(1.5*M_PI - 2*M_PI * x / w, M_PI - 2*M_PI * y / h);
      updateLight();
    } else if (isCtrl) {
      fov = oldFov + (y - my) / 10.0;
      if (fov < 0) {
        fov = 0.001;
      }
      printf("Fov = %f\n", fov);
    } else {
      cameraDistance = oldCameraDistance + (y-my) / 100.;
    }
  } else if (mouseButton = GLUT_MIDDLE_BUTTON) {
    //ox = oox + (x - mx) / 20.0;
    //oy = ooy + (y - my) / 20;
  }
  updateProjection();
}

void Viz::mouseClicked(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      mx = x;
      my = y;
      oldAng[0] = ang[0];
      oldAng[1] = ang[1];
    } 
    mouseButton = GLUT_LEFT_BUTTON;
  } else if(button == GLUT_RIGHT_BUTTON) {
    isShift = glutGetModifiers() & GLUT_ACTIVE_SHIFT;
    isCtrl = glutGetModifiers() & GLUT_ACTIVE_CTRL;
    if(isShift) {
      oldAng[0] = lang[0];
      oldAng[1] = lang[1];
    } else if (isCtrl) {
      oldFov = fov;
    } else {
      oldCameraDistance = cameraDistance;
    } 
    mx = x;
    my = y;
    mouseButton = GLUT_RIGHT_BUTTON;
  } else if (button == GLUT_MIDDLE_BUTTON) {
    if (state == GLUT_DOWN) {
      mx = x;
      my = y;
      oox = ox;
      ooy = oy;
      ooz = oz;
    } 
    mouseButton = GLUT_MIDDLE_BUTTON;
  }
}

void Viz::startWindow(int w_, int h_, int x_, int y_) {
  x = x_;
  y = y_;
  w = w_;
  h = h_;
  if(vw == -1) { 
    vw = w;
    vh = h;
  }
  int tmp = 0;
  char **tmpv;
  glutInit (&tmp, tmpv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
  glutInitWindowSize (w, h);
  glutInitWindowPosition (x, y);
  glutCreateWindow ("Viz :)");
  
  glutDisplayFunc (display);
  glutIdleFunc (display);
  glutMouseFunc(mouseClicked);
  glutMotionFunc(mouseMoved);
	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboardUp);
  glutReshapeFunc (reshape);
  initOpenGL();
  if (initFunc)
    initFunc();
  if (textureMat.rows) 
    Viz::bindTexture(textureMat);
  glutMainLoop ();
  
}

void Viz::findNearFarPlane(vector<GLfloat> &v, GLdouble *mv, float &near, float &far) {
  float mn = 1e10, mx = -1e10;
  for(int i=0; i<v.size(); i+=3) {
    float z = mv[2] * v[i] + mv[6] * v[i+1] + mv[10] * v[i+2] + mv[14];
    if(z < mn) mn = z;
    if(z > mx) mx = z;
  }
  near = -mx;
  far = -mn;
  float range = far - near;
  if (range < 1) range = 1;
  far += range;
  near -= range;
}

vector<GLfloat> *_vertices;
vector<GLubyte> *_colors;
vector<GLuint> *_indices;

void vertex_callback(ply::float32 x) {
  _vertices->push_back(x);
}
void vertex_color_callback(ply::uint8 x) {
  _colors->push_back(x);
}

template  <typename ScalarType>
std::tr1::function <void (ScalarType)> scalar_property_definition_callback(const std::string& element_name, const std::string& property_name);

template  <>
std::tr1::function <void (ply::float32)> scalar_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
  
  if (element_name == "vertex") {
    if (property_name == "x" || property_name == "y" || property_name == "z") {
      return vertex_callback;
    } 
  }
  return 0;
}

template  <>
std::tr1::function <void (ply::uint8)> scalar_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
  
  if (element_name == "vertex") {
    if (property_name == "red" || property_name == "green" || property_name == "blue") {
      return vertex_color_callback;
    }
  }
  return 0;
}

void face_vertex_indices_begin(ply::uint8 size) {
  if (size != 3) {
    printf("face doesn't have 3 vertices\n");
    exit(0);
  }
}

void face_vertex_indices_element(ply::int32 a) {
  _indices->push_back(a);
}

void face_vertex_indices_end() {
}

template <typename SizeType, typename ScalarType> std::tr1::tuple<std::tr1::function<void (SizeType)>, std::tr1::function<void (ScalarType)>, std::tr1::function<void ()> > list_property_definition_callback(const std::string& element_name, const std::string& property_name);

template <>
std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> > list_property_definition_callback(const std::string& element_name, const std::string& property_name)
{
  if ((element_name == "face") && (property_name == "vertex_indices")) {
    return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> >(
      face_vertex_indices_begin,
      face_vertex_indices_element,
      face_vertex_indices_end
    );
  }
  else {
    return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> >(0, 0, 0);
  }
}

void Viz::loadPly(string filename, vector<GLfloat> &vertices, vector<GLuint> &indices, vector<GLubyte> &colors) {
  vertices.clear();
  colors.clear();
  indices.clear();
  _vertices = &vertices;
  _colors = &colors;
  _indices = &indices;
  ply::ply_parser ply_parser;

  ply::ply_parser::scalar_property_definition_callbacks_type scalar_property_definition_callbacks;
  ply::at <ply::float32>(scalar_property_definition_callbacks) = scalar_property_definition_callback <ply::float32>;
  ply::at<ply::uint8>(scalar_property_definition_callbacks) = scalar_property_definition_callback<ply::uint8>;
  ply_parser.scalar_property_definition_callbacks(scalar_property_definition_callbacks);

  ply::ply_parser::list_property_definition_callbacks_type list_property_definition_callbacks;
  ply::at<ply::uint8, ply::int32>(list_property_definition_callbacks) = list_property_definition_callback<ply::uint8, ply::int32>;
  ply_parser.list_property_definition_callbacks(list_property_definition_callbacks);
  ply_parser.parse(filename);
}

void Viz::loadPly(string filename) {
  Viz::loadPly(filename, Viz::vertices, Viz::indices, Viz::colors);
}

void Viz::loadModelXYZ(Mat &model) {
  Mat id;
  loadModelXYZ(model, Mat(), Viz::vertices, Viz::indices, id);
}

void Viz::loadModelXYZ(Mat &model, Mat &id) {
  loadModelXYZ(model, Mat(), Viz::vertices, Viz::indices, id);
}

void Viz::loadModelXYZTexture(Mat &model, Mat &texture) {
  Mat id;
  loadModelXYZ(model, texture, Viz::vertices, Viz::indices, id);
}

void Viz::loadModelXYZ(Mat &model, Mat _texture, vector<GLfloat> &vertices, vector<GLuint> &indices, Mat &id) {
  
  Mat texture;
  if (_texture.rows) {
    if (_texture.type() == CV_64FC3) {
      _texture.convertTo(texture, CV_8UC3, 255);
    } else if (_texture.type() == CV_8UC3) {
      texture = _texture;
    }
  }

  colors.clear();
  vertices.clear();
  indices.clear();
  Mat mask = Mat::zeros(model.size(), CV_8UC1);
  forMat (i, j, mask) {
    mask.at<uchar>(i, j) = model.at<Vec3d>(i, j) != Vec3d(0, 0, 0); 
  }
  id = Mat::zeros(mask.size(), CV_32S);
  int idCount = 0;
  forMat (i, j, mask) {
    if (mask.at<uchar>(i, j)) {
      id.at<int>(i, j) = idCount++;
      vertices.push_back(model.at<Vec3d>(i, j)[0]);
      vertices.push_back(model.at<Vec3d>(i, j)[1]);
      vertices.push_back(model.at<Vec3d>(i, j)[2]);
    } 
  }

  if (_texture.rows) {
    colors.resize(vertices.size());
    forMat (i, j, texture) {
      if (mask.at<uchar>(i, j)) {
        colors[id.at<int>(i, j) * 3 + 2] = texture.at<Vec3b>(i, j)[0];
        colors[id.at<int>(i, j) * 3 + 1] = texture.at<Vec3b>(i, j)[1];
        colors[id.at<int>(i, j) * 3 + 0] = texture.at<Vec3b>(i, j)[2];
      }
    }
  }

  forMat (i, j, mask) {
    if(mask.at<uchar>(i, j) && i+1 < mask.rows && j+1 < mask.cols && 
        mask.at<uchar>(i + 1, j + 1) &&
        mask.at<uchar>(i, j + 1) &&
        mask.at<uchar>(i + 1, j)) {
      indices.push_back(id.at<int>(i, j));
      indices.push_back(id.at<int>(i+1, j+1));
      indices.push_back(id.at<int>(i, j+1));

      indices.push_back(id.at<int>(i, j));
      indices.push_back(id.at<int>(i+1, j));
      indices.push_back(id.at<int>(i+1, j+1));
    }
  }
}

void Viz::normalizeSurface(vector<GLfloat> &v) {
  float c[3][2];
  float r[3];
  float mr;
  c[0][0] = c[1][0] = c[2][0] = 1.0e10; 
  c[0][1] = c[1][1] = c[2][1] = -1.0e10; 
  for(int i = v.size()-1; i >= 0; i--) {
    if(v[i] < c[i%3][0]) c[i%3][0] = v[i];
    if(v[i] > c[i%3][1]) c[i%3][1] = v[i];
  }
  r[0] = c[0][1] - c[0][0];
  r[1] = c[1][1] - c[1][0];
  r[2] = c[2][1] - c[2][0];

  mr = r[0];
  if(r[1] > mr) mr = r[1];
  if(r[2] > mr) mr = r[2];

  c[0][0] = c[1][0] = c[2][0] = 1.0e10; 
  c[0][1] = c[1][1] = c[2][1] = -1.0e10; 
  Viz::normScale = 1 / mr;
  for(int i = v.size()-1; i >= 0 ;i--) {
    v[i] *= Viz::normScale;
    if(v[i] < c[i%3][0]) c[i%3][0] = v[i];
    if(v[i] > c[i%3][1]) c[i%3][1] = v[i];
  }

  for (int i = 0; i < 3; i++) {
    Viz::normTranslate[i] = -0.5 * (c[i][0] + c[i][1]);
  }
  for(int i = v.size()-1; i >= 0; i--) {
    v[i] += Viz::normTranslate[i%3];
  }
}

void Viz::normalizeSurface() {
  Viz::normalizeSurface(Viz::vertices);
}

#ifdef MAIN
int main(int argc, char **argv) {
  Viz::loadPly("/home/supasorn/Downloads/PublicMM1/03_scans_ply/00001_20061015_00418_neutral_face05.ply");
  Viz::issueNormal();
  Viz::normalizeSurface();
  Viz::startWindow(500, 500);
}
#endif

