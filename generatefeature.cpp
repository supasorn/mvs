#define UTIL_CERES
#include "utilities.h"
#include "viz.h"
using namespace cv;
using namespace std;
DEFINE_string(dataset, "/home/supasorn/mvs/datasets/templeRing", "dataset");
DEFINE_bool(load, false, "load points");

int n;

struct Frame {
  string name;
  Mat k, r, t, rt;  
  Mat p; // Projection matrix

  Mat img;
  Mat camera; // camera center

  Mat fourCorners;
};
vector<Frame> frames;


vector<Point3f> testPoints;
vector<Point3f> testColors;
vector<pair<int, int> > viscam;

void getEpiLine(vector<Vec3f> &lines, Frame &f0, Frame &f1, vector<Point2f> &ps) {
  Mat r0_ = f0.r.t();

  Mat r = f1.r * r0_;
  Mat t = f1.t - r * f0.t;

  Mat E = (Mat_<double>(3, 3) << 
      0, -t.at<double>(2, 0), t.at<double>(1, 0),
      t.at<double>(2, 0), 0, -t.at<double>(0, 0),
      -t.at<double>(1, 0), t.at<double>(0, 0), 0
      ) * r;

  Mat F = f1.k.inv().t() * E * f0.k.inv();
  computeCorrespondEpilines(ps, 1, F, lines);
}

int getLineSegment(Point2f &p0, Point2f &p1, Vec3f &line, Size sz) {
  int top = 0;
  Point2f vtop;
  int bottom = 0;
  Point2f vbottom;

  if (line[0] != 0) {
    vtop.x = -line[2] / line[0];
    vtop.y = 0;
    top = vtop.x >= 0 && vtop.x <= sz.width;
    vbottom.x = -(line[2] + line[1] * sz.height) / line[0];
    vbottom.y = sz.height;
    bottom = vbottom.x >= 0 && vbottom.x <= sz.width;
  }

  int left = 0;
  Point2f vleft;
  int right = 0;
  Point2f vright;

  if (line[1] != 0) {
    vleft.x = 0;
    vleft.y = -line[2] / line[1];
    left = vleft.y >= 0 && vleft.y <= sz.height;

    vright.x = sz.width;
    vright.y = -(line[2] + line[0] * sz.width) / line[1];
    right = vright.y >= 0 && vright.y <= sz.height;
  }

  int flat[4] = {left, right, top, bottom};
  Point2f *vflat[4] = {&vleft, &vright, &vtop, &vbottom};

  int first = -1;
  float mx = 0;
  int second = -1;
  for (int i = 0; i < 4; i++) {
    float tmp;
    if (first == -1) {
      if (flat[i]) 
        first = i;
    } else if (flat[i] && (tmp = norm(*vflat[first] - *vflat[i])) > mx) {
      mx = tmp;
      second = i;
    }
  }
  if (second != -1) {
    p0 = *vflat[first];
    p1 = *vflat[second];
    return 1;
  }
  return 0;
}


void detectHarris(Mat &im, vector<KeyPoint> &kp) {
  //Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(2000, 0.0001, 1, 3, true, 0.0001);
  Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(5000, 0.0001, 1, 3, true, 0.0001);
  Mat gray;
  cvtColor(im, gray, CV_BGR2GRAY);
  detector->detect(gray, kp);
}

int getPatch(Mat &patch, Mat &im, Point2f &center, int patchSize) {
  if (patchSize % 2)
    patchSize--;
  if (center.x - patchSize / 2 >= 0 && center.y - patchSize / 2 >= 0 &&
      center.x + patchSize / 2 < im.cols && center.y + patchSize / 2 < im.rows) {
    patch = im(Rect(center.x - patchSize / 2, center.y - patchSize / 2, patchSize, patchSize));
    return 1;
  }
  return 0;
}

void testGeneratePoints(int i0, int i1) {
  vector<KeyPoint> kp[2];
  detectHarris(frames[i0].img, kp[0]);
  detectHarris(frames[i1].img, kp[1]);


  vector<Point2f> ps;
  for (int i = 0; i < kp[0].size(); i++) {
    ps.push_back(kp[0][i].pt);
  }

  vector<Vec3f> lines;
  getEpiLine(lines, frames[i0], frames[i1], ps);
  Point2f p0, p1;

  int n = 500;
  int patchSize = 39;

  vector<Point2f> c0, c1;

  for (int i = 0; i < lines.size(); i++) {
    printf("%d\n", i);
    Mat patch;
    if (!getPatch(patch, frames[i0].img, kp[0][i].pt, patchSize)) continue;

    if (getLineSegment(p0, p1, lines[i], frames[i0].img.size())) {

      double best = 1e10;
      Point2f bestPoint;
      // Stupid plane sweep algo
      for (int j = 0; j < n; j++) {
        float t = 1.0 * j / (n - 1); 
        Point2f in = p1 * t + p0 * (1 - t);

        Mat patch2;
        if (!getPatch(patch2, frames[i1].img, in, patchSize)) continue;

        Scalar diff = sum(abs(patch - patch2));
        if (diff[0] + diff[1] + diff[2] < best) {
          best = diff[0] + diff[1] + diff[2];
          bestPoint = in;
        }
      }
      if (best != 1e10) {
        c0.push_back(kp[0][i].pt);
        c1.push_back(bestPoint);
        testColors.push_back(Point3f(frames[i0].img.at<Vec3b>(kp[0][i].pt.y, kp[0][i].pt.x)));
      }
    }
  }
  Mat tp;
  triangulatePoints(frames[i0].p, frames[i1].p, c0, c1, tp);
  
  for (int i = 0; i < tp.cols; i++) {
    testPoints.push_back(Point3f(
          tp.at<float>(0, i) / tp.at<float>(3, i),
          tp.at<float>(1, i) / tp.at<float>(3, i),
          tp.at<float>(2, i) / tp.at<float>(3, i)
          ));
    viscam.push_back(make_pair(i0, i1));
  }
}


void display() {
#define GLM(a) glVertex3f(frames[i].fourCorners.at<double>(0, a),  frames[i].fourCorners.at<double>(1, a),  frames[i].fourCorners.at<double>(2, a))
#define GLC() glVertex3f(frames[i].camera.at<double>(0, 0), frames[i].camera.at<double>(1, 0), frames[i].camera.at<double>(2, 0))

  glBegin(GL_POINTS);
  glColor3ub(0, 255, 0);
  glVertex3f(0, 0, 0);
  glEnd();
  for (int i = 0; i < frames.size(); i++) {
    glBegin(GL_POINTS);
    glColor3ub(255, 0, 0);
    GLC();
    glEnd();

    // Draw camera frustum
    glBegin(GL_LINE_STRIP);
    GLM(2);
    GLM(3);
    GLC();
    GLM(2);
    GLM(0);
    GLC();
    GLM(1);
    GLM(3);
    GLC();
    GLM(0);
    GLM(1);
    glEnd();
  }
#undef GLM
#undef GLC

  for (int i = 0; i < testPoints.size(); i++) {
    auto &tp = testPoints[i];
    auto &tc = testColors[i];
    glBegin(GL_POINTS);
    glColor3ub(tc.z, tc.y, tc.x);
    glVertex3f(tp.x, tp.y, tp.z);
    glEnd();
  }
}

void savePoints() {
  FILE *fo = fopen("point.bin", "wb");
  int n = frames.size(); 
  // Write the number of cameras, and cameras' centers.
  fwrite(&n, sizeof(int), 1, fo);
  for (int i = 0; i < n; i++) {
    fwrite(&frames[i].camera.at<double>(0, 0), sizeof(double), 3, fo);
  }

  n = testPoints.size();
  fwrite(&n, sizeof(int), 1, fo);
  for (int i = 0; i < n; i++) {
    fwrite(&testPoints[i].x, sizeof(float), 1, fo); 
    fwrite(&testPoints[i].y, sizeof(float), 1, fo); 
    fwrite(&testPoints[i].z, sizeof(float), 1, fo); 
    fwrite(&testColors[i].x, sizeof(float), 1, fo); 
    fwrite(&testColors[i].y, sizeof(float), 1, fo); 
    fwrite(&testColors[i].z, sizeof(float), 1, fo); 
    fwrite(&viscam[i].first, sizeof(int), 1, fo);
    fwrite(&viscam[i].second, sizeof(int), 1, fo);
  }
  fclose(fo);
  printf("Done writing\n");
}

void loadPoints() {
  FILE *fi = fopen("point.bin", "rb");
  int n; 
  // Write the number of cameras, and cameras' centers.
  fread(&n, sizeof(int), 1, fi);
  printf("#cameras = %d\n", n);
  for (int i = 0; i < n; i++) {
    fread(&frames[i].camera.at<double>(0, 0), sizeof(double), 3, fi);
  }

  fread(&n, sizeof(int), 1, fi);
  printf("#points = %d\n", n);
  testPoints.resize(n);
  testColors.resize(n);
  viscam.resize(n);
  for (int i = 0; i < n; i++) {
    fread(&testPoints[i].x, sizeof(float), 1, fi); 
    fread(&testPoints[i].y, sizeof(float), 1, fi); 
    fread(&testPoints[i].z, sizeof(float), 1, fi); 
    fread(&testColors[i].x, sizeof(float), 1, fi); 
    fread(&testColors[i].y, sizeof(float), 1, fi); 
    fread(&testColors[i].z, sizeof(float), 1, fi); 
    //printf("%f %f %f\n", testColors[i].x, testColors[i].y, testColors[i].z);
    fread(&viscam[i].first, sizeof(int), 1, fi);
    fread(&viscam[i].second, sizeof(int), 1, fi);
  }
  fclose(fi);
  printf("Done reading\n");
}


void loadDataset() {
  FILE *fi = fopen((FLAGS_dataset + "/templeR_par.txt").c_str(), "r");
  fscanf(fi, "%d", &n);
  frames.resize(n);
  for (int i = 0; i < n; i++) {
    char tmp[200];
    fscanf(fi, " %s ", tmp);
    frames[i].name = tmp;
    frames[i].img = imread(FLAGS_dataset + "/" + tmp);

    frames[i].k = Mat(3, 3, CV_64F);
    frames[i].r = Mat(3, 3, CV_64F);
    frames[i].t = Mat(3, 1, CV_64F);
    frames[i].rt = Mat(3, 4, CV_64F);
    for (int k = 0; k < 9; k++) 
      fscanf(fi, " %lf", &frames[i].k.at<double>(0, 0) + k);
    for (int k = 0; k < 9; k++) 
      fscanf(fi, " %lf", &frames[i].r.at<double>(0, 0) + k);
    for (int k = 0; k < 3; k++) 
      fscanf(fi, " %lf", &frames[i].t.at<double>(0, 0) + k);

    frames[i].camera = - frames[i].r.t() * frames[i].t;
    frames[i].r.copyTo(frames[i].rt(Rect(0, 0, 3, 3)));
    frames[i].t.copyTo(frames[i].rt(Rect(3, 0, 1, 3)));
    frames[i].p = frames[i].k * frames[i].rt;

    // Let projection matrix = P = [F : B] where F is 3 x 3, B 3 x 1;
    Mat F_ = frames[i].r.t() * frames[i].k.inv(); // F-1 = (kr)-1 = r^T k-1
    Mat B = frames[i].k * frames[i].t;
    // Solving Fx + B = [u; v; w] yielding x = F-1 (uvw - B).
    
    float x[4] = {0, 0, frames[i].img.cols - 1, frames[i].img.cols - 1};
    float y[4] = {0, frames[i].img.rows - 1, 0, frames[i].img.rows - 1};
    frames[i].fourCorners = Mat(3, 4, CV_64F);
    for (int j = 0; j < 4; j++) {
      float w = 0.1;
      Mat uvw = (Mat_<double>(3, 1) << x[j] * w, y[j] * w, w);
      Mat xyz = F_ * (uvw - B);
      xyz.copyTo(frames[i].fourCorners(Rect(j, 0, 1, 3)));
    }
  }
  fclose(fi);
  imshow("img0", frames[0].img);
 
  if (FLAGS_load) {
    loadPoints();
  } else {
    for (int i = 0; i < frames.size(); i++) {
      if (i + 1 < frames.size())
        testGeneratePoints(i, i + 1);
    }
    savePoints();
  }

  Viz::setDisplayCallback(display);
  Viz::startWindow(800, 800);

  waitKey();


}
int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);

  loadDataset();
}
