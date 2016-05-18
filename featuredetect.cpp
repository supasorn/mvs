#define UTIL_CERES
#include "utilities.h"
#include "viz.h"
using namespace cv;
using namespace std;
DEFINE_string(dataset, "/home/supasorn/mvs/datasets/templeRing", "a");
int n;

struct Frame {
  string name;
  // Projection Matrix = k * [r : t] = [kr : kt]
  Mat k, r, t, rt;  
  Mat p;

  Mat img;
  Mat camera; // camera center
  Mat u, v; // right and down vector

  Mat fourCorners;
};
vector<Frame> frames;


vector<Point3f> testPoints;
vector<Point3f> testColors;

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

void drawEpiline(Point2f p, int i0, int i1) {
  vector<Point2f> ps;
  ps.push_back(p);

  vector<Vec3f> lines;
  getEpiLine(lines, frames[i0], frames[i1], ps);

  Scalar color(255, 0, 0);

  Mat a = frames[i0].img.clone();
  Mat b = frames[i1].img.clone();
  Point2f p0, p1;
  if (getLineSegment(p0, p1, lines[0], frames[i0].img.size())) {
    line(b, p0, p1, color);
  }

  circle(a, ps[0], 1, color, -1, CV_AA);
  imshow("img0", a);
  imshow("img1", b);
  //waitKey();
}

void detectHarris(Mat &im, vector<KeyPoint> &kp) {
  //Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(5000, 0.0001, 1, 3, true, 0.0001);
  Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(2000, 0.0001, 1, 3, true, 0.0001);
  //Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(2000, 0.0001, 1, 3, true, 0.0001);
  Mat gray;
  cvtColor(im, gray, CV_BGR2GRAY);
  detector->detect(gray, kp);
}

void displayKeyPoint(string name, Mat &im, vector<KeyPoint> &kp) {
  Mat tmp = im.clone();
  for (int i = 0; i < kp.size(); i++) {
    circle(tmp, kp[i].pt, 1, Scalar(0, 0, 255), -1);
  }
  imshow(name, tmp);
}

void matching(Frame &f0, vector<KeyPoint> &kp0, Frame &f1, vector<KeyPoint> &kp1) {
   
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
  //displayKeyPoint("f0", frames[0].img, kp[0]);
  detectHarris(frames[i1].img, kp[1]);
  //displayKeyPoint("f1", frames[1].img, kp[1]);


  vector<Point2f> ps;
  for (int i = 0; i < kp[0].size(); i++) {
    ps.push_back(kp[0][i].pt);
  }
  vector<Vec3f> lines;
  getEpiLine(lines, frames[i0], frames[i1], ps);
  printf("%d\n", lines.size());
  Point2f p0, p1;
  Scalar color(255, 0, 0);

  int n = 500;
  int patchSize = 39;

  vector<Point2f> c0, c1;

  for (int i = 0; i < lines.size(); i++) {
    printf("%d\n", i);
    Mat patch;
    if (!getPatch(patch, frames[i0].img, kp[0][i].pt, patchSize)) continue;

    if (getLineSegment(p0, p1, lines[i], frames[i0].img.size())) {
      //line(frames[1].img, p0, p1, color);

      double best = 1e10;
      Mat bestPatch;
      Point2f bestPoint;
      for (int j = 0; j < n; j++) {
        float t = 1.0 * j / (n - 1); 
        Point2f in = p1 * t + p0 * (1 - t);

        Mat patch2;
        if (!getPatch(patch2, frames[i1].img, in, patchSize)) continue;

        Scalar diff = sum(abs(patch - patch2));
        if (diff[0] + diff[1] + diff[2] < best) {
          best = diff[0] + diff[1] + diff[2];
          //bestPatch = patch2.clone();
          bestPoint = in;
        }

        //circle(frames[1].img, in, 1, color, -1, CV_AA);
        //imshow("f1", frames[1].img);
        //waitKey();
         
      }
      if (best != 1e10) {
        c0.push_back(kp[0][i].pt);
        c1.push_back(bestPoint);
        //Vec3d tmp = frames[0].img.at<Vec3d>(kp[0][i].pt.y, kp[0][i].pt.x);
        //testColors.push_back(Point3f(tmp[0], tmp[1], tmp[2]));
        testColors.push_back(Point3f(frames[i0].img.at<Vec3b>(kp[0][i].pt.y, kp[0][i].pt.x)));
      }

      //imshow("patch", joinH(patch, bestPatch));
      //waitKey();
    }
  }
  //vector<Point3f> tp;
  Mat tp;
  triangulatePoints(frames[i0].p, frames[i1].p, c0, c1, tp);
  //printf("%d %d\n", tp.rows, tp.cols);
  //printType(tp);
  //printf("%f %f %f %f\n", tp.at<float>(0, 0),tp.at<float>(1, 0),tp.at<float>(2, 0),tp.at<float>(3, 0));
  
  for (int i = 0; i < tp.cols; i++) {
    testPoints.push_back(Point3f(
          tp.at<float>(0, i) / tp.at<float>(3, i),
          tp.at<float>(1, i) / tp.at<float>(3, i),
          tp.at<float>(2, i) / tp.at<float>(3, i)
          ));
  }

  //imshow("f1", frames[1].img);
  //waitKey();
}
void mouse(int event, int x, int y, int flags, void* param) {
  drawEpiline(Point2f(x, y), 0, 1);
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

  for (int i = 0; i < testPoints.size(); i++) {
    auto &tp = testPoints[i];
    auto &tc = testColors[i];
    glBegin(GL_POINTS);
    glColor3ub(tc.z, tc.y, tc.x);
    glVertex3f(tp.x, tp.y, tp.z);
    glEnd();
  }
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
    //cout << frames[i].r << endl;
    //cout << frames[i].t << endl;
    //cout << frames[i].rt << endl;
    //cout << frames[i].p << endl;
    //exit(0);


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

    //cout << frames[i].k << endl;
    //imshow("fram", frames[i].img);
    //waitKey();
  }
  fclose(fi);
  imshow("img0", frames[0].img);
  for (int i = 0; i < frames.size(); i++) {
    if (i + 2 < frames.size())
      testGeneratePoints(i, i+2);
  }
  //testGeneratePoints(0, 1);
  //testGeneratePoints(10, 11);
  


  //detectHaris(frames[0].img);


  setMouseCallback("img0", mouse);
  waitKey();
  Viz::setDisplayCallback(display);
  Viz::startWindow(800, 800);

  waitKey();


}
int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);

  loadDataset();
}
