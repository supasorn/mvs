#define UTIL_CERES
#include "utilities.h"
#include "viz.h"
using namespace cv;
using namespace std;
DEFINE_string(dataset, "/home/supasorn/mvs/datasets/templeRing", "a");
int n;

struct Frame {
  string name;
  Mat k, r, t;
  Mat img;
  Mat camera; // camera center
  Mat u, v; // right and down vector
};
vector<Frame> frames;


void drawEpiline(Point2f p, int i0, int i1) {
  Mat r0_ = frames[i0].r.t();

  Mat r = frames[i1].r * r0_;
  Mat t = frames[i1].t - r * frames[i0].t;

  Mat E = (Mat_<double>(3, 3) << 
      0, -t.at<double>(2, 0), t.at<double>(1, 0),
      t.at<double>(2, 0), 0, -t.at<double>(0, 0),
      -t.at<double>(1, 0), t.at<double>(0, 0), 0
      ) * r;

  Mat F = frames[i1].k.inv().t() * E * frames[i0].k.inv();

  vector<Point2f> ps;
  ps.push_back(p);
  std::vector<Vec3f> lines;
  computeCorrespondEpilines(ps, 1, F, lines);

  int top = 0;
  Point2f vtop;
  int bottom = 0;
  Point2f vbottom;
  Size sz = frames[i0].img.size();

  if (lines[0][0] != 0) {
    vtop.x = -lines[0][2] / lines[0][0];
    vtop.y = 0;
    top = vtop.x >= 0 && vtop.x <= sz.width;
    vbottom.x = -(lines[0][2] + lines[0][1] * sz.height) / lines[0][0];
    vbottom.y = sz.height;
    bottom = vbottom.x >= 0 && vbottom.x <= sz.width;
  }

  int left = 0;
  Point2f vleft;
  int right = 0;
  Point2f vright;

  if (lines[0][1] != 0) {
    vleft.x = 0;
    vleft.y = -lines[0][2] / lines[0][1];
    left = vleft.y >= 0 && vleft.y <= sz.height;

    vright.x = sz.width;
    vright.y = -(lines[0][2] + lines[0][0] * sz.width) / lines[0][1];
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
  Scalar color(255, 0, 0);

  Mat a = frames[i0].img.clone();
  Mat b = frames[i1].img.clone();
  if (second != -1) {
    line(b, *vflat[first], *vflat[second], color);
  }

  printf("%d %d %d %d\n", top, bottom, left, right);
  circle(a, ps[0], 1, color, -1, CV_AA);
  imshow("img0", a);
  imshow("img1", b);
  //waitKey();
}

void detectHaris(Mat &im) {
  int blockSize = 2;
  int apertureSize = 3;
  double k = 0.001;
  Mat gray;

  cvtColor(im, gray, CV_BGR2GRAY);
  Mat dst = Mat::zeros(gray.size(), CV_32FC1);
  Mat dst_norm, dst_norm_scaled;
  /// Detecting corners
  cornerHarris( gray, dst, blockSize, apertureSize, k, BORDER_DEFAULT );

  Mat g1, g2, result;
  GaussianBlur(gray, g1, Size(1,1), 0);
  GaussianBlur(gray, g2, Size(3,3), 0);
  result = g1 - g2;
  imshow("dog", abs(result) * 10);

  // NMS
  int nmswindow = 1;
  forMat (i, j, dst) {
    for (int k = -nmswindow; k <= nmswindow; k++) {
      
    }
  }
  //dilate(dst, dst, getElement(2));

  //normalize( dst, dst_norm, 0, 1, NORM_MINMAX, CV_32FC1, Mat() );
  printType(dst);
  dst_norm = dst;
  int count = 0;
  forMat (i, j, dst_norm) {
    if (dst_norm.at<float>(i, j) > 0.00001) {
      dst_norm.at<float>(i, j) = 1;
      count++;
    }
  }
  printf("%d\n", count);
  //convertScaleAbs( dst_norm, dst_norm_scaled );

  //imshow("gray", gray);
  imshow("dst", dst_norm);
  //imshow("s", dst_norm_scaled);
  waitKey();
}

void detectHarris(Mat &im, vector<KeyPoint> &kp) {
  Ptr<FeatureDetector> detector = new GoodFeaturesToTrackDetector(5000, 0.0001, 1, 3, true, 0.0001);
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

void mouse(int event, int x, int y, int flags, void* param) {
  drawEpiline(Point2f(x, y), 0, 4);
}

void display() {
  glBegin(GL_POINTS);
  glColor3ub(0, 255, 0);
  glVertex3f(0, 0, 0);
  for (int i = 0; i < frames.size(); i++) {
    glColor3ub(255, 0, 0);
    glVertex3f(frames[i].camera.at<double>(0, 0), 
               frames[i].camera.at<double>(1, 0), 
               frames[i].camera.at<double>(2, 0));
  }
  glEnd();
}

void loadDataset() {
  FILE *fi = fopen((FLAGS_dataset + "/templeR_par.txt").c_str(), "r");
  printf("%d\n", fi);
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
    for (int k = 0; k < 9; k++) 
      fscanf(fi, " %lf", &frames[i].k.at<double>(0, 0) + k);
    for (int k = 0; k < 9; k++) 
      fscanf(fi, " %lf", &frames[i].r.at<double>(0, 0) + k);
    for (int k = 0; k < 3; k++) 
      fscanf(fi, " %lf", &frames[i].t.at<double>(0, 0) + k);

    frames[i].camera = - frames[i].r.t() * frames[i].t;
    //cout << frames[i].k << endl;
    //imshow("fram", frames[i].img);
    //waitKey();
  }
  fclose(fi);
  imshow("img0", frames[0].img);
  

  Viz::setDisplayCallback(display);
  Viz::startWindow(800, 800);

  //detectHaris(frames[0].img);
  vector<KeyPoint> kp[2];
  detectHarris(frames[0].img, kp[0]);
  displayKeyPoint("f0", frames[0].img, kp[0]);

  detectHarris(frames[1].img, kp[1]);
  displayKeyPoint("f1", frames[1].img, kp[1]);

  waitKey();

  setMouseCallback("img0", mouse);
  waitKey();


}
int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);

  loadDataset();
}
