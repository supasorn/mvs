#define UTIL_CERES
#include "utilities.h"
using namespace cv;
using namespace std;
DEFINE_string(dataset, "/home/supasorn/mvs/datasets/templeRing", "a");
int n;

struct Frame {
  string name;
  Mat k, r, t;
  Mat img;
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

  Scalar color(255, 0, 0);

  Mat a = frames[i0].img.clone();
  Mat b = frames[i1].img.clone();
  if (fabs(lines[0][0]) > fabs(lines[0][1])) {
    line(b,
        Point(-lines[0][2] / lines[0][0], 0),
        Point(-(lines[0][2]+lines[0][1]*a.rows) / lines[0][0] ,a.rows),
        color);
  } else {
    line(b,
        Point(0,-lines[0][2]/lines[0][1]),
        Point(a.cols,-(lines[0][2]+lines[0][0]*a.cols)/lines[0][1]),
        color);
  }
  circle(a, ps[0], 1, color, -1, CV_AA);
  imshow("img0", a);
  imshow("img1", b);
  //waitKey();
}

void mouse(int event, int x, int y, int flags, void* param) {
  drawEpiline(Point2f(x, y), 0, 20);
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

    //cout << frames[i].k << endl;
    //imshow("fram", frames[i].img);
    //waitKey();
  }
  fclose(fi);
  imshow("img0", frames[0].img);
  setMouseCallback("img0", mouse);
  waitKey();


}
int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);

  loadDataset();
}
