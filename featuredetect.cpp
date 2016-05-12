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

struct Feature {
  float r, c;
  int type;
};

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

void displayKeyPoint(Mat &im, vector<KeyPoint> &kp) {
  Mat tmp = im.clone();
  for (int i = 0; i < kp.size(); i++) {
    circle(tmp, kp[i].pt, 1, Scalar(0, 0, 255), -1);
  }
  imshow("keypoints", tmp);
  waitKey();
}

void mouse(int event, int x, int y, int flags, void* param) {
  drawEpiline(Point2f(x, y), 0, 4);
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

  //detectHaris(frames[0].img);
  vector<KeyPoint> kp;
  detectHarris(frames[0].img, kp);
  displayKeyPoint(frames[0].img, kp);

  setMouseCallback("img0", mouse);
  waitKey();


}
int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  google::InitGoogleLogging(argv[0]);

  loadDataset();
}
