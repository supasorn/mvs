vector<Vec3d> _points;
vector<Vec3i> _faces;
void addCell(const Cell_handle &ch) {
  int base = _points.size() + 1;
  for (int i = 0; i < 4; i++) {
    auto v = ch->vertex(i)->point();
    _points.push_back(Vec3d(v.x(), v.y(), v.z()));
  }
  _faces.push_back(Vec3i(base, base + 2, base + 1));
  _faces.push_back(Vec3i(base + 1, base + 2, base + 3));
  _faces.push_back(Vec3i(base, base + 1, base + 3));
  _faces.push_back(Vec3i(base, base + 3, base + 2));
}

void addArrow(KPoint &a, KPoint &b) {
  int base = _points.size() + 1;
  _points.push_back(Vec3d(a.x(), a.y(), a.z()));
  _points.push_back(Vec3d(a.x() + 0.001, a.y(), a.z()));
  _points.push_back(Vec3d(b.x(), b.y(), b.z()));
  _points.push_back(Vec3d(b.x(), b.y(), b.z()));
  _faces.push_back(Vec3i(base, base+1, base+2));
}

void addFacet(DelaunayMesh::Facet &f, int flag) {
  int base = _points.size() + 1;
  for (int i = 1; i <= 3; i++) {
    auto v = f.first->vertex((f.second + i) % 4)->point();
    _points.push_back(Vec3d(v.x(), v.y(), v.z()));
  }
  if ((f.second % 2) ^ flag)
    _faces.push_back(Vec3i(base, base+2, base+1));
  else
    _faces.push_back(Vec3i(base, base+1, base+2));
}

void outObj() {
  FILE *fo = fopen("debug.obj", "w");
  for (int i = 0; i < _points.size(); i++) {
    fprintf(fo, "v %lf %lf %lf\n", _points[i][0], _points[i][1], _points[i][2]);
  }
  for (int i = 0; i < _faces.size(); i++) { 
    fprintf(fo, "f %d %d %d\n", _faces[i][0], _faces[i][1], _faces[i][2]);
  }
  fclose(fo);
}
