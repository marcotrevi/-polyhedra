// global variables //<>//
int W = 1800;
int H = 900;

float[][] Adj_vertices; // adjacency matrix for the vertices
float[][] Adj_faces; // adjacency matrix for the faces
float[][] Adj_edges; // adjacency matrix for the edges
float[][] Adj_mst; // adjacency matrix for the vertices of the minimum spanning tree

// graphics and math modules
math_module math = new math_module();
graphics_framework graphics = new graphics_framework();

// test
int n_faces = 8;
int n_vertices;
int n_edges;

face[] faces = new face[n_faces];
vertex3D[] vertices;
edge[] edges;

// main
void setup() {
  size(1800, 900);
  //  graphics.cube();
  graphics.squad(6); // n_faces: +2
  //  generate(0);
  math.MST();
}

void draw() {
  if (graphics.window3D) {
    background(0);
    if (graphics.rotation) {
      math.rotate_camera();
    }
    graphics.display_faces();
    graphics.display_edges();
    //graphics.display_MST();
    graphics.display_vertices();
  } else {
    background(255);
    stroke(0);
    graphics.display_diagram();
  }
}

// event handlers
void keyPressed() {

  if (key == '1') {
    int nn = floor(random(3, 6));
    n_faces = nn+2+5; // 5 more faces
    faces = new face[n_faces];
    graphics.squad(nn); // n_faces: +2
    generate(nn+2);
    math.MST();
    for (int i=0; i<n_faces; i++) {
      faces[i].marked = false;
    }
    math.build_tree(faces[0]);
  }
  graphics.UserInput();
}

/*
void recursive_display(face F) {
 // displays open view of the polyhedron.
 // Given a face F, find all faces connected to F according to the complement of the MST.
 // mark all visited faces.
 // repeat until all faces are marked.
 mark(F);
 display(F);
 for (int i=0; i<F.sisters; i++) {
 if (F.sister[i] is not marked) {
 recursive_display(F.sister[i]);
 }
 }
 }
 */


void generate(int index) { // generates an "index + n_faces" polyhedron

  for (int i=index; i<n_faces; i++) {
    faces[i] = new face(random(-PI, PI), random(-PI, PI), random(0.5, 1));
    faces[i].ID = i;
  }

  n_vertices = n_faces*(n_faces-1)*(n_faces-2)/6;
  vertex3D[] buffer_pts = new vertex3D[n_vertices]; 

  int count = 0;
  vertex3D bufferPoint = new vertex3D(0, 0, 0);
  for (int i=0; i<n_faces; i++) {
    for (int j=i+1; j<n_faces; j++) {
      for (int k=j+1; k<n_faces; k++) {
        bufferPoint = math.intersect(faces[i], faces[j], faces[k]);
        if (bufferPoint.ID != -1) {
          buffer_pts[count] = bufferPoint;
          buffer_pts[count].ID = count;
          math.setparents(buffer_pts[count], bufferPoint);
          count = count+1;
        }
      }
    }
  }
  // removing outer vertices (vertices that lie outside the "origin cell")
  int count2 = count;
  for (int i=0; i<count; i++) {
    float xi = buffer_pts[i].x;
    float yi = buffer_pts[i].y;
    float zi = buffer_pts[i].z;
    for (int j=0; j<n_faces; j++) {
      float aj = faces[j].a;
      float bj = faces[j].b;
      float cj = faces[j].c;
      float kj = faces[j].k;
      if ((aj*xi+bj*yi+cj*zi - kj)*kj > 0.00001) {
        buffer_pts[i].ID = -1; // flag for removal
        count2 = count2 - 1;
        j = n_faces; // exit cycle
      }
    }
  }

  n_vertices = count2;
  // populating vertices array
  int ii = 0;
  vertices = new vertex3D[n_vertices];
  for (int i=0; i < count; i++) {
    if (buffer_pts[i].ID != -1) {
      vertices[ii] = new vertex3D(buffer_pts[i].x, buffer_pts[i].y, buffer_pts[i].z);
      vertices[ii].ID = ii;
      math.setparents(vertices[ii], buffer_pts[i]);
      ii = ii+1;
    }
  }

  n_edges = n_vertices + n_faces;
  edges = new edge[n_edges]; 
  math.setAdj_vertices(); // populates also edge list
  math.setAdj_faces();
  Adj_edges = new float[n_edges][n_edges];
  for (int i=0; i<n_faces; i++) {
    faces[i].find_vertices();
    faces[i].find_neighbors();
  }
}

void displayFace(int k, int i1, int i2, int j1, int j2) {
  // displays face k such that:
  // - vertex with index i1 goes to vertex with index j1
  // - vertex with index i2 goes to vertex with index j2
  float[] u1 = math.coordinates2D(k, i1);
  float[] u2 = math.coordinates2D(k, i2);
  float[] v1 = math.coordinates2D(k, j1);
  float[] v2 = math.coordinates2D(k, j2);
  math.set_R(u1, u2, v1, v2);
}