class graphics_framework {
  boolean show_vertices = true; // toggle vertices display
  boolean show_normals = false; // toggle normal vectors display
  boolean show_edges = true; // toggle edge display
  boolean rotation = false; // toggle rotation
  boolean show_labels = false;
  boolean window3D = true;
  float zoom = 80;
  float xc = 0;
  float yc = 0;

  void connect(vertex3D A, vertex3D B, color c) {
    stroke(c);
    line(X(A), Y(A), X(B), Y(B));
  }

  float X(vertex3D P) {
    return W/2+math.coords(P)[0]*math.zoom;
  }

  float Y(vertex3D P) {
    return H/2-math.coords(P)[1]*math.zoom;
  }

  //*********************************************** DISPLAY METHODS

  void display_vertices() {
    noStroke();
    fill(255);
    for (int i=0; i<n_vertices; i++) {
      ellipse(X(vertices[i]), Y(vertices[i]), 3, 3);
      if (show_labels) {
        text(str(vertices[i].ID), X(vertices[i]), Y(vertices[i])-5);
      }
    }
  }

  void display_edges() {
    for (int i=0; i<n_vertices; i++) {
      for (int j=i+1; j<n_vertices; j++) {
        if (Adj_vertices[i][j]>0) {
          stroke(255);
          line(X(vertices[i]), Y(vertices[i]), X(vertices[j]), Y(vertices[j]));
        }
      }
    }
  }

  void display_faces() {
    for (int i=0; i<n_faces; i++) {
      faces[i].display();
    }
  }

  void display_faces2D() {
    for (int i=0; i<n_faces; i++) {
      faces[i].display2D();
      if (show_labels) {
        fill(0);
        text(str(faces[i].ID), W/2, H/2+50);
      }
    }
  }

  void display_diagram() {
    for (int i=0; i<n_faces; i++) {
      faces[i].transformed_display2D(zoom, xc, yc);
    }
  }

  void display_MST() {
    for (int i=0; i<n_vertices; i++) {
      for (int j=i+1; j<n_vertices; j++) {
        if (Adj_mst[i][j]>0) {
          stroke(255, 0, 0);
          line(X(vertices[i]), Y(vertices[i]), X(vertices[j]), Y(vertices[j]));
        }
      }
    }
  }

  //*********************************************** PRINT METHODS

  void printAdj_vertices() {
    println("vertices adjacency matrix:");
    println();
    for (int i=0; i<n_vertices; i++) {
      for (int j=0; j<n_vertices; j++) {
        print(round(Adj_vertices[i][j])+" ");
      }
      println();
    }
  }

  void printAdj_faces() {
    println("faces adjacency matrix:");
    println();
    for (int i=0; i<n_faces; i++) {
      for (int j=0; j<n_faces; j++) {
        print(round(Adj_faces[i][j])+" ");
      }
      println();
    }
  }

  void print_edges() {
    for (int i=0; i<n_edges; i++) {
      //edges[i].printme();
    }
  }

  void print_faces() {
    for (int i=0; i<n_faces; i++) {
      faces[i].printme();
    }
  }

  void print_vertices() {
    for (int i=0; i<n_vertices; i++) {
      vertices[i].printme();
    }
  }

  void print_info() {
    println("n. vertices: "+n_vertices);
    println("n. edges: "+n_edges);
    println("n. faces: "+n_faces);
    println();
  }

  //*********************************************** PRIMITIVES

  void cube() { // cube
    faces[0] = new face(0, 0, 1);
    faces[0].ID = 0;
    faces[1] = new face(PI/2, 0, 1);
    faces[1].ID = 1;
    faces[2] = new face(0, PI/2, 1);
    faces[2].ID = 2;
    faces[3] = new face(-PI/2, 0, 1);
    faces[3].ID = 3;
    faces[4] = new face(0, -PI/2, 1);
    faces[4].ID = 4;
    faces[5] = new face(PI, 0, 1);
    faces[5].ID = 5;
  }

  void prism(int n) { // regular prism
    // n+2 faces
    faces[0] = new face(0, -PI/2, 1);
    faces[0].ID = 0;
    for (int i=0; i<n; i++) {
      faces[i+1] = new face(2*PI/n*i, 0, 1);
      faces[i+1].ID = i+1;
    }
    faces[n+1] = new face(0, PI/2, 1);
    faces[n+1].ID = n+1;
  }

  void squad(int n) { // irregular prism
    // n+2 faces
    faces[0] = new face(0, -PI/2+random(-PI/10, PI/10), random(0.1, 3));
    faces[0].ID = 0;
    for (int i=0; i<n; i++) {
      faces[i+1] = new face(2*PI/n*i+random(-PI/10, PI/10), random(-PI/50, PI/50), random(0.1, 0.5));
      faces[i+1].ID = i+1;
    }
    faces[n+1] = new face(0, PI/2+random(-PI/10, PI/10), random(0.1, 3));
    faces[n+1].ID = n+1;
  }

  void dotted_line(float x1, float y1, float x2, float y2) {
    float ds = 5; // subline length
    float u1 = x2-x1;
    float u2 = y2-y1;
    float norm = sqrt(u1*u1+u2*u2);
    u1 = u1/norm;
    u2 = u2/norm;
    int n = floor(norm/(2*ds));
    println(n);
    for (int i=0; i<n; i++) {
      stroke(0);
      line(x1+2*i*u1*ds, y1+2*i*u2*ds, x1+(2*i+1)*u1*ds, y1+(2*i+1)*u2*ds);
      stroke(255);
      line(x1+(2*i+1)*u1*ds, y1+(2*i+1)*u2*ds, x1+(2*i+2)*u1*ds, y1+(2*i+2)*u2*ds);
    }
  }


  //*********************************************** EVENT HANDLERS

  void UserInput() {
    if (window3D) {
      if (key == CODED) {
        if (keyCode == DOWN) {
          if (math.phi_index_camera+1==math.n_angles) {
            math.phi_index_camera=0;
          } else {
            math.phi_index_camera = math.phi_index_camera+1;
          }
        }
        if (keyCode == UP) {
          if (math.phi_index_camera-1==-1) {
            math.phi_index_camera = math.n_angles-1;
          } else {
            math.phi_index_camera = math.phi_index_camera - 1;
          }
        }
        if (keyCode == LEFT) {
          if (math.theta_index_camera+1==math.n_angles) {
            math.theta_index_camera=0;
          } else {
            math.theta_index_camera = math.theta_index_camera+1;
          }
        }
        if (keyCode == RIGHT) {
          if (math.theta_index_camera-1==-1) {
            math.theta_index_camera=math.n_angles-1;
          } else {
            math.theta_index_camera = math.theta_index_camera-1;
          }
        }
      }
      if (key == 'z') {
        math.zoom = math.zoom + 10;
      }
      if (key == 'x') {
        math.zoom = math.zoom - 10;
      }
      if (key == 'r') {
        rotation = !rotation;
      }
      if (key == 'n') {
        show_normals = !show_normals;
      }
      if (key == 'l') {
        show_labels = !show_labels;
      }
    } else { // 2D window
      if (key == CODED) {
        if (keyCode == DOWN) {
          yc = yc - 10;
        }
        if (keyCode == UP) {
          yc = yc + 10;
        }
        if (keyCode == LEFT) {
          xc = xc - 10;
        }
        if (keyCode == RIGHT) {
          xc = xc + 10;
        }
      }
      if (key == 'z') {
        graphics.zoom = graphics.zoom + 10;
      }
      if (key == 'x') {
        graphics.zoom = graphics.zoom - 10;
      }
    }
    if (key == 'v') {
      print_vertices();
    }
    if (key == 'e') {
      print_edges();
    }
    if (key == 'f') {
      print_faces();
    }
    if (key == 'q') {
      generate(8);
      math.MST();
    }
    if (key == 't') {
      for (int i=0; i<n_faces; i++) {
        faces[i].marked = false;
      }
      math.build_tree(faces[0]);
    }
    if (key == 'w') {
      window3D = !window3D;
    }
    if (key == 'i') {
      print_info();
    }
    math.update_camera();
  }
}