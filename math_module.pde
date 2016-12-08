class math_module {

  //************************************ DISPLAY VARIABLES

  int n_angles = 200;
  float focal = 0.5; // focus distance from camera plane.
  float distance = 10; // camera pinhole distance from origin.
  int theta_index_camera = 0;
  int phi_index_camera = 0;
  float zoom = 50; // zoom on views

  //************************************ GEOMETRY VARIABLES

  vertex3D _n_ = new vertex3D(1, 0, 0); // default unit vectors
  vertex3D _u_ = new vertex3D(0, 1, 0); 
  vertex3D _v_ = new vertex3D(0, 0, 1);
  float cos_phi, sin_phi, cos_theta, sin_theta;
  float[] sine_table = new float[n_angles];
  float[] cosine_table = new float[n_angles];

  //************************************ MATH VARIABLES

  // variables for LU decomposition
  int n = 3; // space dimension
  float[][] P = new float[n][n];
  float[][] L = new float[n][n];
  float[][] U = new float[n][n];
  boolean valid = true;
  float[][] R = new float[2][2];

  math_module() {
    for (int i=0; i<n_angles; i++) {
      sine_table[i] = sin(i*2*PI/n_angles);
      cosine_table[i] = cos(i*2*PI/n_angles);
    }
  }

  //************************************ DISPLAY METHODS

  void update_camera() {
    cos_phi = cosine_table[phi_index_camera];
    sin_phi = sine_table[phi_index_camera];
    cos_theta = cosine_table[theta_index_camera];
    sin_theta = sine_table[theta_index_camera];

    // update unit vectors

    _n_.x = cos_phi*cos_theta;
    _n_.y = cos_phi*sin_theta;
    _n_.z = sin_phi;

    _u_.x = -sin_theta;
    _u_.y = cos_theta;
    _u_.z = 0;

    _v_.x = -sin_phi*cos_theta;
    _v_.y = -sin_phi*sin_theta;
    _v_.z = cos_phi;
  }

  void rotate_camera() {
    if (theta_index_camera + 1 < n_angles) {
      theta_index_camera++;
    } else {
      theta_index_camera = 0;
    }
    update_camera();
  }

  vertex3D camera_center() {
    vertex3D CCenter = new vertex3D(0, 0, 0);
    CCenter.x = (distance)*cosine_table[phi_index_camera]*cosine_table[theta_index_camera];
    CCenter.y = (distance)*cosine_table[phi_index_camera]*sine_table[theta_index_camera];
    CCenter.z = (distance)*sine_table[phi_index_camera];
    return CCenter;
  }

  vertex3D vertex_point(vertex3D P) { // returns Q, point of vertex between point P and image plane
    // standard pinhole model
    vertex3D Q = new vertex3D(0, 0, 0);
    float h = scalar(P, _n_);
    float D = h - distance;
    float dh = distance*h; 
    float df = 2*(distance+focal);
    Q.x = (dh*_n_.x+df*P.x)/D;
    Q.y = (dh*_n_.y+df*P.y)/D;
    Q.z = (dh*_n_.z+df*P.z)/D;
    return Q;
  }

  float[] coords(vertex3D P) {
    vertex3D Q = vertex_point(P);
    float[] coordinates = new float[2];
    coordinates[0] = -scalar(Q, _u_);
    coordinates[1] = -scalar(Q, _v_);
    return coordinates;
  }

  //************************************ MATH METHODS

  float _norm(vertex3D u) {
    return sqrt(scalar(u, u));
  }

  float distance(float[] u, float[] v) {
    float d = 0;
    for (int i=0; i<n; i++) {
      d = d+(u[i]-v[i])*(u[i]-v[i]);
    }
    return sqrt(d);
  }

  float scalar(vertex3D u, vertex3D v) {
    return u.x*v.x + u.y*v.y + u.z*v.z;
  }

  vertex3D outer_product(vertex3D u, vertex3D v) {
    vertex3D w = new vertex3D(0, 0, 0);
    w.x = u.y*v.z - u.z*v.y;
    w.y = u.z*v.x - u.x*v.z;
    w.z = u.x*v.y - u.y*v.x;
    return w;
  }

  vertex3D matrixproduct(float[][] A, vertex3D v) {
    vertex3D w = new vertex3D(0, 0, 0);
    w.x = A[0][0]*v.x+A[0][1]*v.y+A[0][2]*v.z; 
    w.y = A[1][0]*v.x+A[1][1]*v.y+A[1][2]*v.z; 
    w.z = A[2][0]*v.x+A[2][1]*v.y+A[2][2]*v.z; 
    return w;
  }

  float[] prod_1(float[][] M, float[] v) {
    float[] u = new float[n];
    u[0]=0;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        u[i] = u[i]+M[i][j]*v[j];
      }
    }
    return u;
  }

  vertex3D intersect(face p1, face p2, face p3) {
    vertex3D P = new vertex3D(0, 0, 0);
    float[][] A = new float[3][3];
    float[] b = new float[3];
    float[] x = new float[3];
    A[0][0] = p1.a;
    A[0][1] = p1.b;
    A[0][2] = p1.c;
    A[1][0] = p2.a;
    A[1][1] = p2.b;
    A[1][2] = p2.c;
    A[2][0] = p3.a;
    A[2][1] = p3.b;
    A[2][2] = p3.c;
    b[0] = p1.k;
    b[1] = p2.k;
    b[2] = p3.k;
    x = solve(A, b);
    P.x = x[0];
    P.y = x[1];
    P.z = x[2];
    P.parents[p1.ID] = 1;
    P.parents[p2.ID] = 1;
    P.parents[p3.ID] = 1;
    if (!valid) {
      P.ID = -1;
    }
    return P;
  }

  float[] backsubstitution(float[][] U, float[] b) {
    // U has to be a nonsingular upper triangular matrix for this to work.
    // solves Ux = b.
    float[] x = new float[n];
    float s;
    x[n-1] = b[n-1]/U[n-1][n-1];
    for (int i=n-2; i>=0; i--) {
      s = 0; 
      for (int j=i+1; j<n; j++) {
        s = s+U[i][j]*x[j];
      }
      x[i] = (b[i]-s)/U[i][i];
    }
    return x;
  }

  float[] fwdsubstitution(float[][] L, float[] b) {
    // L has to be a nonsingular lower triangular matrix for this to work.
    // solves Lx = b.
    float[] x = new float[n];
    float s;
    x[0] = b[0]/L[0][0];
    for (int i=0; i<n; i++) {
      s = 0; 
      for (int j=0; j<i; j++) {
        s = s+L[i][j]*x[j];
      }
      x[i] = (b[i]-s)/L[i][i];
    }
    return x;
  }

  int maxindex(float[][] A, int j) {
    // returns position of max entry in column j of matrix A
    // beginning from element A[j][j] downwards.
    int index = j;
    float max = abs(A[j][j]);
    for (int i=j; i<n; i++) {
      if (abs(A[i][j])>max) {
        max = abs(A[i][j]);
        index = i;
      }
    }
    return index;
  }

  void swap(int i, int j, float[][] A) {
    // swaps row i with row j in matrix A.
    float a = 0;
    for (int k=0; k<n; k++) {
      a = A[i][k];
      A[i][k] = A[j][k];
      A[j][k] = a;
    }
  }

  float[][] identity(int n) {
    float[][] A = new float[n][n];
    for (int i=0; i<n; i++) {
      for (int j=i+1; j<n; j++) {
        A[i][j] = 0;
        A[j][i] = 0;
      }
      A[i][i] = 1;
    }
    return A;
  }

  float[] solve(float[][] A, float[] b) { // LU decompossition
    // solves Ax = b
    valid = true;
    float[] x = new float[n];
    float[] y = new float[n];
    // P,L,U initialization
    P = identity(n);
    L = identity(n);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        U[i][j] = A[i][j];
      }
    }
    // for each column find the pivot:
    int index = 0;
    float r = 0;
    for (int j=0; j<n; j++) {
      index = maxindex(U, j); // pivot
      if (index != j) {
        swap(j, index, U);
        swap(j, index, P);
        if (j>0) {
          // swapping subdiagonal entries of L
          float l = 0;
          for (int k=0; k<j; k++) {
            l = L[j][k];
            L[j][k] = L[index][k];
            L[index][k] = l;
          }
        }
      }
      // setting to zero all elements below U[j][j]
      for (int i=j+1; i<n; i++) {
        r = U[i][j]/U[j][j];
        L[i][j] = r;
        for (int k=j; k<n; k++) {
          U[i][k] = U[i][k] - r*U[j][k];
        }
      }
    }
    // Now we have the matrices P, L, U such that PA = LU.
    // Checking if L or U are singular.
    float d = L[0][0]*L[1][1]*L[2][2]*U[0][0]*U[1][1]*U[2][2];
    if (abs(d) < 0.00001 || Float.isNaN(d)) {
      valid = false;
    } else {
      // Solving first Ly = Pb
      y = fwdsubstitution(L, prod_1(P, b));
      // Now solving Ux = y
      x = backsubstitution(U, y);
    }
    return x;
  }

  float[] coordinates2D(int index, int k) {
    // returns coordinates of vertex with ID "k" in base (u,v) of face with ID "index" 
    float[] coords = new float[2];
    float x = vertices[k].x - faces[index].O.x;
    float y = vertices[k].y - faces[index].O.y;
    float z = vertices[k].z - faces[index].O.z;
    coords[0] = x*faces[index].u.x+y*faces[index].u.y+z*faces[index].u.z;
    coords[1] = x*faces[index].v.x+y*faces[index].v.y+z*faces[index].v.z;
    return coords;
  }

  void setparents(vertex3D P, vertex3D Q) {
    for (int i=0; i<n_faces; i++) {
      P.parents[i] = Q.parents[i];
    }
  }

  boolean AreConnected(int i, int j) {
    float r = 0;
    for (int k=0; k<n_faces; k++) {
      r = r+vertices[i].parents[k]*vertices[j].parents[k];
    }
    if (r > 1) {
      return true;
    } else {
      return false;
    }
  }

  //************************************ SETTING ADJACENCY MATRICES

  void setAdj_vertices() {
    Adj_vertices = new float[n_vertices][n_vertices];
    for (int i=0; i<n_vertices; i++) {
      for (int j=i; j<n_vertices; j++) {
        Adj_vertices[i][j] = 0;
        Adj_vertices[j][i] = 0;
      }
    }

    float[] P = new float[3];
    float[] Q = new float[3];
    int count = 0;
    for (int i=0; i<n_vertices; i++) {
      for (int j=i+1; j<n_vertices; j++) {
        if (AreConnected(i, j)) {
          P[0] = vertices[i].x; 
          P[1] = vertices[i].y; 
          P[2] = vertices[i].z; 
          Q[0] = vertices[j].x; 
          Q[1] = vertices[j].y; 
          Q[2] = vertices[j].z; 
          Adj_vertices[i][j] = distance(P, Q);
          Adj_vertices[j][i] = distance(P, Q);
          // populating edge list
          edges[count] = new edge(i, j);
          edges[count].ID = count;
          count = count + 1;
        }
      }
    }
  }

  void setAdj_faces() {
    Adj_faces = new float[n_faces][n_faces];
    for (int i=0; i<n_faces; i++) {
      for (int j=i; j<n_faces; j++) {
        Adj_faces[i][j] = 0;
        Adj_faces[j][i] = 0;
      }
    }

    for (int i=0; i<n_faces; i++) {
      for (int j=i+1; j<n_faces; j++) {
        int count = 0;
        for (int k=0; k<n_vertices; k++) {
          count = count + vertices[k].parents[i]*vertices[k].parents[j];// if point belongs to both planes
        }
        if (count == 2) {
          Adj_faces[i][j] = 1;
          Adj_faces[j][i] = 1;
        }
      }
    }
  }

  void MST() {
    // minimum spanning tree (Kruskal's way)
    // builds the MST's adjacency matrix.
    // creating ordered list of edges
    Adj_mst = new float[n_vertices][n_vertices];
    int n = n_edges;
    int[][] edgelist = new int[n][3];
    // edgelist[k][0] and edgelist[k][1] are the edge's nodes
    // edgelist[k][2] is the subtree label needed to calculate the MST
    int k=0; 
    for (int i=0; i<n_vertices; i++) {
      for (int j=i+1; j<n_vertices; j++) {
        if (Adj_vertices[i][j]>0) {
          edgelist[k][0] = i;
          edgelist[k][1] = j;
          edgelist[k][2] = 0; // defaul label
          k= k+1;
        }
      }
    }
    // sorting edgelist with StupidSort
    int temp0 = -1;
    int temp1 = -1;
    int i0, i1;
    int j0, j1;
    for (int i=0; i<n; i++) {
      i0 = edgelist[i][0];
      i1 = edgelist[i][1];
      for (int j=i+1; j<n; j++) {
        j0 = edgelist[j][0];
        j1 = edgelist[j][1];
        if (Adj_vertices[i0][i1] > Adj_vertices[j0][j1]) {
          // swap edges
          temp0 = edgelist[i][0];
          temp1 = edgelist[i][1];
          edgelist[i][0] = edgelist[j][0];
          edgelist[i][1] = edgelist[j][1];
          edgelist[j][0] = temp0;
          edgelist[j][1] = temp1;
          i0 = edgelist[i][0];
          i1 = edgelist[i][1];
          j0 = edgelist[j][0];
          j1 = edgelist[j][1];
        }
      }
    }
    // edgelist is now sorted
    // adding one edge at each step while checking that no cycle is formed
    int newlabel = 0; // subtree labels. MST will have label 1. 0 is not a valid label
    for (k=0; k<n; k++) {
      if (edgelist[k][2] == 0) {
        // found an unlabeled edge
        // cases:
        // 1) both nodes belong only to unlabeled edges => label edge with newlabel and update Adj_mst
        // 2) both nodes belong to different labeled edges => label edge and all edges of both trees with smallest label and update Adj_mst
        // 3) only one node belongs to a labeled edge => label edge with found label and update Adj_mst
        int v0_k = edgelist[k][0];
        int v1_k = edgelist[k][1];
        int v0_k_label = 0;
        int v1_k_label = 0;
        for (int i=0; i<n; i++) {
          if (edgelist[i][2] != 0) {
            // searching among labeled edges
            int v0_i = edgelist[i][0];
            int v1_i = edgelist[i][1];
            if (v0_k == v0_i || v0_k == v1_i) {
              v0_k_label = edgelist[i][2];
            }
            if (v1_k == v0_i || v1_k == v1_i) {
              v1_k_label = edgelist[i][2];
            }
          }
        }
        // case 1
        if (v0_k_label == 0 && v1_k_label == 0) {
          newlabel = newlabel + 1;
          edgelist[k][2] = newlabel;
          Adj_mst[v0_k][v1_k] = 1;
          Adj_mst[v1_k][v0_k] = 1;
        }
        // case 2
        if (v0_k_label !=0 && v1_k_label != 0 && v0_k_label != v1_k_label) {
          int minlabel, maxlabel;  
          if (v0_k_label < v1_k_label) {
            minlabel = v0_k_label;
            maxlabel = v1_k_label;
          } else {
            minlabel = v1_k_label;
            maxlabel = v0_k_label;
          }
          for (int i=0; i<n; i++) {
            if (edgelist[i][2] == maxlabel) {
              edgelist[i][2] = minlabel;
            }
          }
          edgelist[k][2] = minlabel;
          Adj_mst[v0_k][v1_k] = 1;
          Adj_mst[v1_k][v0_k] = 1;
        }
        // case 3
        if (v0_k_label*v1_k_label == 0 && v0_k_label + v1_k_label > 0) {
          int label;
          if (v0_k_label > 0) {
            label = v0_k_label;
          } else {
            label = v1_k_label;
          }
          edgelist[k][2] = label;
          Adj_mst[v0_k][v1_k] = 1;
          Adj_mst[v1_k][v0_k] = 1;
        }
      }
    }
  }


  void edgecount() {
    int c = 0;
    for (int i=0; i<n_vertices; i++) {
      for (int j=i+1; j<n_vertices; j++) {
        if (Adj_vertices[i][j] > 0) {
          c = c+1;
        }
      }
    }
    n_edges = c;
  }

  float[][] R(float[] u, float[] v) {
    // u and v are unit norm vectors
    // rotation matrix that sends u on v
    float[][] R = new float[2][2];
    R[0][0] = u[0]*v[0]+u[1]*v[1];
    R[1][1] = u[0]*v[0]+u[1]*v[1];
    R[0][1] = u[1]*v[0]-u[0]*v[1];
    R[1][0] = u[0]*v[1]-u[1]*v[0];
    return R;
  }

  void set_R(float[] C, float[] D, float[] A, float[] B) {
    // calculates rotation matrix R such as to send segment CD to segment AB.
    // the two segments must have the same length.
    float[] u = new float[2];
    float[] v = new float[2];

    u[0] = D[0]-C[0];
    u[1] = D[1]-C[1];
    float u_norm = sqrt(u[0]*u[0]+u[1]*u[1]);
    u[0] = u[0]/u_norm;
    u[1] = u[1]/u_norm;

    v[0] = B[0]-A[0];
    v[1] = B[1]-A[1];
    float v_norm = sqrt(v[0]*v[0]+v[1]*v[1]);
    v[0] = v[0]/v_norm;
    v[1] = v[1]/v_norm;

    R = R(u, v);
  }

  float[] transform(float[] X, float[] A, float[] X0, float[][] M) {
    // returns A+M.(X-X0)
    float[] Y = new float[2];
    Y[0] = X[0]-X0[0];
    Y[1] = X[1]-X0[1];
    float y0 = M[0][0]*Y[0]+M[0][1]*Y[1];
    float y1 = M[1][0]*Y[0]+M[1][1]*Y[1];
    Y[0] = A[0]+y0;
    Y[1] = A[1]+y1;
    return Y;
  }

  int[] get_edge(int F1, int F2) {
    // returns he IDs of the two vertices shared by face F1 and face F2
    int[] edge = new int[2];
    int count = 0;
    for (int i=0; i<n_vertices; i++) {
      if (vertices[i].parents[F1]*vertices[i].parents[F2] == 1) {
        edge[count] = i;
        count = count+1;
      }
    }
    return edge;
  }

  void display_stack(face F1, face F2) {
    // displays face F2 stacked on face F1 in 2D window.
    // faces must be neighbors.
    F1.display2D();
    int[] edge = new int[2];
    edge = get_edge(F1.ID, F2.ID);
    float[] A = coordinates2D(F1.ID, edge[0]);
    float[] B = coordinates2D(F1.ID, edge[1]);
    float[] C = coordinates2D(F2.ID, edge[0]);
    float[] D = coordinates2D(F2.ID, edge[1]);
    set_R(C, D, A, B);
    F2.display_transform2D(A, C, R);
  }

  void build_tree(face F) {
    F.marked = true;
    int[] edge = new int[2];
    for (int k=0; k<F.n_neighbors; k++) {
      int index = F.neighbors[k];
      edge = get_edge(F.ID, index);
      if (faces[index].marked == false && Adj_mst[edge[0]][edge[1]] != 1) {
        // if the face has not been visited yet and is not a "cut-face"

        int n_local_vertices = faces[index].n_local_vertices;
        int index0 = -1;
        int index1 = -1;
        for (int i=0; i<F.n_local_vertices; i++) {
          if (F.local_vertices[i] == edge[0]) {
            index0 = i;
          }
          if (F.local_vertices[i] == edge[1]) {
            index1 = i;
          }
        }

        // get new coordinates of the edge vertices of this face
        float[] A = new float[2];
        A = F.new_coords[index0];
        float[] B = new float[2];
        B = F.new_coords[index1];

        // get coordinates of the edge vertices in neighbor's base
        float[] C = new float[2];
        C = coordinates2D(index, edge[0]);
        float[] D = new float[2];
        D = coordinates2D(index, edge[1]);        

        // transforming coordinates
        set_R(C, D, A, B);
        for (int i=0; i<n_local_vertices; i++) {
          faces[index].new_coords[i] = math.transform(faces[index].coords[i], A, C, R);
        }

        build_tree(faces[index]);
      }
    }
  }
}