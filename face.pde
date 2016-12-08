class face {
  int ID = 0;
  int[] local_vertices;
  int[] neighbors;
  int n_local_vertices;
  int n_neighbors;
  vertex3D n = new vertex3D(0, 0, 1);
  vertex3D u = new vertex3D(1, 0, 0);
  vertex3D v = new vertex3D(0, 1, 0);
  vertex3D O = new vertex3D(0, 0, 0);
  float[][] coords;
  float[][] new_coords; // for 2D tree diagram display
  float theta, phi, d, a, b, c, k;
  float cos_theta, sin_theta, cos_phi, sin_phi;
  int sides = 36; // number of sides of polygon
  color col;
  boolean marked = false;

  face(float theta, float phi, float d) { // normal vector to the face
    // equation ax + by + cz = k
    this.theta = theta;
    this.phi = phi;
    this.d = d;
    col = color(random(50, 255), random(50, 255), random(50, 255), 50);
    update();
  }

  void update() {
    // update reference frame and coefficients
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    n.x = cos_phi*cos_theta;
    n.y = cos_phi*sin_theta;
    n.z = sin_phi;

    u.x = -sin_theta;
    u.y = cos_theta;
    u.z = 0;

    v.x = -sin_phi*cos_theta;
    v.y = -sin_phi*sin_theta;
    v.z = cos_phi;

    O.x = d*n.x;
    O.y = d*n.y;
    O.z = d*n.z;

    a = n.x;
    b = n.y;
    c = n.z;
    k = O.x*n.x+O.y*n.y+O.z*n.z;
  }

  void printme() {
    println("face "+ID);
    print("vertices: ");
    for (int i=0; i<n_local_vertices; i++) {
      print(local_vertices[i]+", ");
    }
    println();
    println();
    print("all neighbors: ");
    println();
    int[] edge = new int[2];
    for (int i=0; i<n_neighbors; i++) {
      edge = math.get_edge(ID, neighbors[i]);
      println("face "+neighbors[i]+" on edge ("+edge[0]+", "+edge[1]+")");
    }
    println();
    print("non-cut neighbors: ");
    println();
    for (int i=0; i<n_neighbors; i++) {
      edge = math.get_edge(ID, neighbors[i]);
      if (Adj_mst[edge[0]][edge[1]] != 1) {
        println("face "+neighbors[i]+" on edge ("+edge[0]+", "+edge[1]+")");
      }
    }
    println();
  }

  void display_plane() {
    update();
    vertex3D[] vertices = new vertex3D[sides];
    float cost = 0;
    float sint = 0;
    float r = 2;
    for (int i=0; i<sides; i++) {
      cost = cos(2*PI/sides*i);
      sint = sin(2*PI/sides*i);
      vertices[i] = new vertex3D(d*n.x+u.x*cost*r+v.x*sint*r, d*n.y+u.y*cost*r+v.y*sint*r, d*n.z+u.z*cost*r+v.z*sint*r);
    }
    for (int i=0; i<sides-1; i++) {
      noStroke();
      fill(col);
      triangle(
        graphics.X(O), graphics.Y(O), 
        graphics.X(vertices[i]), graphics.Y(vertices[i]), 
        graphics.X(vertices[i+1]), graphics.Y(vertices[i+1]) 
        );
    }
    triangle(
      graphics.X(O), graphics.Y(O), 
      graphics.X(vertices[sides-1]), graphics.Y(vertices[sides-1]), 
      graphics.X(vertices[0]), graphics.Y(vertices[0]) 
      );
  }

  void display() {
    if (graphics.show_labels) {
      fill(0, 255, 0);
      text(str(ID), graphics.X(O), graphics.Y(O));
    }
    if (graphics.show_normals) {
      stroke(0, 255, 0);
      vertex3D normal = new vertex3D(O.x+n.x, O.y+n.y, O.z+n.z);
      line(graphics.X(O), graphics.Y(O), graphics.X(normal), graphics.Y(normal));
    }
  }

  void find_neighbors() {
    n_neighbors = n_local_vertices;
    neighbors = new int[n_neighbors];
    int index = 0;
    for (int i=0; i<n_faces; i++) {
      if (Adj_faces[ID][i] == 1 && ID != i) {
        neighbors[index] = i;
        index = index+1;
      }
    }
  }

  void find_vertices() {
    int[] temp = new int[n_vertices];
    int count = 0;
    for (int i=0; i<n_vertices; i++) {
      if (vertices[i].parents[ID] == 1) {
        temp[count] = i;
        count = count +1;
      }
    }
    n_local_vertices = count;
    local_vertices = new int[n_local_vertices];
    for (int i=0; i<n_local_vertices; i++) {
      local_vertices[i] = temp[i];
    }
    arrange_vertices();

    coords = new float[n_local_vertices][2];
    new_coords = new float[n_local_vertices][2];

    for (int i=0; i<n_local_vertices; i++) {
      coords[i] = math.coordinates2D(ID, local_vertices[i]);
      new_coords[i] = coords[i];
    }
  }

  void arrange_vertices() {
    // rearranging local_vertices to form a chain of vertices
    int index1, index2, count, temp; 
    boolean found;
    for (int i=0; i<n_local_vertices-1; i++) {
      index1 = local_vertices[i];
      index2 = local_vertices[i+1];
      if (Adj_vertices[index1][index2] == 0) {
        // vertices are not connected
        found = false;
        count = i+2;
        while (!found && count < n_local_vertices) {
          index2 = local_vertices[count];
          if (Adj_vertices[index1][index2] > 0) {
            // vertex found
            found = true;
            temp = local_vertices[i+1];
            local_vertices[i+1] = index2;
            local_vertices[count] = temp;
          } else {          
            count = count + 1;
          }
        }
      }
    }
  }

  void display_face() {
    int count = 0;
    for (int i=0; i<n_vertices; i++) {
      if (vertices[i].parents[ID] == 1) { 
        count = count +1;
      }
    }
    float[][] v = new float[count][2];
    int t = 0;
    for (int i=0; i<n_vertices; i++) {
      if (vertices[i].parents[ID] == 1) { 
        v[t][0] = graphics.X(vertices[i]);
        v[t][1] = graphics.Y(vertices[i]);
        t = t+1;
      }
    }
    for (int i=0; i<count-1; i++) {
      stroke(0);
      line(v[i][0], v[i][1], v[i+1][0], v[i+1][1]);
    }
    line(v[count-1][0], v[count-1][1], v[0][0], v[0][1]);
  }

  void display_transform2D(float[] A, float[] C, float[][] R) {
    // each point P gets transformed to A+R.(P-C)
    float[][] new_coords = new float[n_local_vertices][2];
    for (int i=0; i<n_local_vertices; i++) {
      new_coords[i] = math.transform(coords[i], A, C, R);
    }
    float r = 20;
    for (int i=0; i<n_local_vertices-1; i++) {
      line(W/2+r*new_coords[i][0], H/2-r*new_coords[i][1], W/2+r*new_coords[i+1][0], H/2-r*new_coords[i+1][1]);
    }
    line(W/2+r*new_coords[n_local_vertices-1][0], H/2-r*new_coords[n_local_vertices-1][1], W/2+r*new_coords[0][0], H/2-r*new_coords[0][1]);
  }

  void display2D() {
    float r = 20;
    for (int i=0; i<n_local_vertices-1; i++) {
      line(W/2+r*coords[i][0], H/2-r*coords[i][1], W/2+r*coords[i+1][0], H/2-r*coords[i+1][1]);
    }
    line(W/2+r*coords[n_local_vertices-1][0], H/2-r*coords[n_local_vertices-1][1], W/2+r*coords[0][0], H/2-r*coords[0][1]);
  }

  void transformed_display2D(float r, float xc, float yc) {
    for (int i=0; i<n_local_vertices-1; i++) {
      if (Adj_mst[local_vertices[i]][local_vertices[i+1]] != 1) {
        graphics.dotted_line(W/2+xc+r*new_coords[i][0], H/2-yc-r*new_coords[i][1], W/2+xc+r*new_coords[i+1][0], H/2-yc-r*new_coords[i+1][1]);
      } else {
        stroke(0);
        line(W/2+xc+r*new_coords[i][0], H/2-yc-r*new_coords[i][1], W/2+xc+r*new_coords[i+1][0], H/2-yc-r*new_coords[i+1][1]);
      }
      if (Adj_mst[local_vertices[n_local_vertices-1]][local_vertices[0]] != 1) {
        graphics.dotted_line(W/2+xc+r*new_coords[n_local_vertices-1][0], H/2-yc-r*new_coords[n_local_vertices-1][1], W/2+xc+r*new_coords[0][0], H/2-yc-r*new_coords[0][1]);
      } else {
        stroke(0);
        line(W/2+xc+r*new_coords[n_local_vertices-1][0], H/2-yc-r*new_coords[n_local_vertices-1][1], W/2+xc+r*new_coords[0][0], H/2-yc-r*new_coords[0][1]);
      }
    }
  }
}