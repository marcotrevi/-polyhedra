class vertex3D {
  int ID = 0;
  float x, y, z;
  int[] parents = new int[n_faces]; // parent faces

  vertex3D(float x, float y, float z) {
    this.x=x;
    this.y=y;
    this.z=z;
    for (int i=0; i<n_faces; i++) {
      parents[i] = 0;
    }
  }
  float x() {
    return x;
  }
  float y() {
    return y;
  }
  float z() {
    return z;
  }

  void printme() {
    println("vertex "+ID);
    print("parent faces: ");
    for (int i=0; i<n_faces; i++) {
      if (parents[i] == 1) {
        print(i+", ");
      }
    }
    println();
    print("neighbors: ");
    for (int i=0; i<n_vertices; i++) {
      if (Adj_vertices[ID][i] > 0 && i != ID) {
        print(i+", ");
      }
    }
    println();
    println();
  }
}