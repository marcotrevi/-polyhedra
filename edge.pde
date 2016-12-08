class edge {
  int ID = 0;
  int v0_ID, v1_ID;
  int F0_ID = -1;
  int F1_ID = -1;

  boolean locked = false;

  edge(int v0_ID, int v1_ID) {
    this.v0_ID = v0_ID;
    this.v1_ID = v1_ID;
  }
  
  void printme(){
    println("edge "+ID);
    println("endpoints: "+v0_ID+", "+v1_ID);
    println("faces: "+F0_ID+", "+F1_ID);
    println();
  }
}