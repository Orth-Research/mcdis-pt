/*
 class Boundaries.h defines arrays that include (nearest, face-diagonal, etc) neighboring sites of a site i. 
*/

class Boundaries {
 public:
  Boundaries(const int ddim, const int LL, const int LLz, int *nneighbors_direct, int *nneighbors_face, int *nneighbors_body, int *nneighbors_direct_nn);

 public:

 private: 
  // parameters passed to class
  int dim; 
  int L;
  int Lz;
  int *neighbors_direct;
  int *neighbors_face;
  int *neighbors_body;
  int *neighbors_direct_nn;
  // private class variables initialized at construction
  int N2D;
  int N3D;
};

Boundaries::Boundaries(const int ddim, const int LL, const int LLz, int *nneighbors_direct, int *nneighbors_face, int *nneighbors_body, int *nneighbors_direct_nn) : dim(ddim), L(LL), Lz(LLz), neighbors_direct(nneighbors_direct), neighbors_face(nneighbors_face), neighbors_body(nneighbors_body), neighbors_direct_nn(nneighbors_direct_nn), N2D(LL*LL), N3D(LL*LL*LLz)
{
  if ((dim == 2) || (Lz == 1)) {
    std::cout << "dim = " << dim << ", N2D = " << N2D << ", Lz = " << Lz << std::endl;
    // nearest-neighbors on the 2D square lattice
    for (int i=0; i < N2D; i++) {
      neighbors_direct[4*i] = (i + 1) % N2D;
      neighbors_direct[4*i+1] = (i + N2D - 1) % N2D;
      neighbors_direct[4*i+2] = (i + L) % N2D;
      neighbors_direct[4*i+3] = (i + N2D - L) % N2D;
      if (i%L==0) neighbors_direct[4*i+1] = i + L -1;
      if ((i+1)%L==0) neighbors_direct[4*i] = i - L + 1;

// for now set all further neighbor sites to zero (can be fixed later if needed)
      neighbors_face[4*i] = 0;
      neighbors_face[4*i+1] = 0;
      neighbors_face[4*i+2] = 0;
      neighbors_face[4*i+3] = 0;

      neighbors_body[4*i] = 0;
      neighbors_body[4*i+1] = 0;
      neighbors_body[4*i+2] = 0;
      neighbors_body[4*i+3] = 0;

      neighbors_direct_nn[4*i] = 0;
      neighbors_direct_nn[4*i+1] = 0;
      neighbors_direct_nn[4*i+2] = 0;
      neighbors_direct_nn[4*i+3] = 0;
    }
  } // if dim ==2

  else if (dim == 3) {
    // nearest-neighbors on the 3D cubic lattice with PBC along all three directions
    // fill in nearest-neighbor array
    for (int lz = 0; lz < Lz; lz++) {
      for (int i = 0; i < N2D; i++) {
        neighbors_direct[lz*6*N2D + 6*i] = lz*N2D + (i + 1) % N2D;
        neighbors_direct[lz*6*N2D + 6*i + 1] = lz*N2D + (i + N2D - 1) % N2D;
        neighbors_direct[lz*6*N2D + 6*i + 2] = lz*N2D + (i + L) % N2D;
        neighbors_direct[lz*6*N2D + 6*i + 3] = lz*N2D + (i + N2D - L) % N2D;
        neighbors_direct[lz*6*N2D + 6*i + 4] = ((lz + 1)%Lz)*N2D + i;
        neighbors_direct[lz*6*N2D + 6*i + 5] = ((lz + Lz - 1)%Lz)*N2D + i;

        if (i%L == 0) {
          neighbors_direct[lz*6*N2D + 6*i + 1] = lz*N2D + (i + N2D - 1 + L)%N2D;
        }
        if ((i+1)%L == 0) {
          neighbors_direct[lz*6*N2D + 6*i] = lz*N2D + (i + N2D - L + 1)%N2D;
        }         
      } // for i running within one layer
    } // for lz

    // face diagonal neighbors of 3D cubic lattice with PBC along all three directions
    for (int lz = 0; lz < Lz; lz++) {
      for (int i = 0; i < N2D; i++) {
        // same lz plane
        neighbors_face[lz*12*N2D + 12*i] = lz*N2D + (i + L + 1) % N2D; // right up
        neighbors_face[lz*12*N2D + 12*i + 1] = lz*N2D + (i + L - 1) % N2D; // left up
        neighbors_face[lz*12*N2D + 12*i + 2] = lz*N2D + (i + N2D - L + 1) % N2D; // right down
        neighbors_face[lz*12*N2D + 12*i + 3] = lz*N2D + (i + N2D - L - 1) % N2D; // left down

        // lz up
        neighbors_face[lz*12*N2D + 12*i + 4] = ((lz + 1)%Lz)*N2D + (i + 1) % N2D; // right 
        neighbors_face[lz*12*N2D + 12*i + 5] = ((lz + 1)%Lz)*N2D + (i + N2D - 1) % N2D; // left
        neighbors_face[lz*12*N2D + 12*i + 6] = ((lz + 1)%Lz)*N2D + (i + L) % N2D; // up
        neighbors_face[lz*12*N2D + 12*i + 7] = ((lz + 1)%Lz)*N2D + (i + N2D - L) % N2D; // down
        
        // lz down 
        neighbors_face[lz*12*N2D + 12*i + 8] = ((lz + Lz - 1)%Lz)*N2D + (i + 1) % N2D; // right
        neighbors_face[lz*12*N2D + 12*i + 9] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - 1) % N2D; // left
        neighbors_face[lz*12*N2D + 12*i + 10] = ((lz + Lz - 1)%Lz)*N2D + (i + L) % N2D; // up
        neighbors_face[lz*12*N2D + 12*i + 11] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - L) % N2D; // down

        if (i%L == 0) { // + L
          neighbors_face[lz*12*N2D + 12*i + 1] = lz*N2D + (i + L - 1 + L) % N2D; // left up
          neighbors_face[lz*12*N2D + 12*i + 3] = lz*N2D + (i + N2D - 1) % N2D; // left down
          neighbors_face[lz*12*N2D + 12*i + 5] = ((lz + 1)%Lz)*N2D + (i + N2D - 1 + L) % N2D; // left
          neighbors_face[lz*12*N2D + 12*i + 9] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - 1 + L) % N2D; // left
        }
        if ((i+1)%L == 0) { // - L
          neighbors_face[lz*12*N2D + 12*i] = lz*N2D + (i + 1) % N2D; // right up
          neighbors_face[lz*12*N2D + 12*i + 2] = lz*N2D + (i + N2D - L + 1 - L) % N2D; // right down
          neighbors_face[lz*12*N2D + 12*i + 4] = ((lz + 1)%Lz)*N2D + (i + 1 - L) % N2D; // right 
          neighbors_face[lz*12*N2D + 12*i + 8] = ((lz + Lz - 1)%Lz)*N2D + (i + 1 - L) % N2D; // right
        }

      }
    }

    // body diagonal neighbors of 3D cubic lattice with PBC along all three directions
    for (int lz = 0; lz < Lz; lz++) {
      for (int i = 0; i < N2D; i++) {
        //  lz up
        neighbors_body[lz*8*N2D + 8*i] = ((lz+1)%Lz)*N2D + (i + L + 1) % N2D; // right up
        neighbors_body[lz*8*N2D + 8*i + 1] = ((lz+1)%Lz)*N2D + (i + L - 1) % N2D; // left up
        neighbors_body[lz*8*N2D + 8*i + 2] = ((lz+1)%Lz)*N2D + (i + N2D - L + 1) % N2D; // right down
        neighbors_body[lz*8*N2D + 8*i + 3] = ((lz+1)%Lz)*N2D + (i + N2D - L - 1) % N2D; // left down
        // lz down
        neighbors_body[lz*8*N2D + 8*i + 4] = ((lz + Lz -1)%Lz)*N2D + (i + L + 1) % N2D; // right up
        neighbors_body[lz*8*N2D + 8*i + 5] = ((lz + Lz - 1)%Lz)*N2D + (i + L - 1) % N2D; // left up
        neighbors_body[lz*8*N2D + 8*i + 6] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - L + 1) % N2D; // right down
        neighbors_body[lz*8*N2D + 8*i + 7] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - L - 1) % N2D; // left down

        if (i%L == 0) { // + L 
          neighbors_body[lz*8*N2D + 8*i + 1] = ((lz+1)%Lz)*N2D + (i + L - 1 + L) % N2D; // left up
          neighbors_body[lz*8*N2D + 8*i + 3] = ((lz+1)%Lz)*N2D + (i + N2D - 1) % N2D; // left down
          neighbors_body[lz*8*N2D + 8*i + 5] = ((lz + Lz - 1)%Lz)*N2D + (i + L - 1 + L) % N2D; // left up
          neighbors_body[lz*8*N2D + 8*i + 7] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - 1) % N2D; // left down
        }

        if ((i+1)%L == 0) { // - L
          neighbors_body[lz*8*N2D + 8*i] = ((lz+1)%Lz)*N2D + (i + 1) % N2D; // right up
          neighbors_body[lz*8*N2D + 8*i + 2] = ((lz+1)%Lz)*N2D + (i + N2D - L + 1 - L) % N2D; // right down
          neighbors_body[lz*8*N2D + 8*i + 4] = ((lz + Lz -1)%Lz)*N2D + (i + 1) % N2D; // right up
          neighbors_body[lz*8*N2D + 8*i + 6] = ((lz + Lz - 1)%Lz)*N2D + (i + N2D - L + 1 - L) % N2D; // right down
        }

      }
    }

    // 2nd direct neighbors on the 3D cubic lattice with PBC along all three directions
    for (int lz = 0; lz < Lz; lz++) {
      for (int i = 0; i < N2D; i++) {
        neighbors_direct_nn[lz*6*N2D + 6*i] = lz*N2D + (i + 2) % N2D;
        neighbors_direct_nn[lz*6*N2D + 6*i + 1] = lz*N2D + (i + N2D - 2) % N2D;
        neighbors_direct_nn[lz*6*N2D + 6*i + 2] = lz*N2D + (i + 2*L) % N2D;
        neighbors_direct_nn[lz*6*N2D + 6*i + 3] = lz*N2D + (i + N2D - 2*L) % N2D;
        neighbors_direct_nn[lz*6*N2D + 6*i + 4] = ((lz + 2)%Lz)*N2D + i;
        neighbors_direct_nn[lz*6*N2D + 6*i + 5] = ((lz + Lz - 2)%Lz)*N2D + i;
        if (i%L == 0) { // + L
          neighbors_direct_nn[lz*6*N2D + 6*i + 1] = lz*N2D + (i + N2D - 2 + L)%N2D;
        }
        if ((i+1)%L == 0) { // - L
          neighbors_direct_nn[lz*6*N2D + 6*i] = lz*N2D + (i + 2 - L + N2D)%N2D;
        }         
      } // for i running within one layer
    } // for lz

  } // if dim == 3
  else {
    std::cout << "You have specified dimension different from 2 or 3. This is not supported by this program!" << std::endl;
  }

} // closes class constructor
