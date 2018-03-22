# gxbert
Modularization and Vectorization of the Geant4 Bertini Cascade Model

### Checkout gxbert

    cd ${SRC_DIR}
    git clone https://github.com/gxbert/gxbert.git 

### how to build tests
- gxbert

  cd ${BUILD_DIR}
  cmake $SRC_DIR/gxbert -DCMAKE_BUILD_TYPE=Release 
  make 
  ./gxbertTest proton Pb 1.5GeV 100

- g4bert

  cd ${BUILD_DIR}/g4bert
  cmake $SRC_DIR/gxbert/test/g4bert -DGeant4_DIR=/path_to_geant4_insall_dir
  make
  \#setup Geant4 envs
  ./g4bert proton Pb 1.5GeV 100
