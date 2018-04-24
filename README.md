# gxbert
Modularization and Vectorization of the Geant4 Bertini Cascade Model

### Checkout gxbert

    cd ${SRC_DIR}
    git clone https://github.com/gxbert/gxbert.git 

### how to build and test gxbert (gxbertTest)

    cd ${BUILD_DIR}
    cmake $SRC_DIR/gxbert -DCMAKE_BUILD_TYPE=Release 
    make 
    ./gxbertTest proton Pb 1.5GeV 100

### how to build and test g4bert

    cd ${BUILD_DIR}/g4bert
    cmake $SRC_DIR/gxbert/test/g4bert -DGeant4_DIR=/path_to_geant4_insall_dir
    make
    ./g4bert proton Pb 1.5GeV 100

### how to build gvbertTest (the main routine to test SIMD vectorization of a 
    simplified chain of gxbert)

    export VCROOT=${INSTALL_DIR}/vc_1_3_0
    export VECCORE=${INSTALL_DIR}/VecCore

    cd ${BUILD_DIR}
    cmake $SRC_DIR/gxbert -DCMAKE_BUILD_TYPE=Release -DGXBERT_VECTEST=ON \
          -DBACKEND=Vc   -DCMAKE_PREFIX_PATH="${VCROOT};${VECCORE}"
    make 

### how to build benchmark (unit tests for SIMD vectorization and GPU)

    cd ${BUILD_DIR}
    cmake $SRC_DIR/gxbert -DCMAKE_BUILD_TYPE=Release -DGXBERT_VECTEST=ON \
      -DBACKEND=Vc -DCMAKE_PREFIX_PATH="${VCROOT};${VECCORE}" -DDGXBERT_CUDA=ON
    make 
    ./benchmark [ntracks] [nrepetitions] [minEnergy] [maxEnergy]
