#.. _Testing:

Testing
-------

Testing of `MARBLES` can be performed using CTest, which is included in the CMake build system. With CMake, this can be enabled with the following options ::

  $ cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
          -DCMAKE_CXX_COMPILER:STRING=clang++ \
          -DCMAKE_C_COMPILER:STRING=clang \
          -DCMAKE_BUILD_TYPE:STRING=Release \
          -DMARBLES_DIM:STRING=3 \
          -DMARBLES_ENABLE_MPI:BOOL=OFF \
          -DMARBLES_ENABLE_CPPCHECK:BOOL=OFF \
          -DMARBLES_ENABLE_CLANG_TIDY:BOOL=OFF \
          ..

To perform the tests and compare to previously generated gold files, use the following additional options: `-DMARBLES_TEST_WITH_FCOMPARE:BOOL=ON` and `-DMARBLES_REFERENCE_GOLDS_DIRECTORY:STRING=${HOME}/marbles/Build/golds/current`.
