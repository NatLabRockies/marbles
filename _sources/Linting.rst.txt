#.. _Linting:

Linting
-------

Linting of `MARBLES` can be performed using cmake options to configure cppcheck, clang-tidy, and codespell. For example, the continuous integration tools through Github Actions uses the following for clang-tidy ::

  $ cmake -DCMAKE_BUILD_TYPE:STRING=Debug \
          -DCMAKE_CXX_COMPILER:STRING=clang++ \
          -DCMAKE_C_COMPILER:STRING=clang \
          -DMARBLES_DIM:STRING=3 \
          -DMARBLES_ENABLE_MPI:BOOL=OFF \
          -DMARBLES_TEST_WITH_FCOMPARE:BOOL=OFF \
          -DMARBLES_ENABLE_ALL_WARNINGS:BOOL=ON \
          -DMARBLES_ENABLE_CPPCHECK:BOOL=OFF \
          -DMARBLES_ENABLE_CLANG_TIDY:BOOL=ON \
          ..
  $ cmake --build . | tee -a clang-tidy-full-report.txt
