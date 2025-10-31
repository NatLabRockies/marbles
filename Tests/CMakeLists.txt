configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

if(MARBLES_TEST_WITH_FCOMPARE)
  if("${MARBLES_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, MARBLES_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    set(GOLD_FILES_DIRECTORY ${MARBLES_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
    message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
  endif()
endif()

if(MARBLES_SAVE_GOLDS)
  if("${MARBLES_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, MARBLES_SAVED_GOLDS_DIRECTORY must be set and the directory exist")
  else()
    if(EXISTS ${MARBLES_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${MARBLES_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${MARBLES_SAVED_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
macro(setup_test)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    set(RUNTIME_OPTIONS "max_step=10 amr.plot_file=plt amr.plot_int=10 amr.chk_int=-1 amrex.the_arena_is_managed=0 amrex.abort_on_out_of_gpu_memory=1")
    if(MARBLES_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(MARBLES_ENABLE_MPI)
      if(MARBLES_ENABLE_CUDA)
        set(TEST_NP 2)
      else()
        set(TEST_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(TEST_NP 1)
      unset(MPI_COMMANDS)
    endif()
    if(MARBLES_ENABLE_CUDA OR (MARBLES_ENABLE_HYPRE AND (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")))
      set(FCOMPARE_TOLERANCE "-r 1e-10 --abs_tol 1.0e-12")
    endif()
    if(MARBLES_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME}/")
    endif()
    if(MARBLES_TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE_EXE} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME)
    setup_test()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${marbles_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 5400
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "regression"
                         ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME)
    add_test_r(${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no_ci")
endfunction(add_test_re)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(channel_regular)
add_test_r(tg)
add_test_r(single_rotated_box)

if(MARBLES_DIM EQUAL 2)
  add_test_r(channel_cylinder_amr_2d)
endif()

if(MARBLES_DIM EQUAL 3)
  add_test_r(single_cylinder)
  add_test_r(channel_cylinder)
  add_test_r(cylinder_turek_2d2)
endif()

add_test_re(sod)
add_test_re(sod_amr)

#=============================================================================
# Regression tests excluded from CI
#=============================================================================
#add_test_re(something)
