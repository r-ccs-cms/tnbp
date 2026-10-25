# Generic graph contraction and bond relocation, with optional CUDA compilation.
add_executable(graph_contraction_regression regression.cc)
target_link_libraries(graph_contraction_regression PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME graph_contraction_regression_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:graph_contraction_regression> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(graph_contraction_regression_${ranks}_rank PROPERTIES TIMEOUT 90 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()

# Shape-only evaluator is standalone: no tensor backend or MPI dependency.
add_executable(graph_contraction_cost_shapes cost_shapes.cc)
target_compile_features(graph_contraction_cost_shapes PRIVATE cxx_std_17)
target_include_directories(graph_contraction_cost_shapes PRIVATE "${TNBP_ROOT}/include")
add_executable(graph_contraction_cost_tensor cost_tensor.cc)
target_link_libraries(graph_contraction_cost_tensor PRIVATE tnbp_test_backend)
if(BUILD_TESTING)
    add_test(NAME graph_contraction_cost_shapes COMMAND graph_contraction_cost_shapes)
    set_tests_properties(graph_contraction_cost_shapes PROPERTIES TIMEOUT 30)
    if(TNBP_TEST_BACKEND STREQUAL "cpu")
        foreach(ranks IN ITEMS 1 2 3)
            add_test(NAME graph_contraction_cost_tensor_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:graph_contraction_cost_tensor> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(graph_contraction_cost_tensor_${ranks}_rank PROPERTIES
                TIMEOUT 90 PROCESSORS ${ranks} ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endif()
endif()
