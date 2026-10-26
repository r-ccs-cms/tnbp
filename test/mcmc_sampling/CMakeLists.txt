# CUDA executables are launched explicitly inside a GPU allocation.
add_executable(mcmc_sampling ../../apps/mcmc_sampling/main.cc)
target_link_libraries(mcmc_sampling PRIVATE tnbp_test_backend)
add_executable(mcmc_sampling_complex ../../apps/mcmc_sampling/main.cc)
target_compile_definitions(mcmc_sampling_complex PRIVATE USE_COMPLEX)
target_link_libraries(mcmc_sampling_complex PRIVATE tnbp_test_backend)
if(TNBP_TEST_BACKEND STREQUAL "cuda")
    target_compile_definitions(mcmc_sampling PRIVATE USE_TCAPI_CUDA)
    target_compile_definitions(mcmc_sampling_complex PRIVATE USE_TCAPI_CUDA)
endif()
add_subdirectory(heavyhex51)
add_executable(mcmc_regression main.cc)
target_link_libraries(mcmc_regression PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    add_executable(mcmc_options_regression options_regression.cc)
    target_compile_features(mcmc_options_regression PRIVATE cxx_std_17)
    add_test(NAME mcmc_options_regression COMMAND mcmc_options_regression)
    foreach(kind IN ITEMS real complex)
        add_executable(mcmc_checkpoint_fixture_${kind} checkpoint_fixture.cc)
        target_link_libraries(mcmc_checkpoint_fixture_${kind} PRIVATE tnbp_test_backend)
        if(kind STREQUAL "complex")
            target_compile_definitions(mcmc_checkpoint_fixture_${kind} PRIVATE USE_COMPLEX)
            set(sampler mcmc_sampling_complex)
        else()
            set(sampler mcmc_sampling)
        endif()
        foreach(ranks IN ITEMS 1 2)
            set(fixture checkpoint_${kind}_${ranks})
            set(data_dir "${CMAKE_CURRENT_BINARY_DIR}/${fixture}")
            add_test(NAME ${fixture}_setup COMMAND mcmc_checkpoint_fixture_${kind} "${data_dir}" ${ranks})
            set_tests_properties(${fixture}_setup PROPERTIES FIXTURES_SETUP ${fixture})
            add_test(NAME mcmc_${fixture} COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:${sampler}> ${MPIEXEC_POSTFLAGS}
                --inputedge "${data_dir}/graph.edges" --linefile "${data_dir}/graph.lines"
                --inputname "${data_dir}/state" --num_samples 8 --num_batches 2
                --sampling_method 1 --bmps_bond_dim 16)
            set_tests_properties(mcmc_${fixture} PROPERTIES FIXTURES_REQUIRED ${fixture}
                TIMEOUT 30 PROCESSORS ${ranks} PASS_REGULAR_EXPRESSION "sample 7.1: 002"
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endforeach()
    foreach(ranks IN ITEMS 1 2)
        foreach(help_target IN ITEMS mcmc_sampling mcmc_sampling_complex)
            if(help_target STREQUAL "mcmc_sampling")
                set(help_flag -h)
            else()
                set(help_flag --help)
            endif()
            add_test(NAME ${help_target}_help_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:${help_target}> ${MPIEXEC_POSTFLAGS}
                --inputedge does-not-exist.edges ${help_flag})
            set_tests_properties(${help_target}_help_${ranks}_rank PROPERTIES
                TIMEOUT 30 PROCESSORS ${ranks}
                PASS_REGULAR_EXPRESSION "Usage: mcmc_sampling"
                FAIL_REGULAR_EXPRESSION "start program;sampling:;sample 0"
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
        add_test(NAME mcmc_heat_bath_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} heat-bath)
        set_tests_properties(mcmc_heat_bath_${ranks}_rank PROPERTIES TIMEOUT 180 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_cli_heat_bath_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_sampling_complex> ${MPIEXEC_POSTFLAGS}
            --inputedge "${CMAKE_CURRENT_SOURCE_DIR}/data/diagonal03.edges"
            --linefile "${CMAKE_CURRENT_SOURCE_DIR}/data/two.lines"
            --groupfile "${CMAKE_CURRENT_SOURCE_DIR}/data/paired.groups"
            --p_bond_dim 2 --v_bond_dim 2 --bmps_bond_dim 16 --num_samples 16 --sampling_method 1)
        set_tests_properties(mcmc_cli_heat_bath_${ranks}_rank PROPERTIES TIMEOUT 30 PROCESSORS ${ranks}
            PASS_REGULAR_EXPRESSION "sample 15.0: [01][01][01][01]"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_automatic_groups_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} automatic-groups)
        set_tests_properties(mcmc_automatic_groups_${ranks}_rank PROPERTIES TIMEOUT 180 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_shift_absorption_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} shift-absorption)
        set_tests_properties(mcmc_shift_absorption_${ranks}_rank PROPERTIES TIMEOUT 180 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_fixed_sites_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} fixed-sites)
        set_tests_properties(mcmc_fixed_sites_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_cli_fixed_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_sampling_complex> ${MPIEXEC_POSTFLAGS}
            --inputedge "${CMAKE_CURRENT_SOURCE_DIR}/data/diagonal03.edges"
            --linefile "${CMAKE_CURRENT_SOURCE_DIR}/data/two.lines"
            --groupfile "${CMAKE_CURRENT_SOURCE_DIR}/data/paired.groups"
            --num_samples 4)
        set_tests_properties(mcmc_cli_fixed_${ranks}_rank PROPERTIES TIMEOUT 30 PROCESSORS ${ranks}
            PASS_REGULAR_EXPRESSION "sample 3.0: 0000"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_singleton_ends_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} singleton-ends)
        set_tests_properties(mcmc_singleton_ends_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(mcmc_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_manual_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_regression> ${MPIEXEC_POSTFLAGS} manual)
        set_tests_properties(mcmc_manual_${ranks}_rank PROPERTIES TIMEOUT 120 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_cli_automatic_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_sampling_complex> ${MPIEXEC_POSTFLAGS}
            --inputedge "${CMAKE_CURRENT_SOURCE_DIR}/data/diagonal03.edges"
            --linefile "${CMAKE_CURRENT_SOURCE_DIR}/data/two.lines"
            --p_bond_dim 2 --v_bond_dim 2 --bmps_bond_dim 16
            --num_samples 16 --num_batches 2 --sampling_method 1)
        set_tests_properties(mcmc_cli_automatic_${ranks}_rank PROPERTIES TIMEOUT 30 PROCESSORS ${ranks}
            PASS_REGULAR_EXPRESSION "sample 15.1: [01][01][01][01]"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME mcmc_cli_manual_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:mcmc_sampling_complex> ${MPIEXEC_POSTFLAGS}
            --inputedge "${CMAKE_CURRENT_SOURCE_DIR}/data/diagonal03.edges"
            --linefile "${CMAKE_CURRENT_SOURCE_DIR}/data/two.lines"
            --groupfile "${CMAKE_CURRENT_SOURCE_DIR}/data/paired.groups"
            --p_bond_dim 2 --v_bond_dim 2 --bmps_bond_dim 16 --num_samples 16 --num_batches 2)
        set_tests_properties(mcmc_cli_manual_${ranks}_rank PROPERTIES TIMEOUT 30 PROCESSORS ${ranks}
            PASS_REGULAR_EXPRESSION "sample 15.1: [01][01][01][01]"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()

# Generate native CUDA checkpoints inside a one-GPU scheduler step.
if(TNBP_TEST_BACKEND STREQUAL "cuda")
    foreach(kind IN ITEMS real complex)
        add_executable(mcmc_checkpoint_fixture_${kind} checkpoint_fixture.cc)
        target_link_libraries(mcmc_checkpoint_fixture_${kind} PRIVATE tnbp_test_backend)
        target_compile_definitions(mcmc_checkpoint_fixture_${kind} PRIVATE USE_TCAPI_CUDA)
        if(kind STREQUAL "complex")
            target_compile_definitions(mcmc_checkpoint_fixture_${kind} PRIVATE USE_COMPLEX)
        endif()
    endforeach()
endif()

add_executable(mcmc_full_heat_bath full_heat_bath.cc)
target_link_libraries(mcmc_full_heat_bath PRIVATE tnbp_test_backend)
if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2)
        foreach(mode IN ITEMS kernel integration)
            if(mode STREQUAL "kernel")
                set(command_args $<TARGET_FILE:mcmc_full_heat_bath>)
            else()
                set(command_args $<TARGET_FILE:mcmc_regression> full-heat-bath)
            endif()
            add_test(NAME mcmc_full_heat_bath_${mode}_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS} ${command_args})
            set_tests_properties(mcmc_full_heat_bath_${mode}_${ranks}_rank PROPERTIES
                TIMEOUT 180 PROCESSORS ${ranks} ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
            if(mode STREQUAL "integration")
                # Eight 4096-sample distribution cases; accuracy thresholds stay unchanged.
                set_tests_properties(mcmc_full_heat_bath_${mode}_${ranks}_rank PROPERTIES TIMEOUT 300)
            endif()
        endforeach()
        foreach(width IN ITEMS 1 2)
            add_test(NAME mcmc_cli_full_heat_bath_${width}_sites_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:mcmc_sampling_complex> --inputedge "${CMAKE_CURRENT_SOURCE_DIR}/data/diagonal03.edges"
                --linefile "${CMAKE_CURRENT_SOURCE_DIR}/data/two.lines"
                --p_bond_dim 3 --v_bond_dim 2 --bmps_bond_dim 16 --num_samples 16 --num_batches 2
                --update_sites ${width})
            set_tests_properties(mcmc_cli_full_heat_bath_${width}_sites_${ranks}_rank PROPERTIES
                TIMEOUT 30 PROCESSORS ${ranks} PASS_REGULAR_EXPRESSION "sample 15.1: [012][012][012][012]"
                ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endforeach()
endif()

if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    add_test(NAME mcmc_full_heat_bath_kernel_3_rank COMMAND ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mcmc_full_heat_bath>)
    set_tests_properties(mcmc_full_heat_bath_kernel_3_rank PROPERTIES TIMEOUT 180 PROCESSORS 3
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()


if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME mcmc_honeycomb_kernel_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mcmc_full_heat_bath> honeycomb)
        set_tests_properties(mcmc_honeycomb_kernel_${ranks}_rank PROPERTIES TIMEOUT 180 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
    foreach(ranks IN ITEMS 1 2)
        add_test(NAME mcmc_honeycomb_integration_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mcmc_regression> honeycomb)
        set_tests_properties(mcmc_honeycomb_integration_${ranks}_rank PROPERTIES TIMEOUT 300 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()

if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME mcmc_honeycomb_zigzag_${ranks}_rank COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mcmc_full_heat_bath> zigzag)
        set_tests_properties(mcmc_honeycomb_zigzag_${ranks}_rank PROPERTIES TIMEOUT 180 PROCESSORS ${ranks}
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
    foreach(ranks IN ITEMS 1 2)
        foreach(layout IN ITEMS automatic paired)
            foreach(width IN ITEMS 1 2)
                add_test(NAME mcmc_zigzag_distribution_${layout}_${width}_sites_${ranks}_rank
                    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:mcmc_regression> zigzag-distribution ${layout} ${width})
                set_tests_properties(mcmc_zigzag_distribution_${layout}_${width}_sites_${ranks}_rank
                    PROPERTIES TIMEOUT 300 PROCESSORS ${ranks}
                    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
            endforeach()
        endforeach()
    endforeach()
endif()

if(BUILD_TESTING AND TNBP_TEST_BACKEND STREQUAL "cpu")
    foreach(ranks IN ITEMS 1 2 3)
        add_test(NAME mcmc_heat_bath_fatal_${ranks}_rank COMMAND ${CMAKE_COMMAND}
            -DMPIEXEC=${MPIEXEC_EXECUTABLE} -DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG} -DRANKS=${ranks}
            -DEXECUTABLE=$<TARGET_FILE:mcmc_full_heat_bath> -DCASES=0,1,2,3,4,5,6,7,8,9
            -P ${CMAKE_CURRENT_SOURCE_DIR}/../common/run_fatal_cases.cmake)
        set_tests_properties(mcmc_heat_bath_fatal_${ranks}_rank PROPERTIES TIMEOUT 90
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()
