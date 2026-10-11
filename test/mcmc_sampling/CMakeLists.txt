# Sampling migration is currently validated on the CPU backend.
if(TNBP_TEST_BACKEND STREQUAL "cpu")
    add_executable(mcmc_sampling ../../apps/mcmc_sampling/main.cc)
    target_link_libraries(mcmc_sampling PRIVATE tnbp_test_backend)
    add_executable(mcmc_sampling_complex ../../apps/mcmc_sampling/main.cc)
    target_compile_definitions(mcmc_sampling_complex PRIVATE USE_COMPLEX)
    target_link_libraries(mcmc_sampling_complex PRIVATE tnbp_test_backend)
    add_subdirectory(heavyhex51)
    add_executable(mcmc_regression main.cc)
    target_link_libraries(mcmc_regression PRIVATE tnbp_test_backend)
    if(BUILD_TESTING)
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
endif()
