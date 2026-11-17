foreach(example basic_usage price_european convergence_study
        variance_reduction_study binomial_convergence
        mc_greeks_study)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE pricer::pricer pricer_warnings)
endforeach()
