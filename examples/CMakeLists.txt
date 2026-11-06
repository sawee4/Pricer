foreach(example basic_usage price_european)
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} PRIVATE pricer::pricer pricer_warnings)
endforeach()
