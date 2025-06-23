set(SOURCE_LIST 
    # TODO
)

target_sources(${TARGET_NAME} PRIVATE ${SOURCE_LIST})
target_link_libraries(${TARGET_NAME} PUBLIC gtest gtest_main)
target_link_libraries(${TARGET_NAME} PUBLIC lbmLib)
