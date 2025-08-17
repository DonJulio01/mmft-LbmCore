set(LBM_SOURCES
    lbm.cpp
)

target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(${TARGET_NAME} PRIVATE ${LBM_SOURCES})
