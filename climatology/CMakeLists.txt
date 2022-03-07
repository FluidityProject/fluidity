add_executable(create_climatology_atlas
  EXCLUDE_FROM_ALL
  create_climatology_atlas.cpp
)
target_link_libraries(create_climatology_atlas
  fluidity
)
