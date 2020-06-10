add_executable(i-adhore threadPool.cpp higherLevel.cpp levelTwo.cpp alignComp.cpp SvgWriter.cpp AlignmentDrawer.cpp parallel.cpp SynthenicCloud.cpp BaseCluster.cpp Cluster.cpp DataSet.cpp GHM.cpp GHMProfile.cpp Gene.cpp GeneFamily.cpp GeneList.cpp GenePairs.cpp ListElement.cpp Multiplicon.cpp Profile.cpp Settings.cpp iADHoRe.cpp hpmath.cpp util.cpp)
target_link_libraries(i-adhore bmp alignment ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${PNG_LIBRARIES})
install(TARGETS i-adhore RUNTIME DESTINATION bin)

#add_executable(i-align threadPool.cpp higherLevel.cpp levelTwo.cpp alignComp.cpp SvgWriter.cpp AlignmentDrawer.cpp parallel.cpp AlignDataSet.cpp SynthenicCloud.cpp BaseCluster.cpp Cluster.cpp DataSet.cpp GHM.cpp GHMProfile.cpp Gene.cpp GeneFamily.cpp GeneList.cpp GenePairs.cpp ListElement.cpp Multiplicon.cpp Profile.cpp Settings.cpp iALIGN.cpp hpmath.cpp util.cpp)
#target_link_libraries(i-align bmp alignment ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${PNG_LIBRARIES})
#install(TARGETS i-align RUNTIME DESTINATION bin)

add_executable(i-visualize PostProcessor.cpp AlignmentVisualizer.cpp threadPool.cpp higherLevel.cpp levelTwo.cpp alignComp.cpp SvgWriter.cpp AlignmentDrawer.cpp parallel.cpp SynthenicCloud.cpp BaseCluster.cpp Cluster.cpp DataSet.cpp GHM.cpp GHMProfile.cpp Gene.cpp GeneFamily.cpp GeneList.cpp GenePairs.cpp ListElement.cpp Multiplicon.cpp Profile.cpp Settings.cpp hpmath.cpp util.cpp)
target_link_libraries(i-visualize bmp alignment ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${PNG_LIBRARIES})
install(TARGETS i-visualize RUNTIME DESTINATION bin)

add_subdirectory(alignment)
add_subdirectory(datastructures)
add_subdirectory(bmp)
