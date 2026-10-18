CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(TESTIVOLTABLE)

#---------------------------------------------------------

ENABLE_TESTING()

IF (NOT DEFINED IVOL_DIR)
  GET_FILENAME_COMPONENT(IVOL_ABSOLUTE_PATH "../../.." ABSOLUTE)
  SET(IVOL_DIR ${IVOL_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED IVOL_DIR)

IF (NOT DEFINED IVOL_DATA_DIR)
  CMAKE_PATH(APPEND IVOL_DATA_DIRA "${IVOL_DIR}" "data" "")
  CMAKE_PATH(APPEND IVOL_DATA_DIRB "$ENV{HOME}" "Data/Nrrd" "")
  CMAKE_PATH(APPEND IVOL_DATA_DIRC "/c/users/$ENV{USERNAME}/Nrrd" "" "")
  IF (EXISTS ${IVOL_DATA_DIRA})
    SET(IVOL_DATA_DIR ${IVOL_DATA_DIRA})
  ELSEIF (EXISTS ${IVOL_DATA_DIRB})
    SET(IVOL_DATA_DIR ${IVOL_DATA_DIRB})
    ELSEIF (EXISTS ${IVOL_DATA_DIRC})
    SET(IVOL_DATA_DIR ${IVOL_DATA_DIRC})
  ELSE()
    MESSAGE( FATAL_ERROR "Unable to find data directory. Directories ${IVOL_DATA_DIRA} and ${IVOL_DATA_DIRB} and ${IVOL_DATA_DIRC} do not exist.")
  ENDIF()
  
ENDIF (NOT DEFINED IVOL_DATA_DIR)


# Compare two programs, ${NEW_PROGRAM} and ${PREV_PROGRAM}.
IF (NOT DEFINED NEW_PROGRAM)
  SET(NEW_PROGRAM "../ivoldual")
ENDIF()

IF (NOT DEFINED PREV_PROGRAM)
  SET(PREV_PROGRAM "../ivoldual-ctest-prototype")
ENDIF()

MESSAGE(STATUS "This ctest compares new program \"${NEW_PROGRAM}\" against existing program \"${PREV_PROGRAM}\".\n")
MESSAGE(STATUS "  Use \"cmake -D NO_COMMENTS=1\" to ignore .off file comments in comparing program outputs.\n")


#---------------------------------------------------------
# MACROS
#---------------------------------------------------------

# Add comparison test. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_X tname cargs isovalue1 isovalue2 infile outfile_suffix)

  IF (DEFINED NO_COMMENTS)
    # Set and use NO_COMMENTS to avoid unused variable warning.
    SET(NO_COMMENTS "-no_comments")
    SET(arg_list "${cargs} ${NO_COMMENTS} ${isovalue1} ${isovalue2} ${infile}")
  ELSE()
    SET(arg_list "${cargs} ${isovalue1} ${isovalue2} ${infile}")
  ENDIF()

ADD_TEST( ${tname} ${CMAKE_COMMAND} -Dprog0=${NEW_PROGRAM} -Dprog1=${PREV_PROGRAM} -Dtest_arg=${arg_list} -Doutput0=out0.${outfile_suffix} -Doutput1=out1.${outfile_suffix} -Dflag_write_to_file=true -P ${CMAKE_SOURCE_DIR}/compare_prog.cmake )
  
ENDMACRO()


# Add standard comparison test. Output to .off files.
MACRO(ADD_COMPARISON_TEST tname cargs isovalue infile)
  ADD_COMPARISON_TEST_X(${tname} "${cargs}" ${isovalue} ${infile} "off")
ENDMACRO()


# Add comparison test on tdata.
MACRO(ADD_COMPARISON_TEST_ON_DATA_X tdata test_suffix cargs output_suffix)

  IF ( "${test_suffix}" STREQUAL "")
    SET( tname "${tdata}" )
  ELSE()
    SET( tname ${tdata}_${test_suffix} )
  ENDIF()

  ADD_COMPARISON_TEST_X( ${tname} "${cargs}"
                         "${${tdata}_isovalueI}" "${${tdata}_isovalueII}"
                         "${${tdata}_file}" "${output_suffix}" )
  SET_PROPERTY( TEST ${tname} PROPERTY LABELS ${${tdata}_speed} )

  IF ( "${test_suffix}" STREQUAL "")
    SET_PROPERTY( TEST ${tname} APPEND PROPERTY LABELS "basic" )
  ENDIF()

ENDMACRO()


# Set FLAG_SKIP_TEST to TRUE for certain combinations of tests and cargs. 
FUNCTION(QUERY_SKIP_TEST tdata cargs)

  # Default
  SET(FLAG_SKIP_TEST FALSE PARENT_SCOPE)

   IF ( "${cargs}" MATCHES ".*-supersample.*" )
     IF ( NOT ("${${tdata}_speed}" MATCHES "^veryfast" ) )
       SET(FLAG_SKIP_TEST TRUE PARENT_SCOPE)
    ENDIF()
  ENDIF()

ENDFUNCTION()


# Add comparison teston all data. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_ON_ALL_DATA_X test_suffix cargs output_suffix)

  FOREACH (tdata ${testdata_list})

    QUERY_SKIP_TEST("${tdata}" "${cargs}") 
    IF (${FLAG_SKIP_TEST})
      CONTINUE()
    ENDIF()

    ADD_COMPARISON_TEST_ON_DATA_X( "${tdata}" "${test_suffix}" 
    				   "${cargs}" "${output_suffix}" )

  ENDFOREACH(tdata)

ENDMACRO()


MACRO(ADD_COMPARISON_TEST_ON_ALL_DATA test_suffix cargs)
  ADD_COMPARISON_TEST_ON_ALL_DATA_X( "${test_suffix}" "${cargs}" "off" )
ENDMACRO()


# Add comparison test on all 3D data. Specify output file suffix.
MACRO(ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X test_suffix cargs output_suffix)

  FOREACH (tdata ${test3Ddata_list})
    QUERY_SKIP_TEST("${tdata}" "${cargs}") 
    IF (${FLAG_SKIP_TEST})
      CONTINUE()
    ENDIF()

    IF (NOT(${${tdata}_dimension} EQUAL 3))
      CONTINUE()
    ENDIF()

    ADD_COMPARISON_TEST_ON_DATA_X( "${tdata}" "${test_suffix}" 
    				   "${cargs}" "${output_suffix}" )

  ENDFOREACH(tdata)

ENDMACRO()


# Add comparison test on all 3D data.
MACRO(ADD_COMPARISON_TEST_ON_ALL_3D_DATA test_suffix cargs)
  ADD_COMPARISON_TEST_ON_ALL_3D_DATA_X("${test_suffix}" "${cargs}" "off")
ENDMACRO()


MACRO(SET_DATA_PROP tdata isovalueI isovalueII dimension tfile tspeed)
  SET(${tdata}_isovalueI ${isovalueI})
  SET(${tdata}_isovalueII ${isovalueII})
  SET(${tdata}_file ${tfile})
  SET(${tdata}_dimension ${dimension})
  SET(${tdata}_speed ${tspeed})
ENDMACRO()

MACRO(SET_3D_DATA_PROP tdata isovalueI isovalueII tfile tspeed)
  SET_DATA_PROP("${tdata}" ${isovalueI} ${isovalueII} 3 "${tfile}" ${tspeed})
ENDMACRO()


#---------------------------------------------------------
# SET VARIABLES
#---------------------------------------------------------

SET(volvis_dir "${IVOL_DATA_DIR}/3D/volvis")
SET(vollib_dir "${IVOL_DATA_DIR}/3D/vollib")

# 3D test sets, integer isovalues
SET_3D_DATA_PROP( Silicium_80_160 80 160 ${volvis_dir}/silicium.nhdr "veryfastII" )
SET_3D_DATA_PROP( Fuel_70_150 70 170 ${volvis_dir}/fuel.nhdr "veryfastII" )
SET_3D_DATA_PROP( Neghip_70_200 70 200 ${volvis_dir}/neghip.nhdr "veryfastII" )
SET_3D_DATA_PROP( MarschnerLobb_80_120 80 120 ${volvis_dir}/marschnerlobb.nhdr "veryfastII" )
SET_3D_DATA_PROP( Engine_70_140 70 140 "${volvis_dir}/engine.nhdr" "fast" )
SET_3D_DATA_PROP( Aneurysm_30_255 30 255 "${volvis_dir}/aneurysm.nhdr" "fastII" )
SET_3D_DATA_PROP( Frog_45_100 45 100 "${vollib_dir}/Frog.nhdr" "fastII" )
SET_3D_DATA_PROP( Lobster_20_150 20 150 "${volvis_dir}/lobster.nhdr" "fast" )

SET(test3DdataI_list Silicium_80_160 Fuel_70_150 Neghip_70_200 MarschnerLobb_80_120
    Engine_70_140 Aneurysm_30_255 Frog_45_100 Lobster_20_150 )
  
# 3D test sets, non-integer isovalues
SET_3D_DATA_PROP( Silicium_80p5_160p5 80.5 160.5 ${volvis_dir}/silicium.nhdr "veryfastII" )
SET_3D_DATA_PROP( Fuel_70p5_150p5 70.5 150.5 ${volvis_dir}/fuel.nhdr "veryfastII" )
SET_3D_DATA_PROP( Neghip_70p5_200p5 80.5 200.5 ${volvis_dir}/neghip.nhdr "veryfastII" )
SET_3D_DATA_PROP( MarschnerLobb_80p5_120p5 80.5 120.5
                  "${volvis_dir}/marschnerlobb.nhdr" "veryfastII" )
SET_3D_DATA_PROP( Engine_70p5_140p5 70.5 140.5
                  "${volvis_dir}/engine.nhdr" "fastII" )
SET_3D_DATA_PROP( Aneurysm_30p5_255p5 30.5 255.5
                  "${volvis_dir}/aneurysm.nhdr" "fastII" )
SET_3D_DATA_PROP( Frog_45p5_100p5 45.5 100.5
                  "${vollib_dir}/Frog.nhdr" "fastII" )
SET_3D_DATA_PROP( Lobster_20p5_150p5 20.5 150.5
                  "${volvis_dir}/lobster.nhdr" "fast" )
SET(test3DdataNI_list Silicium_20p5_100p5 Fuel_70p5_150p5 Neghip_70p5_200p5
                      MarschnerLobb_80p5_120p5 Engine_70p5_140p5
                      Aneurysm_30p5_255p5 Frog_45p5_100p5
                      Lobster_20p5_150p5 )  

# All test sets
SET(test3Ddata_list ${test3DdataI_list} ${test3DdataNI_list})
SET(testdata_list ${test3Ddata_list})


#---------------------------------------------------------
# ADD TESTS
#---------------------------------------------------------

ADD_COMPARISON_TEST_ON_ALL_DATA("" "")
ADD_COMPARISON_TEST_ON_ALL_DATA("subsample" "-subsample 2")
ADD_COMPARISON_TEST_ON_ALL_DATA("subdivide" "-subdivide")
ADD_COMPARISON_TEST_ON_ALL_DATA("split_hex" "-subdivide -split_hex")
ADD_COMPARISON_TEST_ON_ALL_DATA("add_outer_layer" "-subdivide -add_outer_layer")


