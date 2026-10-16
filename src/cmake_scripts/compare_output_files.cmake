# Compare the output of two programs.

#---------------------------------------------------------
# MACROS
#---------------------------------------------------------

MACRO(RUN_COMMAND commandX commandX_args)

SET (arg_list "${commandX_args}")
SEPARATE_ARGUMENTS(arg_list)

EXECUTE_PROCESS(
  COMMAND ${commandX} ${arg_list}
  RESULT_VARIABLE commandX_failed
)

IF ( no_fail_on_return1 )
  IF ( commandX_failed AND NOT ${commandX_failed} EQUAL 1)
    message( SEND_ERROR "Command ${commandX} failed on arguments ${arg_list}." )
  ENDIF ( commandX_failed AND NOT ${commandX_failed} EQUAL 1)
ELSE ( no_fail_on_return1 )
  IF ( commandX_failed )
    message( SEND_ERROR "Command ${commandX} failed on arguments ${arg_list}." )
  ENDIF ( commandX_failed )
ENDIF( no_fail_on_return1 )

ENDMACRO()


MACRO(GET_OFF_FILE_PREFIX off_filename file_prefix)
  STRING(REGEX REPLACE "(.)*/" "" file_prefix2 ${off_filename})                
  STRING(REGEX REPLACE ".off$" "" file_prefix ${file_prefix2})
ENDMACRO()

# Store output files in strings.
MACRO(STORE_OUTPUT_FILES file_prefix suffix_list prefix)

  FOREACH(suffix ${suffix_list})

    SET(out_file ${file_prefix}.${suffix})
    IF (NOT EXISTS ${out_file})
      MESSAGE(SEND_ERROR "Output file ${out_file} does not exist.")
    ENDIF()

    FILE(STRINGS ${out_file} ${prefix}_${suffix})

  ENDFOREACH()

ENDMACRO()


# Remove output files in strings.
MACRO(REMOVE_OUTPUT_FILES file_prefix suffix_list)

  FOREACH(suffix ${suffix_list})

    SET(out_file ${file_prefix}.${suffix})
    FILE(REMOVE ${out_file})
  ENDFOREACH()

ENDMACRO()


# Compare strings representing output files.
MACRO(COMPARE_OUTPUT_FILE_STRINGS prefix0 prefix1 suffix_list)

  FOREACH(suffix ${suffix_list})
    IF (NOT (${prefix0}_${suffix} STREQUAL ${prefix1}_${suffix}))
      message(SEND_ERROR "Output files .${suffix} of ${prog0} and ${prog1} do not match.")
    ENDIF()

  ENDFOREACH()

ENDMACRO()


#---------------------------------------------------------

# Check arguments
# prog0 is the first program
IF ( NOT prog0 )
  message( FATAL_ERROR "Variable prog0 not defined." )
ENDIF ( NOT prog0 )

# prog1 is the second program
IF ( NOT prog1 )
  message( FATAL_ERROR "Variable prog1 not defined." )
ENDIF ( NOT prog1 )

# infile is the input file
IF ( NOT infile )
  message( FATAL_ERROR "Variable infile not defined." )
ENDIF ( NOT infile )

# suffix_list is a list of output file suffixes
IF ( NOT suffix_list )
  message( FATAL_ERROR "Variable suffix_list not defined." )
ENDIF ( NOT suffix_list )

# test_arg is a string containing the command arguments
IF ( NOT test_arg )
  SET ( test_arg "" )
ENDIF ( NOT test_arg )

# Change " " in suffix_list back to ";"
SEPARATE_ARGUMENTS(suffix_list)

# Eliminate extra spaces
STRING(REGEX REPLACE "[ ][ ]+" " " test_arg "${test_arg}")
STRING(REGEX REPLACE "^[ ]+" "" test_arg "${test_arg}")

IF (NOT (DEFINED file_prefix))
  GET_OFF_FILE_PREFIX(${infile} file_prefix)
ENDIF()

# Remove any files with output file filenames.
REMOVE_OUTPUT_FILES(${file_prefix} "${suffix_list}")

RUN_COMMAND(${prog0} "${test_arg} ${infile}")
STORE_OUTPUT_FILES(${file_prefix} "${suffix_list}" "out0")
REMOVE_OUTPUT_FILES(${file_prefix} "${suffix_list}")

RUN_COMMAND(${prog1} "${test_arg} ${infile}")
STORE_OUTPUT_FILES(${file_prefix} "${suffix_list}" "out1")
REMOVE_OUTPUT_FILES(${file_prefix} "${suffix_list}")

COMPARE_OUTPUT_FILE_STRINGS("out0" "out1" "${suffix_list}")

