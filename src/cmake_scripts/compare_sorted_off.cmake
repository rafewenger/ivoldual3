# Compare sorted off file created by two programs.
# Sort using ijkmesh sort.

#---------------------------------------------------------
# MACROS
#---------------------------------------------------------

MACRO(RUN_COMMAND commandX commandX_args commandX_outfile)

IF (flag_write_to_file)
  SET (arg_list "-o ${commandX_outfile} ${commandX_args}")
  SEPARATE_ARGUMENTS(arg_list)

  EXECUTE_PROCESS(
  COMMAND ${commandX} ${arg_list}
  RESULT_VARIABLE commandX_failed
)
ELSE (flag_write_to_file)
  SET (arg_list "${commandX_args}")
  SEPARATE_ARGUMENTS(arg_list)

  EXECUTE_PROCESS(
    COMMAND ${commandX} ${arg_list}
    RESULT_VARIABLE commandX_failed
    OUTPUT_FILE ${commandX_outfile}
  )
ENDIF (flag_write_to_file)

IF ( commandX_failed )
  message( SEND_ERROR "Command ${commandX} failed on arguments ${arg_list}." )
ENDIF ( commandX_failed )

# Sort output off file.
SET (ijkmesh_arg_list "sort ${commandX_outfile}")
SEPARATE_ARGUMENTS(ijkmesh_arg_list)

EXECUTE_PROCESS(
  COMMAND ijkmesh ${ijkmesh_arg_list}
  RESULT_VARIABLE sort_failed
)

IF (sort_failed)
  message( SEND_ERROR "Error sorting off file ${commandX_outfile}." )
ENDIF ( sort_failed )

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

# test_arg is a string containing the command arguments
IF ( NOT test_arg )
  SET ( test_arg "" )
ENDIF ( NOT test_arg )

# output0 contains the name of the prog0 output file
IF ( NOT output0 )
  message( FATAL_ERROR "Variable output0 not defined." )
ENDIF ( NOT output0 )

# output1 contains the name of the prog1 output file
IF ( NOT output1 )
  message( FATAL_ERROR "Variable output1 not defined." )
ENDIF ( NOT output1 )

# Eliminate extra spaces
STRING(REGEX REPLACE "[ ][ ]+" " " test_arg "${test_arg}")
STRING(REGEX REPLACE "^[ ]+" "" test_arg "${test_arg}")

RUN_COMMAND(${prog0} "${test_arg}" "${output0}")
RUN_COMMAND(${prog1} "${test_arg}" "${output1}")

EXECUTE_PROCESS( 
  COMMAND ${CMAKE_COMMAND} -E compare_files ${output0} ${output1}
  RESULT_VARIABLE output_mismatch
)

IF ( output_mismatch )
  message( SEND_ERROR "Output of ${prog1} does not match output of ${prog0}." )
ENDIF ( output_mismatch )
