m4_include([m4/ac_cxx_namespaces.m4])
m4_include([m4/ac_cxx_have_stl.m4])

# SYNOPSIS
#
# AX_ARG_WITH(lib_name, headers ..., paths_to_check ..., help_string)
#
# DESCRIPTION
#
# This macro defines the helper argument "--with-$libname" using AC_ARG_WITH and
# probes the given (as 2nd argument) headers first in system-wide locations and
# then for the specified space-separated locations (given as 3rd argument). If
# probing succeeds it adds the headers location to $CPPFLAGS and library locations
# to $LDFLAGS. The special path "auto" means that the library will be autoprobed
# in user's $HOME. If headers have been located this macro defines the variable
# $ac_lib_{lib_name} (which is set to "yes") and also in case of "auto" location
# appends the found dirs with headers to $CPPFLAGS and appends the directories
# with library binaries to $LDFLAGS.
#
# TODO: If directories contain spaces this will cause problems (fixing $IFS will
# cause problems in other places).
#
AC_DEFUN([AX_PROBE_LIBRARY], [
	AC_ARG_WITH([$1], [AC_HELP_STRING([--with-$1], [$4 (default: "$3")])], [], [with_$1="$3"])

	dnl Testing for default locations, ignoring the optional locations:   
	AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

	dnl Testing the specified locations:   
	if test "${ac_lib_$1}" != "yes" -a "${with_$1}" != ""
	then
		AX_RESET_HEADERS_CACHE([$2])

		for ac_test_location in ${with_$1}
		do
			dnl Probing the library in user's $HOME:
			if test "${ac_test_location}" = "auto"
			then
				dnl Read the directory entries by mask sorted alphabetically in reverse order:
				for ac_location in `ls -1d $HOME/$1-* 2>/dev/null | tac` 
				do
					if test -d "${ac_location}"
					then
						dnl Save the current state
						ax_probe_library_save_LDFLAGS=${LDFLAGS}
						ax_probe_library_save_CPPFLAGS=${CPPFLAGS}
						
						dnl Compose the list of unique locations of headers:
						for ac_inc_location in `find "${ac_location}" -iname '*.h' |
							while read ac_include_location; do dirname "${ac_include_location}"; done |
								sort -u`
						do
							CPPFLAGS="-I${ac_inc_location} ${CPPFLAGS}"
						done
									

						dnl Compose the list of unique locations of libraries (standard library extensions are taken from autoconf/libs.m4:185):
						for ac_lib_location in `find "${ac_location}" -iname '*.so' -o -iname '*.sl' -o -iname '*.dylib' -o -iname '*.a' -o -iname '*.dll' |
							while read ac_library_location; do dirname "${ac_library_location}"; done |
								sort -u`
						do
							LDFLAGS="-L${ac_lib_location} ${LDFLAGS}"
						done
						
						AC_MSG_CHECKING([$1 for $2 in ${ac_location}])
						AS_ECHO()
						_AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}" and LDFLAGS="${LDFLAGS}"])
						
						AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])
						
						dnl We have found the location, leave the loop:
						if test "${ac_lib_$1}" = "yes"
						then
							break 2;
						fi
						
						dnl Restore the state to original in case of unsuccessful attempt
						LDFLAGS=${ax_probe_library_save_LDFLAGS}
						CPPFLAGS=${ax_probe_library_save_CPPFLAGS}
						AX_RESET_HEADERS_CACHE([$2])
					fi
				done
			else
				dnl Save the current state
				ax_probe_library_save_CPPFLAGS=${CPPFLAGS}

				CPPFLAGS="-I${ac_test_location} $CPPFLAGS"
				
				AC_MSG_CHECKING([$1 for $2 in ${ac_test_location}])
				AS_ECHO()
				_AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}"])

				AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

				dnl We have found the location, leave the loop:
				if test "${ac_lib_$1}" = "yes"
				then
					break;
				fi

				dnl Restore the state to original in case of unsuccessful attempt
				CPPFLAGS=${ax_probe_library_save_CPPFLAGS}
				AX_RESET_HEADERS_CACHE([$2])
			fi
		done
	fi

	if test "${ac_lib_$1}" != "yes"
	then
		AC_MSG_ERROR([$2 header(s) is missing. Check the default/listed above headers locations.])
	fi
]) # AX_PROBE_LIBRARY

# SYNOPSIS
#
# AX_RESET_HEADERS_CACHE(headers ...)
#
# DESCRIPTION
#
# This macro invalidates the headers cache variables created by previous AC_CHECK_HEADER/AC_CHECK_HEADERS checks.
#
AC_DEFUN([AX_RESET_HEADERS_CACHE], [
	AS_FOR([AX_var], [ax_var], [$1], [
		dnl You can replace "ac_cv_header_" with any prefix from http://www.gnu.org/software/autoconf/manual/html_node/Cache-Variable-Index.html
		AS_VAR_PUSHDEF([ax_Var], [ac_cv_header_${ax_var}])
		AS_UNSET([ax_Var])
		AS_VAR_POPDEF([ax_Var])
	])
]) # AX_RESET_HEADERS_CACHE

# SYNOPSIS
#
# AX_TRY_LINK(library, includes, function-body [, action-if-true [, action-if-false]])
#
# DESCRIPTION
#
# This macro is a combination of autoconf's AC_TRY_LINK/AC_CHECK_LIB that checks the given given C++ program (3rd argument) successfully compiles
# and adds the library (1st argument) to the $LIBS list (keeping this list unique).
#
AC_DEFUN([AX_TRY_LINK], [
	dnl Below logic is a workaround for the limitation, that variables may not allow
	dnl symbols like "+" or "-". See AC_CHECK_LIB source comments for more information.
	m4_ifval([$4], , [AH_CHECK_LIB([$1])])
	AS_LITERAL_IF([$1],
		[AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_$1_$2])],
		[AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_$1''_$2])])

	AC_CACHE_CHECK([for -l$1], [ac_Lib], [
		dnl Save the current state
		AC_LANG_SAVE
		AC_LANG_CPLUSPLUS
		ax_try_link_save_LIBS=${LIBS}
		LIBS="-l$1 ${LIBS}"

		AC_TRY_LINK([$2], [$3], [AS_VAR_SET([ac_Lib], [yes])], [AS_VAR_SET([ac_Lib], [no])])

		dnl Restore the state to original regardless to the result
		LIBS=${ax_try_link_save_LIBS}
		AC_LANG_RESTORE
	])

	dnl If the variable is set, we define a constant and push library to $LIBS by default or execute $4, otherwise execute $5.
	AS_VAR_IF([ac_Lib], [yes],
		[m4_default([$4], [
			AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_LIB$1))
			dnl Do not prepend a library, if it is already in the list:
			(echo "${LIBS}" | grep -q -- "-l$1 ") || LIBS="-l$1 ${LIBS}"
		])],
		[$5]
	)
	AS_VAR_POPDEF([ac_Lib])
]) # AX_TRY_LINK
