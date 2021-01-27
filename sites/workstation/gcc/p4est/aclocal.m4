# generated automatically by aclocal 1.16.2 -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
m4_if(m4_defn([AC_AUTOCONF_VERSION]), [2.69],,
[m4_warning([this file was generated for autoconf 2.69.
You have another version of autoconf.  It may work, but is not guaranteed to.
If you have problems, you may need to regenerate the build system entirely.
To do so, use the procedure documented by the package, typically 'autoreconf'.])])

# Copyright (C) 2002-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
# (This private macro should not be called outside this file.)
AC_DEFUN([AM_AUTOMAKE_VERSION],
[am__api_version='1.16'
dnl Some users find AM_AUTOMAKE_VERSION and mistake it for a way to
dnl require some minimum version.  Point them to the right macro.
m4_if([$1], [1.16.2], [],
      [AC_FATAL([Do not call $0, use AM_INIT_AUTOMAKE([$1]).])])dnl
])

# _AM_AUTOCONF_VERSION(VERSION)
# -----------------------------
# aclocal traces this macro to find the Autoconf version.
# This is a private macro too.  Using m4_define simplifies
# the logic in aclocal, which can simply ignore this definition.
m4_define([_AM_AUTOCONF_VERSION], [])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION and AM_AUTOMAKE_VERSION so they can be traced.
# This function is AC_REQUIREd by AM_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
[AM_AUTOMAKE_VERSION([1.16.2])dnl
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
_AM_AUTOCONF_VERSION(m4_defn([AC_AUTOCONF_VERSION]))])

# Copyright (C) 2011-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_AR([ACT-IF-FAIL])
# -------------------------
# Try to determine the archiver interface, and trigger the ar-lib wrapper
# if it is needed.  If the detection of archiver interface fails, run
# ACT-IF-FAIL (default is to abort configure with a proper error message).
AC_DEFUN([AM_PROG_AR],
[AC_BEFORE([$0], [LT_INIT])dnl
AC_BEFORE([$0], [AC_PROG_LIBTOOL])dnl
AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([ar-lib])dnl
AC_CHECK_TOOLS([AR], [ar lib "link -lib"], [false])
: ${AR=ar}

AC_CACHE_CHECK([the archiver ($AR) interface], [am_cv_ar_interface],
  [AC_LANG_PUSH([C])
   am_cv_ar_interface=ar
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([[int some_variable = 0;]])],
     [am_ar_try='$AR cru libconftest.a conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
      AC_TRY_EVAL([am_ar_try])
      if test "$ac_status" -eq 0; then
        am_cv_ar_interface=ar
      else
        am_ar_try='$AR -NOLOGO -OUT:conftest.lib conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
        AC_TRY_EVAL([am_ar_try])
        if test "$ac_status" -eq 0; then
          am_cv_ar_interface=lib
        else
          am_cv_ar_interface=unknown
        fi
      fi
      rm -f conftest.lib libconftest.a
     ])
   AC_LANG_POP([C])])

case $am_cv_ar_interface in
ar)
  ;;
lib)
  # Microsoft lib, so override with the ar-lib wrapper script.
  # FIXME: It is wrong to rewrite AR.
  # But if we don't then we get into trouble of one sort or another.
  # A longer-term fix would be to have automake use am__AR in this case,
  # and then we could set am__AR="$am_aux_dir/ar-lib \$(AR)" or something
  # similar.
  AR="$am_aux_dir/ar-lib $AR"
  ;;
unknown)
  m4_default([$1],
             [AC_MSG_ERROR([could not determine $AR interface])])
  ;;
esac
AC_SUBST([AR])dnl
])

# AM_AUX_DIR_EXPAND                                         -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to '$srcdir/foo'.  In other projects, it is set to
# '$srcdir', '$srcdir/..', or '$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is '.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

AC_DEFUN([AM_AUX_DIR_EXPAND],
[AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
# Expand $ac_aux_dir to an absolute path.
am_aux_dir=`cd "$ac_aux_dir" && pwd`
])

# AM_CONDITIONAL                                            -*- Autoconf -*-

# Copyright (C) 1997-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[AC_PREREQ([2.52])dnl
 m4_if([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
       [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])dnl
AC_SUBST([$1_FALSE])dnl
_AM_SUBST_NOTMAKE([$1_TRUE])dnl
_AM_SUBST_NOTMAKE([$1_FALSE])dnl
m4_define([_AM_COND_VALUE_$1], [$2])dnl
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([[conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.]])
fi])])

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# There are a few dirty hacks below to avoid letting 'AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "OBJC", "OBJCXX", "UPC", or "GJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

m4_if([$1], [CC],   [depcc="$CC"   am_compiler_list=],
      [$1], [CXX],  [depcc="$CXX"  am_compiler_list=],
      [$1], [OBJC], [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
      [$1], [OBJCXX], [depcc="$OBJCXX" am_compiler_list='gcc3 gcc'],
      [$1], [UPC],  [depcc="$UPC"  am_compiler_list=],
      [$1], [GCJ],  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                    [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named 'D' -- because '-MD' means "put the output
  # in D".
  rm -rf conftest.dir
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir
  # We will build objects and dependencies in a subdirectory because
  # it helps to detect inapplicable dependency modes.  For instance
  # both Tru64's cc and ICC support -MD to output dependencies as a
  # side effect of compilation, but ICC will put the dependencies in
  # the current directory while Tru64 will put them in the object
  # directory.
  mkdir sub

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  am__universal=false
  m4_case([$1], [CC],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac],
    [CXX],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac])

  for depmode in $am_compiler_list; do
    # Setup a source with many dependencies, because some compilers
    # like to wrap large dependency lists on column 80 (with \), and
    # we should not choose a depcomp mode which is confused by this.
    #
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    : > sub/conftest.c
    for i in 1 2 3 4 5 6; do
      echo '#include "conftst'$i'.h"' >> sub/conftest.c
      # Using ": > sub/conftst$i.h" creates only sub/conftst1.h with
      # Solaris 10 /bin/sh.
      echo '/* dummy */' > sub/conftst$i.h
    done
    echo "${am__include} ${am__quote}sub/conftest.Po${am__quote}" > confmf

    # We check with '-c' and '-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle '-M -o', and we need to detect this.  Also, some Intel
    # versions had trouble with output in subdirs.
    am__obj=sub/conftest.${OBJEXT-o}
    am__minus_obj="-o $am__obj"
    case $depmode in
    gcc)
      # This depmode causes a compiler race in universal mode.
      test "$am__universal" = false || continue
      ;;
    nosideeffect)
      # After this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested.
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    msvc7 | msvc7msys | msvisualcpp | msvcmsys)
      # This compiler won't grok '-c -o', but also, the minuso test has
      # not run yet.  These depmodes are late enough in the game, and
      # so weak that their functioning should not be impacted.
      am__obj=conftest.${OBJEXT-o}
      am__minus_obj=
      ;;
    none) break ;;
    esac
    if depmode=$depmode \
       source=sub/conftest.c object=$am__obj \
       depfile=sub/conftest.Po tmpdepfile=sub/conftest.TPo \
       $SHELL ./depcomp $depcc -c $am__minus_obj sub/conftest.c \
         >/dev/null 2>conftest.err &&
       grep sub/conftst1.h sub/conftest.Po > /dev/null 2>&1 &&
       grep sub/conftst6.h sub/conftest.Po > /dev/null 2>&1 &&
       grep $am__obj sub/conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      # icc doesn't choke on unknown options, it will just issue warnings
      # or remarks (even with -Werror).  So we grep stderr for any message
      # that says an option was ignored or not supported.
      # When given -MP, icc 7.0 and 7.1 complain thusly:
      #   icc: Command line warning: ignoring option '-M'; no argument required
      # The diagnosis changed in icc 8.0:
      #   icc: Command line remark: option '-MP' not supported
      if (grep 'ignoring option' conftest.err ||
          grep 'not supported' conftest.err) >/dev/null 2>&1; then :; else
        am_cv_$1_dependencies_compiler_type=$depmode
        break
      fi
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
AM_CONDITIONAL([am__fastdep$1], [
  test "x$enable_dependency_tracking" != xno \
  && test "$am_cv_$1_dependencies_compiler_type" = gcc3])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES.
AC_DEFUN([AM_SET_DEPDIR],
[AC_REQUIRE([AM_SET_LEADING_DOT])dnl
AC_SUBST([DEPDIR], ["${am__leading_dot}deps"])dnl
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE([dependency-tracking], [dnl
AS_HELP_STRING(
  [--enable-dependency-tracking],
  [do not reject slow dependency extractors])
AS_HELP_STRING(
  [--disable-dependency-tracking],
  [speeds up one-time build])])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
  am__nodep='_no'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])dnl
_AM_SUBST_NOTMAKE([AMDEPBACKSLASH])dnl
AC_SUBST([am__nodep])dnl
_AM_SUBST_NOTMAKE([am__nodep])dnl
])

# Generate code to set up dependency tracking.              -*- Autoconf -*-

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[{
  # Older Autoconf quotes --file arguments for eval, but not when files
  # are listed without --file.  Let's play safe and only enable the eval
  # if we detect the quoting.
  # TODO: see whether this extra hack can be removed once we start
  # requiring Autoconf 2.70 or later.
  AS_CASE([$CONFIG_FILES],
          [*\'*], [eval set x "$CONFIG_FILES"],
          [*], [set x $CONFIG_FILES])
  shift
  # Used to flag and report bootstrapping failures.
  am_rc=0
  for am_mf
  do
    # Strip MF so we end up with the name of the file.
    am_mf=`AS_ECHO(["$am_mf"]) | sed -e 's/:.*$//'`
    # Check whether this is an Automake generated Makefile which includes
    # dependency-tracking related rules and includes.
    # Grep'ing the whole file directly is not great: AIX grep has a line
    # limit of 2048, but all sed's we know have understand at least 4000.
    sed -n 's,^am--depfiles:.*,X,p' "$am_mf" | grep X >/dev/null 2>&1 \
      || continue
    am_dirpart=`AS_DIRNAME(["$am_mf"])`
    am_filepart=`AS_BASENAME(["$am_mf"])`
    AM_RUN_LOG([cd "$am_dirpart" \
      && sed -e '/# am--include-marker/d' "$am_filepart" \
        | $MAKE -f - am--depfiles]) || am_rc=$?
  done
  if test $am_rc -ne 0; then
    AC_MSG_FAILURE([Something went wrong bootstrapping makefile fragments
    for automatic dependency tracking.  If GNU make was not used, consider
    re-running the configure script with MAKE="gmake" (or whatever is
    necessary).  You can also try re-running configure with the
    '--disable-dependency-tracking' option to at least be able to build
    the package (albeit without support for automatic dependency tracking).])
  fi
  AS_UNSET([am_dirpart])
  AS_UNSET([am_filepart])
  AS_UNSET([am_mf])
  AS_UNSET([am_rc])
  rm -f conftest-deps.mk
}
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking is enabled.
# This creates each '.Po' and '.Plo' makefile fragment that we'll need in
# order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" MAKE="${MAKE-make}"])])

# Do all the work for Automake.                             -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This macro actually does too much.  Some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

dnl Redefine AC_PROG_CC to automatically invoke _AM_PROG_CC_C_O.
m4_define([AC_PROG_CC],
m4_defn([AC_PROG_CC])
[_AM_PROG_CC_C_O
])

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_PREREQ([2.65])dnl
dnl Autoconf wants to disallow AM_ names.  We explicitly allow
dnl the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl
AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
AC_REQUIRE([AC_PROG_INSTALL])dnl
if test "`cd $srcdir && pwd`" != "`pwd`"; then
  # Use -I$(srcdir) only when $(srcdir) != ., so that make's output
  # is not polluted with repeated "-I."
  AC_SUBST([am__isrc], [' -I$(srcdir)'])_AM_SUBST_NOTMAKE([am__isrc])dnl
  # test to see if srcdir already configured
  if test -f $srcdir/config.status; then
    AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
  fi
fi

# test whether we have cygpath
if test -z "$CYGPATH_W"; then
  if (cygpath --version) >/dev/null 2>/dev/null; then
    CYGPATH_W='cygpath -w'
  else
    CYGPATH_W=echo
  fi
fi
AC_SUBST([CYGPATH_W])

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[AC_DIAGNOSE([obsolete],
             [$0: two- and three-arguments forms are deprecated.])
m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
dnl Diagnose old-style AC_INIT with new-style AM_AUTOMAKE_INIT.
m4_if(
  m4_ifdef([AC_PACKAGE_NAME], [ok]):m4_ifdef([AC_PACKAGE_VERSION], [ok]),
  [ok:ok],,
  [m4_fatal([AC_INIT should be called with package and version arguments])])dnl
 AC_SUBST([PACKAGE], ['AC_PACKAGE_TARNAME'])dnl
 AC_SUBST([VERSION], ['AC_PACKAGE_VERSION'])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED([PACKAGE], ["$PACKAGE"], [Name of package])
 AC_DEFINE_UNQUOTED([VERSION], ["$VERSION"], [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG([ACLOCAL], [aclocal-${am__api_version}])
AM_MISSING_PROG([AUTOCONF], [autoconf])
AM_MISSING_PROG([AUTOMAKE], [automake-${am__api_version}])
AM_MISSING_PROG([AUTOHEADER], [autoheader])
AM_MISSING_PROG([MAKEINFO], [makeinfo])
AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
AC_REQUIRE([AM_PROG_INSTALL_STRIP])dnl
AC_REQUIRE([AC_PROG_MKDIR_P])dnl
# For better backward compatibility.  To be removed once Automake 1.9.x
# dies out for good.  For more background, see:
# <https://lists.gnu.org/archive/html/automake/2012-07/msg00001.html>
# <https://lists.gnu.org/archive/html/automake/2012-07/msg00014.html>
AC_SUBST([mkdir_p], ['$(MKDIR_P)'])
# We need awk for the "check" target (and possibly the TAP driver).  The
# system "awk" is bad on some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl
AC_REQUIRE([AM_SET_LEADING_DOT])dnl
_AM_IF_OPTION([tar-ustar], [_AM_PROG_TAR([ustar])],
	      [_AM_IF_OPTION([tar-pax], [_AM_PROG_TAR([pax])],
			     [_AM_PROG_TAR([v7])])])
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_CC],
		  [_AM_DEPENDENCIES([CC])],
		  [m4_define([AC_PROG_CC],
			     m4_defn([AC_PROG_CC])[_AM_DEPENDENCIES([CC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
		  [_AM_DEPENDENCIES([CXX])],
		  [m4_define([AC_PROG_CXX],
			     m4_defn([AC_PROG_CXX])[_AM_DEPENDENCIES([CXX])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJC],
		  [_AM_DEPENDENCIES([OBJC])],
		  [m4_define([AC_PROG_OBJC],
			     m4_defn([AC_PROG_OBJC])[_AM_DEPENDENCIES([OBJC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJCXX],
		  [_AM_DEPENDENCIES([OBJCXX])],
		  [m4_define([AC_PROG_OBJCXX],
			     m4_defn([AC_PROG_OBJCXX])[_AM_DEPENDENCIES([OBJCXX])])])dnl
])
AC_REQUIRE([AM_SILENT_RULES])dnl
dnl The testsuite driver may need to know about EXEEXT, so add the
dnl 'am__EXEEXT' conditional if _AM_COMPILER_EXEEXT was seen.  This
dnl macro is hooked onto _AC_COMPILER_EXEEXT early, see below.
AC_CONFIG_COMMANDS_PRE(dnl
[m4_provide_if([_AM_COMPILER_EXEEXT],
  [AM_CONDITIONAL([am__EXEEXT], [test -n "$EXEEXT"])])])dnl

# POSIX will say in a future version that running "rm -f" with no argument
# is OK; and we want to be able to make that assumption in our Makefile
# recipes.  So use an aggressive probe to check that the usage we want is
# actually supported "in the wild" to an acceptable degree.
# See automake bug#10828.
# To make any issue more visible, cause the running configure to be aborted
# by default if the 'rm' program in use doesn't match our expectations; the
# user can still override this though.
if rm -f && rm -fr && rm -rf; then : OK; else
  cat >&2 <<'END'
Oops!

Your 'rm' program seems unable to run without file operands specified
on the command line, even when the '-f' option is present.  This is contrary
to the behaviour of most rm programs out there, and not conforming with
the upcoming POSIX standard: <http://austingroupbugs.net/view.php?id=542>

Please tell bug-automake@gnu.org about your system, including the value
of your $PATH and any error possibly output before this message.  This
can help us improve future automake versions.

END
  if test x"$ACCEPT_INFERIOR_RM_PROGRAM" = x"yes"; then
    echo 'Configuration will proceed anyway, since you have set the' >&2
    echo 'ACCEPT_INFERIOR_RM_PROGRAM variable to "yes"' >&2
    echo >&2
  else
    cat >&2 <<'END'
Aborting the configuration process, to ensure you take notice of the issue.

You can download and install GNU coreutils to get an 'rm' implementation
that behaves properly: <https://www.gnu.org/software/coreutils/>.

If you want to complete the configuration process using your problematic
'rm' anyway, export the environment variable ACCEPT_INFERIOR_RM_PROGRAM
to "yes", and re-run configure.

END
    AC_MSG_ERROR([Your 'rm' program is bad, sorry.])
  fi
fi
dnl The trailing newline in this macro's definition is deliberate, for
dnl backward compatibility and to allow trailing 'dnl'-style comments
dnl after the AM_INIT_AUTOMAKE invocation. See automake bug#16841.
])

dnl Hook into '_AC_COMPILER_EXEEXT' early to learn its expansion.  Do not
dnl add the conditional right here, as _AC_COMPILER_EXEEXT may be further
dnl mangled by Autoconf and run in a shell conditional statement.
m4_define([_AC_COMPILER_EXEEXT],
m4_defn([_AC_COMPILER_EXEEXT])[m4_provide([_AM_COMPILER_EXEEXT])])

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  The stamp files are numbered to have different names.

# Autoconf calls _AC_AM_CONFIG_HEADER_HOOK (when defined) in the
# loop where config.status creates the headers, so we can generate
# our stamp files there.
AC_DEFUN([_AC_AM_CONFIG_HEADER_HOOK],
[# Compute $1's index in $config_headers.
_am_arg=$1
_am_stamp_count=1
for _am_header in $config_headers :; do
  case $_am_header in
    $_am_arg | $_am_arg:* )
      break ;;
    * )
      _am_stamp_count=`expr $_am_stamp_count + 1` ;;
  esac
done
echo "timestamp for $_am_arg" >`AS_DIRNAME(["$_am_arg"])`/stamp-h[]$_am_stamp_count])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.
AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
if test x"${install_sh+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    install_sh="\${SHELL} '$am_aux_dir/install-sh'" ;;
  *)
    install_sh="\${SHELL} $am_aux_dir/install-sh"
  esac
fi
AC_SUBST([install_sh])])

# Copyright (C) 2003-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# Check whether the underlying file-system supports filenames
# with a leading dot.  For instance MS-DOS doesn't.
AC_DEFUN([AM_SET_LEADING_DOT],
[rm -rf .tst 2>/dev/null
mkdir .tst 2>/dev/null
if test -d .tst; then
  am__leading_dot=.
else
  am__leading_dot=_
fi
rmdir .tst 2>/dev/null
AC_SUBST([am__leading_dot])])

# Check to see how 'make' treats includes.	            -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MAKE_INCLUDE()
# -----------------
# Check whether make has an 'include' directive that can support all
# the idioms we need for our automatic dependency tracking code.
AC_DEFUN([AM_MAKE_INCLUDE],
[AC_MSG_CHECKING([whether ${MAKE-make} supports the include directive])
cat > confinc.mk << 'END'
am__doit:
	@echo this is the am__doit target >confinc.out
.PHONY: am__doit
END
am__include="#"
am__quote=
# BSD make does it like this.
echo '.include "confinc.mk" # ignored' > confmf.BSD
# Other make implementations (GNU, Solaris 10, AIX) do it like this.
echo 'include confinc.mk # ignored' > confmf.GNU
_am_result=no
for s in GNU BSD; do
  AM_RUN_LOG([${MAKE-make} -f confmf.$s && cat confinc.out])
  AS_CASE([$?:`cat confinc.out 2>/dev/null`],
      ['0:this is the am__doit target'],
      [AS_CASE([$s],
          [BSD], [am__include='.include' am__quote='"'],
          [am__include='include' am__quote=''])])
  if test "$am__include" != "#"; then
    _am_result="yes ($s style)"
    break
  fi
done
rm -f confinc.* confmf.*
AC_MSG_RESULT([${_am_result}])
AC_SUBST([am__include])])
AC_SUBST([am__quote])])

# Fake the existence of programs that GNU maintainers use.  -*- Autoconf -*-

# Copyright (C) 1997-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])

# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it is modern enough.
# If it is, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([missing])dnl
if test x"${MISSING+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    MISSING="\${SHELL} \"$am_aux_dir/missing\"" ;;
  *)
    MISSING="\${SHELL} $am_aux_dir/missing" ;;
  esac
fi
# Use eval to expand $SHELL
if eval "$MISSING --is-lightweight"; then
  am_missing_run="$MISSING "
else
  am_missing_run=
  AC_MSG_WARN(['missing' script is too old or missing])
fi
])

# Helper functions for option handling.                     -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# --------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), [1])])

# _AM_SET_OPTIONS(OPTIONS)
# ------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[m4_foreach_w([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_CC_C_O
# ---------------
# Like AC_PROG_CC_C_O, but changed for automake.  We rewrite AC_PROG_CC
# to automatically call this.
AC_DEFUN([_AM_PROG_CC_C_O],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([compile])dnl
AC_LANG_PUSH([C])dnl
AC_CACHE_CHECK(
  [whether $CC understands -c and -o together],
  [am_cv_prog_cc_c_o],
  [AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
  # Make sure it works both with $CC and with simple cc.
  # Following AC_PROG_CC_C_O, we do the test twice because some
  # compilers refuse to overwrite an existing .o file with -o,
  # though they will create one.
  am_cv_prog_cc_c_o=yes
  for am_i in 1 2; do
    if AM_RUN_LOG([$CC -c conftest.$ac_ext -o conftest2.$ac_objext]) \
         && test -f conftest2.$ac_objext; then
      : OK
    else
      am_cv_prog_cc_c_o=no
      break
    fi
  done
  rm -f core conftest*
  unset am_i])
if test "$am_cv_prog_cc_c_o" != yes; then
   # Losing compiler, so override with the script.
   # FIXME: It is wrong to rewrite CC.
   # But if we don't then we get into trouble of one sort or another.
   # A longer-term fix would be to have automake use am__CC in this case,
   # and then we could set am__CC="\$(top_srcdir)/compile \$(CC)"
   CC="$am_aux_dir/compile $CC"
fi
AC_LANG_POP([C])])

# For backward compatibility.
AC_DEFUN_ONCE([AM_PROG_CC_C_O], [AC_REQUIRE([AC_PROG_CC])])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_RUN_LOG(COMMAND)
# -------------------
# Run COMMAND, save the exit status in ac_status, and log it.
# (This has been adapted from Autoconf's _AC_RUN_LOG macro.)
AC_DEFUN([AM_RUN_LOG],
[{ echo "$as_me:$LINENO: $1" >&AS_MESSAGE_LOG_FD
   ($1) >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   (exit $ac_status); }])

# Check to make sure that the build environment is sane.    -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Reject unsafe characters in $srcdir or the absolute working directory
# name.  Accept space and tab only in the latter.
am_lf='
'
case `pwd` in
  *[[\\\"\#\$\&\'\`$am_lf]]*)
    AC_MSG_ERROR([unsafe absolute working directory name]);;
esac
case $srcdir in
  *[[\\\"\#\$\&\'\`$am_lf\ \	]]*)
    AC_MSG_ERROR([unsafe srcdir value: '$srcdir']);;
esac

# Do 'set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   am_has_slept=no
   for am_try in 1 2; do
     echo "timestamp, slept: $am_has_slept" > conftest.file
     set X `ls -Lt "$srcdir/configure" conftest.file 2> /dev/null`
     if test "$[*]" = "X"; then
	# -L didn't work.
	set X `ls -t "$srcdir/configure" conftest.file`
     fi
     if test "$[*]" != "X $srcdir/configure conftest.file" \
	&& test "$[*]" != "X conftest.file $srcdir/configure"; then

	# If neither matched, then we have a broken ls.  This can happen
	# if, for instance, CONFIG_SHELL is bash and it inherits a
	# broken ls alias from the environment.  This has actually
	# happened.  Such a system could not be considered "sane".
	AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
  alias in your environment])
     fi
     if test "$[2]" = conftest.file || test $am_try -eq 2; then
       break
     fi
     # Just in case.
     sleep 1
     am_has_slept=yes
   done
   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT([yes])
# If we didn't sleep, we still need to ensure time stamps of config.status and
# generated files are strictly newer.
am_sleep_pid=
if grep 'slept: no' conftest.file >/dev/null 2>&1; then
  ( sleep 1 ) &
  am_sleep_pid=$!
fi
AC_CONFIG_COMMANDS_PRE(
  [AC_MSG_CHECKING([that generated files are newer than configure])
   if test -n "$am_sleep_pid"; then
     # Hide warnings about reused PIDs.
     wait $am_sleep_pid 2>/dev/null
   fi
   AC_MSG_RESULT([done])])
rm -f conftest.file
])

# Copyright (C) 2009-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SILENT_RULES([DEFAULT])
# --------------------------
# Enable less verbose build rules; with the default set to DEFAULT
# ("yes" being less verbose, "no" or empty being verbose).
AC_DEFUN([AM_SILENT_RULES],
[AC_ARG_ENABLE([silent-rules], [dnl
AS_HELP_STRING(
  [--enable-silent-rules],
  [less verbose build output (undo: "make V=1")])
AS_HELP_STRING(
  [--disable-silent-rules],
  [verbose build output (undo: "make V=0")])dnl
])
case $enable_silent_rules in @%:@ (((
  yes) AM_DEFAULT_VERBOSITY=0;;
   no) AM_DEFAULT_VERBOSITY=1;;
    *) AM_DEFAULT_VERBOSITY=m4_if([$1], [yes], [0], [1]);;
esac
dnl
dnl A few 'make' implementations (e.g., NonStop OS and NextStep)
dnl do not support nested variable expansions.
dnl See automake bug#9928 and bug#10237.
am_make=${MAKE-make}
AC_CACHE_CHECK([whether $am_make supports nested variables],
   [am_cv_make_support_nested_variables],
   [if AS_ECHO([['TRUE=$(BAR$(V))
BAR0=false
BAR1=true
V=1
am__doit:
	@$(TRUE)
.PHONY: am__doit']]) | $am_make -f - >/dev/null 2>&1; then
  am_cv_make_support_nested_variables=yes
else
  am_cv_make_support_nested_variables=no
fi])
if test $am_cv_make_support_nested_variables = yes; then
  dnl Using '$V' instead of '$(V)' breaks IRIX make.
  AM_V='$(V)'
  AM_DEFAULT_V='$(AM_DEFAULT_VERBOSITY)'
else
  AM_V=$AM_DEFAULT_VERBOSITY
  AM_DEFAULT_V=$AM_DEFAULT_VERBOSITY
fi
AC_SUBST([AM_V])dnl
AM_SUBST_NOTMAKE([AM_V])dnl
AC_SUBST([AM_DEFAULT_V])dnl
AM_SUBST_NOTMAKE([AM_DEFAULT_V])dnl
AC_SUBST([AM_DEFAULT_VERBOSITY])dnl
AM_BACKSLASH='\'
AC_SUBST([AM_BACKSLASH])dnl
_AM_SUBST_NOTMAKE([AM_BACKSLASH])dnl
])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_STRIP
# ---------------------
# One issue with vendor 'install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in "make install-strip", and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using 'strip' when the user
# run "make install-strip".  However 'strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the 'STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be 'maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# Copyright (C) 2006-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_SUBST_NOTMAKE(VARIABLE)
# ---------------------------
# Prevent Automake from outputting VARIABLE = @VARIABLE@ in Makefile.in.
# This macro is traced by Automake.
AC_DEFUN([_AM_SUBST_NOTMAKE])

# AM_SUBST_NOTMAKE(VARIABLE)
# --------------------------
# Public sister of _AM_SUBST_NOTMAKE.
AC_DEFUN([AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE($@)])

# Check how to create a tarball.                            -*- Autoconf -*-

# Copyright (C) 2004-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_TAR(FORMAT)
# --------------------
# Check how to create a tarball in format FORMAT.
# FORMAT should be one of 'v7', 'ustar', or 'pax'.
#
# Substitute a variable $(am__tar) that is a command
# writing to stdout a FORMAT-tarball containing the directory
# $tardir.
#     tardir=directory && $(am__tar) > result.tar
#
# Substitute a variable $(am__untar) that extract such
# a tarball read from stdin.
#     $(am__untar) < result.tar
#
AC_DEFUN([_AM_PROG_TAR],
[# Always define AMTAR for backward compatibility.  Yes, it's still used
# in the wild :-(  We should find a proper way to deprecate it ...
AC_SUBST([AMTAR], ['$${TAR-tar}'])

# We'll loop over all known methods to create a tar archive until one works.
_am_tools='gnutar m4_if([$1], [ustar], [plaintar]) pax cpio none'

m4_if([$1], [v7],
  [am__tar='$${TAR-tar} chof - "$$tardir"' am__untar='$${TAR-tar} xf -'],

  [m4_case([$1],
    [ustar],
     [# The POSIX 1988 'ustar' format is defined with fixed-size fields.
      # There is notably a 21 bits limit for the UID and the GID.  In fact,
      # the 'pax' utility can hang on bigger UID/GID (see automake bug#8343
      # and bug#13588).
      am_max_uid=2097151 # 2^21 - 1
      am_max_gid=$am_max_uid
      # The $UID and $GID variables are not portable, so we need to resort
      # to the POSIX-mandated id(1) utility.  Errors in the 'id' calls
      # below are definitely unexpected, so allow the users to see them
      # (that is, avoid stderr redirection).
      am_uid=`id -u || echo unknown`
      am_gid=`id -g || echo unknown`
      AC_MSG_CHECKING([whether UID '$am_uid' is supported by ustar format])
      if test $am_uid -le $am_max_uid; then
         AC_MSG_RESULT([yes])
      else
         AC_MSG_RESULT([no])
         _am_tools=none
      fi
      AC_MSG_CHECKING([whether GID '$am_gid' is supported by ustar format])
      if test $am_gid -le $am_max_gid; then
         AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        _am_tools=none
      fi],

  [pax],
    [],

  [m4_fatal([Unknown tar format])])

  AC_MSG_CHECKING([how to create a $1 tar archive])

  # Go ahead even if we have the value already cached.  We do so because we
  # need to set the values for the 'am__tar' and 'am__untar' variables.
  _am_tools=${am_cv_prog_tar_$1-$_am_tools}

  for _am_tool in $_am_tools; do
    case $_am_tool in
    gnutar)
      for _am_tar in tar gnutar gtar; do
        AM_RUN_LOG([$_am_tar --version]) && break
      done
      am__tar="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$$tardir"'
      am__tar_="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$tardir"'
      am__untar="$_am_tar -xf -"
      ;;
    plaintar)
      # Must skip GNU tar: if it does not support --format= it doesn't create
      # ustar tarball either.
      (tar --version) >/dev/null 2>&1 && continue
      am__tar='tar chf - "$$tardir"'
      am__tar_='tar chf - "$tardir"'
      am__untar='tar xf -'
      ;;
    pax)
      am__tar='pax -L -x $1 -w "$$tardir"'
      am__tar_='pax -L -x $1 -w "$tardir"'
      am__untar='pax -r'
      ;;
    cpio)
      am__tar='find "$$tardir" -print | cpio -o -H $1 -L'
      am__tar_='find "$tardir" -print | cpio -o -H $1 -L'
      am__untar='cpio -i -H $1 -d'
      ;;
    none)
      am__tar=false
      am__tar_=false
      am__untar=false
      ;;
    esac

    # If the value was cached, stop now.  We just wanted to have am__tar
    # and am__untar set.
    test -n "${am_cv_prog_tar_$1}" && break

    # tar/untar a dummy directory, and stop if the command works.
    rm -rf conftest.dir
    mkdir conftest.dir
    echo GrepMe > conftest.dir/file
    AM_RUN_LOG([tardir=conftest.dir && eval $am__tar_ >conftest.tar])
    rm -rf conftest.dir
    if test -s conftest.tar; then
      AM_RUN_LOG([$am__untar <conftest.tar])
      AM_RUN_LOG([cat conftest.dir/file])
      grep GrepMe conftest.dir/file >/dev/null 2>&1 && break
    fi
  done
  rm -rf conftest.dir

  AC_CACHE_VAL([am_cv_prog_tar_$1], [am_cv_prog_tar_$1=$_am_tool])
  AC_MSG_RESULT([$am_cv_prog_tar_$1])])

AC_SUBST([am__tar])
AC_SUBST([am__untar])
]) # _AM_PROG_TAR

# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_prefix_config_h.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PREFIX_CONFIG_H [(OUTPUT-HEADER [,PREFIX [,ORIG-HEADER]])]
#
# DESCRIPTION
#
#   Generate an installable config.h.
#
#   A package should not normally install its config.h as a system header,
#   but if it must, this macro can be used to avoid namespace pollution by
#   making a copy of config.h with a prefix added to all the macro names.
#
#   Each "#define SOMEDEF" line of the configuration header has the given
#   prefix added, in the same case as the first character of the macro name.
#
#   Defaults:
#
#     OUTPUT-HEADER = $PACKAGE-config.h
#     PREFIX = $PACKAGE
#     ORIG-HEADER, from AM_CONFIG_HEADER(config.h)
#
#   Your configure.ac script should contain both macros in this order.
#
#   Example:
#
#     AC_INIT(config.h.in)        # config.h.in as created by "autoheader"
#     AM_INIT_AUTOMAKE(testpkg, 0.1.1)    # makes #undef VERSION and PACKAGE
#     AM_CONFIG_HEADER(config.h)          # prep config.h from config.h.in
#     AX_PREFIX_CONFIG_H(mylib/_config.h) # prep mylib/_config.h from it..
#     AC_MEMORY_H                         # makes "#undef NEED_MEMORY_H"
#     AC_C_CONST_H                        # makes "#undef const"
#     AC_OUTPUT(Makefile)                 # creates the "config.h" now
#                                         # and also mylib/_config.h
#
#   If the argument to AX_PREFIX_CONFIG_H would have been omitted then the
#   default output file would have been called simply "testpkg-config.h",
#   but even under the name "mylib/_config.h" it contains prefix-defines
#   like
#
#     #ifndef TESTPKG_VERSION
#     #define TESTPKG_VERSION "0.1.1"
#     #endif
#     #ifndef TESTPKG_NEED_MEMORY_H
#     #define TESTPKG_NEED_MEMORY_H 1
#     #endif
#     #ifndef _testpkg_const
#     #define _testpkg_const _const
#     #endif
#
#   and this "mylib/_config.h" can be installed along with other header
#   files, which is most convenient when creating a shared library (that has
#   some headers) whose functionality depends on features detected at
#   compile-time. No need to invent some "mylib-confdefs.h.in" manually.
#
#   Note that some AC_DEFINEs that end up in the config.h file are actually
#   self-referential - e.g. AC_C_INLINE, AC_C_CONST, and the AC_TYPE_OFF_T
#   say that they "will define inline|const|off_t if the system does not do
#   it by itself". You might want to clean up about these - consider an
#   extra mylib/conf.h that reads something like:
#
#     #include <mylib/_config.h>
#     #ifndef _testpkg_const
#     #define _testpkg_const const
#     #endif
#
#   and then start using _testpkg_const in the header files. That is also a
#   good thing to differentiate whether some library-user has starting to
#   take up with a different compiler, so perhaps it could read something
#   like this:
#
#     #ifdef _MSC_VER
#     #include <mylib/_msvc.h>
#     #else
#     #include <mylib/_config.h>
#     #endif
#     #ifndef _testpkg_const
#     #define _testpkg_const const
#     #endif
#
# LICENSE
#
#   Copyright (c) 2014 Reuben Thomas <rrt@sc3d.org>
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2008 Marten Svantesson
#   Copyright (c) 2008 Gerald Point <Gerald.Point@labri.fr>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 15

AC_DEFUN([AX_PREFIX_CONFIG_H],[dnl
AC_PREREQ([2.62])
AC_BEFORE([AC_CONFIG_HEADERS],[$0])dnl
AC_CONFIG_COMMANDS(m4_default([$1], [$PACKAGE-config.h]),[dnl
AS_VAR_PUSHDEF([_OUT],[ac_prefix_conf_OUT])dnl
AS_VAR_PUSHDEF([_DEF],[ac_prefix_conf_DEF])dnl
AS_VAR_PUSHDEF([_PKG],[ac_prefix_conf_PKG])dnl
AS_VAR_PUSHDEF([_LOW],[ac_prefix_conf_LOW])dnl
AS_VAR_PUSHDEF([_UPP],[ac_prefix_conf_UPP])dnl
AS_VAR_PUSHDEF([_INP],[ac_prefix_conf_INP])dnl
m4_pushdef([_script],[conftest.prefix])dnl
m4_pushdef([_symbol],[m4_cr_Letters[]m4_cr_digits[]_])dnl
_OUT=`echo m4_default([$1], [$PACKAGE-config.h])`
_DEF=`echo _$_OUT | sed -e "y:m4_cr_letters:m4_cr_LETTERS[]:" -e "s/@<:@^m4_cr_Letters@:>@/_/g"`
_PKG=`echo m4_default([$2], [$PACKAGE])`
_LOW=`echo _$_PKG | sed -e "y:m4_cr_LETTERS-:m4_cr_letters[]_:"`
_UPP=`echo $_PKG | sed -e "y:m4_cr_letters-:m4_cr_LETTERS[]_:"  -e "/^@<:@m4_cr_digits@:>@/s/^/_/"`
_INP=`echo "$3" | sed -e 's/ *//'`
if test ".$_INP" = "."; then
   for ac_file in : $CONFIG_HEADERS; do test "_$ac_file" = _: && continue
     case "$ac_file" in
        *.h) _INP=$ac_file ;;
        *)
     esac
     test ".$_INP" != "." && break
   done
fi
if test ".$_INP" = "."; then
   case "$_OUT" in
      */*) _INP=`basename "$_OUT"`
      ;;
      *-*) _INP=`echo "$_OUT" | sed -e "s/@<:@_symbol@:>@*-//"`
      ;;
      *) _INP=config.h
      ;;
   esac
fi
if test -z "$_PKG" ; then
   AC_MSG_ERROR([no prefix for _PREFIX_PKG_CONFIG_H])
else
  if test ! -f "$_INP" ; then if test -f "$srcdir/$_INP" ; then
     _INP="$srcdir/$_INP"
  fi fi
  AC_MSG_NOTICE(creating $_OUT - prefix $_UPP for $_INP defines)
  if test -f $_INP ; then
    AS_ECHO(["s/^@%:@undef  *\\(@<:@m4_cr_LETTERS[]_@:>@\\)/@%:@undef $_UPP""_\\1/"]) > _script
    AS_ECHO(["s/^@%:@undef  *\\(@<:@m4_cr_letters@:>@\\)/@%:@undef $_LOW""_\\1/"]) >> _script
    AS_ECHO(["s/^@%:@def[]ine  *\\(@<:@m4_cr_LETTERS[]_@:>@@<:@_symbol@:>@*\\)\\(.*\\)/@%:@ifndef $_UPP""_\\1\\"]) >> _script
    AS_ECHO(["@%:@def[]ine $_UPP""_\\1\\2\\"]) >> _script
    AS_ECHO(["@%:@endif/"]) >> _script
    AS_ECHO(["s/^@%:@def[]ine  *\\(@<:@m4_cr_letters@:>@@<:@_symbol@:>@*\\)\\(.*\\)/@%:@ifndef $_LOW""_\\1\\"]) >> _script
    AS_ECHO(["@%:@define $_LOW""_\\1\\2\\"]) >> _script
    AS_ECHO(["@%:@endif/"]) >> _script
    # now executing _script on _DEF input to create _OUT output file
    echo "@%:@ifndef $_DEF"      >$tmp/pconfig.h
    echo "@%:@def[]ine $_DEF 1" >>$tmp/pconfig.h
    echo ' ' >>$tmp/pconfig.h
    echo /'*' $_OUT. Generated automatically at end of configure. '*'/ >>$tmp/pconfig.h

    sed -f _script $_INP >>$tmp/pconfig.h
    echo ' ' >>$tmp/pconfig.h
    echo '/* once:' $_DEF '*/' >>$tmp/pconfig.h
    echo "@%:@endif" >>$tmp/pconfig.h
    if cmp -s $_OUT $tmp/pconfig.h 2>/dev/null; then
      AC_MSG_NOTICE([$_OUT is unchanged])
    else
      ac_dir=`AS_DIRNAME(["$_OUT"])`
      AS_MKDIR_P(["$ac_dir"])
      rm -f "$_OUT"
      mv $tmp/pconfig.h "$_OUT"
    fi
  else
    AC_MSG_ERROR([input file $_INP does not exist - skip generating $_OUT])
  fi
  rm -f conftest.*
fi
m4_popdef([_symbol])dnl
m4_popdef([_script])dnl
AS_VAR_POPDEF([_INP])dnl
AS_VAR_POPDEF([_UPP])dnl
AS_VAR_POPDEF([_LOW])dnl
AS_VAR_POPDEF([_PKG])dnl
AS_VAR_POPDEF([_DEF])dnl
AS_VAR_POPDEF([_OUT])dnl
],[PACKAGE="$PACKAGE"])])

# ===========================================================================
#     http://www.gnu.org/software/autoconf-archive/ax_split_version.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_SPLIT_VERSION
#
# DESCRIPTION
#
#   Splits a version number in the format MAJOR.MINOR.POINT into its
#   separate components.
#
#   Sets the variables.
#
# LICENSE
#
#   Copyright (c) 2008 Tom Howard <tomhoward@users.sf.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 9

AC_DEFUN([AX_SPLIT_VERSION],[
    AC_REQUIRE([AC_PROG_SED])
    AX_MAJOR_VERSION=`echo "$VERSION" | $SED 's/\([[^.]][[^.]]*\).*/\1/'`
    AX_MINOR_VERSION=`echo "$VERSION" | $SED 's/[[^.]][[^.]]*.\([[^.]][[^.]]*\).*/\1/'`
    AX_POINT_VERSION=`echo "$VERSION" | $SED 's/[[^.]][[^.]]*.[[^.]][[^.]]*.\(.*\)/\1/'`
    AC_MSG_CHECKING([Major version])
    AC_MSG_RESULT([$AX_MAJOR_VERSION])
    AC_MSG_CHECKING([Minor version])
    AC_MSG_RESULT([$AX_MINOR_VERSION])
    AC_MSG_CHECKING([Point version])
    AC_MSG_RESULT([$AX_POINT_VERSION])
])

dnl This is a modified version of the Teuchos config dir from Trilinos
dnl with the following license.
dnl
dnl ***********************************************************************
dnl
dnl                    Teuchos: Common Tools Package
dnl                 Copyright (2004) Sandia Corporation
dnl
dnl Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
dnl license for use of this work by or on behalf of the U.S. Government.
dnl
dnl This library is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU Lesser General Public License as
dnl published by the Free Software Foundation; either version 2.1 of the
dnl License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
dnl USA
dnl Questions? Contact Michael A. Heroux (maherou@sandia.gov)
dnl
dnl ***********************************************************************
dnl
dnl @synopsis SC_BLAS(PREFIX, DGEMM-FUNCTION,
dnl                   [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by SC_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id: acx_blas.m4,v 1.3 2006/04/21 02:29:27 jmwille Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
dnl Edited by Jim Willenbring on 5-14-2004 to check for dgemm instead of
dnl sgemm.
dnl Edited by Jim Willenbring on 4-17-2006 to stop looking for BLAS if
dnl a specific BLAS library specified by a user cannot be used.

dnl Edited by Carsten Burstedde <carsten@ices.utexas.edu>
dnl Expect the F77_ autoconf macros to be called outside of this file.
dnl Take as argument a mangled DGEMM function to check for.
dnl This way the SC_BLAS macro can be called multiple times
dnl with different Fortran environments to minimize F77 dependencies.
dnl Replaced obsolete AC_TRY_LINK_FUNC macro.
dnl Disabled the PhiPack test since it requires BLAS_LIBS anyway.
dnl Fixed buggy generic Mac OS X library test.

dnl Subroutine to link a program using blas
dnl SC_BLAS_LINK (<added to CHECKING message>)
AC_DEFUN([SC_BLAS_LINK], [
        AC_MSG_CHECKING([for BLAS by linking a C program$1])
        AC_LINK_IFELSE([AC_LANG_PROGRAM(dnl
[[
#ifdef __cplusplus
extern "C"
void $sc_blas_func (char *, char *, int *, int *, int *, double *, double *,
                    int *, double *, int *, double *, double *, int *);
#endif
]], [[
int     i = 1;
double  alpha = 1., beta = 1.;
double  A = 1., B = 1., C = 1.;
$sc_blas_func ("N", "N", &i, &i, &i, &alpha, &A, &i, &B, &i, &beta, &C, &i);
]])],
[AC_MSG_RESULT([successful])],
[AC_MSG_RESULT([failed]); sc_blas_ok=no])
])

dnl The first argument of this macro should be the package prefix.
dnl The second argument of this macro should be a mangled DGEMM function.
AC_DEFUN([SC_BLAS], [
AC_PREREQ(2.50)
dnl Expect this to be called already.
dnl AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
dnl AC_REQUIRE([AC_F77_WRAPPERS])
sc_blas_ok=no
user_spec_blas_failed=no

AC_ARG_WITH([blas], [AS_HELP_STRING([--with-blas=<lib>],
            [change default BLAS library to <lib>
             or specify --without-blas to use no BLAS and LAPACK at all])],,
	     [withval=yes])
case $withval in
	yes | "") ;;
	no) sc_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$withval" ;;
	*) BLAS_LIBS="-l$withval" ;;
esac

dnl Expect the mangled DGEMM function name to be in $2.
sc_blas_func="$2"

sc_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test "x$sc_blas_ok" = xno; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sc_blas_func in $BLAS_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sc_blas_func])],
                       [sc_blas_ok=yes], [user_spec_blas_failed=yes])
	AC_MSG_RESULT($sc_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# If the user specified a blas library that could not be used we will
# halt the search process rather than risk finding a blas library that
# the user did not specify.

if test "x$user_spec_blas_failed" != xyes; then

# BLAS linked to by default?  (happens on some supercomputers)
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_FUNC($sc_blas_func, [sc_blas_ok=yes])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sc_blas_func,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[sc_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
# Disabled since we might want more than sgemm and dgemm.
if test "x$sc_blas_ok" = xno && false ; then
	AC_CHECK_LIB(blas, $dgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[sc_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Intel MKL library?
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(mkl, $sc_blas_func, [sc_blas_ok=yes;BLAS_LIBS="-lmkl"])
fi

# BLAS in Apple vecLib library?
if test "x$sc_blas_ok" = xno; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_CHECK_FUNC($sc_blas_func, [sc_blas_ok=yes;BLAS_LIBS="-framework vecLib"])
	LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(cxml, $sc_blas_func, [sc_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(dxml, $sc_blas_func, [sc_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test "x$sc_blas_ok" = xno; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sc_blas_func,
                                [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 sc_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(scs, $sc_blas_func, [sc_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(complib.sgimath, $sc_blas_func,
		     [sc_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(blas, $sc_blas_func,
		[AC_CHECK_LIB(essl, $sc_blas_func,
			[sc_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic Mac OS X library?
if test "x$sc_blas_ok" = xno; then
	save_LIBS="$LIBS"; LIBS="-framework Accelerate $LIBS"
	AC_CHECK_FUNC($sc_blas_func, [sc_blas_ok=yes
                               BLAS_LIBS="-framework Accelerate"])
	LIBS="$save_LIBS"
fi

# Generic BLAS library?
if test "x$sc_blas_ok" = xno; then
	AC_CHECK_LIB(blas, $sc_blas_func, [sc_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

dnl AC_SUBST(BLAS_LIBS)

fi # If the user specified library wasn't found, we skipped the remaining
   # checks.

LIBS="$sc_blas_save_LIBS"
BLAS_FLIBS=

# Test link a BLAS program
if test "x$sc_blas_ok" = xyes ; then
    dnl Link without FLIBS first
    sc_blas_save_run_LIBS="$LIBS"
    LIBS="$BLAS_LIBS $LIBS"
    SC_BLAS_LINK([ without FLIBS])
    LIBS="$sc_blas_save_run_LIBS"

    if test "x$sc_blas_ok" = xno ; then
        dnl Link with FLIBS it didn't work without
        sc_blas_save_run_LIBS="$LIBS"
        LIBS="$BLAS_LIBS $LIBS $FLIBS"
        sc_blas_ok=yes
        SC_BLAS_LINK([ with FLIBS])
        LIBS="$sc_blas_save_run_LIBS"
        BLAS_FLIBS="$FLIBS"
    fi
fi
dnl Now BLAS_FLIBS may be set or not

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$sc_blas_ok" = xyes ; then
        ifelse([$3],,
               [AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.])],[$3])
        :
elif test "x$sc_blas_ok" != xdisable ; then
        sc_blas_ok=no
        $4
fi

])


dnl sc_builtin.m4 - custom macros for distributing third-party software
dnl
dnl This file is part of the SC Library.
dnl The SC library provides support for parallel scientific applications.
dnl
dnl Copyright (C) 2008,2009 Carsten Burstedde, Lucas Wilcox.

dnl Documentation for macro names: brackets indicate optional arguments

dnl SC_BUILTIN_GETOPT_PREFIX(PREFIX)
dnl This function only activates if PREFIX_WITH_GETOPT is "yes".
dnl This function checks if getopt_long can be compiled.
dnl The shell variable PREFIX_PROVIDE_GETOPT is set to "yes" or "no".
dnl Both a define and automake conditional are set.
dnl
AC_DEFUN([SC_BUILTIN_GETOPT_PREFIX],
[
$1_PROVIDE_GETOPT="no"
AC_MSG_CHECKING([for getopt])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <getopt.h>]], [[
int oi;
struct option lo;
getopt_long (0, 0, "abc:", &lo, &oi);
]])], [AC_MSG_RESULT([successful])], [
    AC_MSG_RESULT([failed])
    AC_MSG_NOTICE([did not find getopt. Activating builtin])
    $1_PROVIDE_GETOPT="yes"
    AC_DEFINE([PROVIDE_GETOPT], 1, [Use builtin getopt])
])
AM_CONDITIONAL([$1_PROVIDE_GETOPT], [test "$$1_PROVIDE_GETOPT" = "yes"])
])
AC_DEFUN([SC_BUILTIN_GETOPT], [SC_BUILTIN_GETOPT_PREFIX([SC])])

dnl SC_BUILTIN_OBSTACK_PREFIX(PREFIX)
dnl This function only activates if PREFIX_WITH_OBSTACK is "yes".
dnl This function checks if a simple obstack program can be compiled.
dnl The shell variable PREFIX_PROVIDE_OBSTACK is set to "yes" or "no".
dnl Both a define and automake conditional are set.
dnl
AC_DEFUN([SC_BUILTIN_OBSTACK_PREFIX],
[
$1_PROVIDE_OBSTACK="no"
AC_MSG_CHECKING([for obstack])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <obstack.h>]], [[
struct obstack ob;
static void        *(*obstack_chunk_alloc) (size_t) = 0;
static void         (*obstack_chunk_free) (void *) = 0;
obstack_init (&ob);
obstack_free (&ob, 0);
]])], [AC_MSG_RESULT([successful])], [
    AC_MSG_RESULT([failed])
    AC_MSG_NOTICE([did not find obstack. Activating builtin])
    $1_PROVIDE_OBSTACK="yes"
    AC_DEFINE([PROVIDE_OBSTACK], 1, [Use builtin obstack])
])
AM_CONDITIONAL([$1_PROVIDE_OBSTACK], [test "$$1_PROVIDE_OBSTACK" = "yes"])
])
AC_DEFUN([SC_BUILTIN_OBSTACK], [SC_BUILTIN_OBSTACK_PREFIX([SC])])

dnl SC_BUILTIN_ALL_PREFIX(PREFIX)
dnl Aggregate all checks from this file for convenience.
dnl
AC_DEFUN([SC_BUILTIN_ALL_PREFIX],
[
SC_BUILTIN_GETOPT_PREFIX([$1])
SC_BUILTIN_OBSTACK_PREFIX([$1])
])
AC_DEFUN([SC_BUILTIN_ALL], [SC_BUILTIN_ALL_PREFIX([SC])])

# ===========================================================================
#            http://autoconf-archive.cryp.to/ax_c_check_flag.html
# and renamed by Carsten Burstedde <carsten@ices.utexas.edu>
# ===========================================================================
#
# SYNOPSIS
#
#   SC_C_CHECK_FLAG(FLAG-TO-CHECK,
#                   [PROLOGUE],[BODY],[ACTION-IF-SUCCESS],[ACTION-IF-FAILURE])
#
# DESCRIPTION
#
#   This macro tests if the C compiler supports the flag FLAG-TO-CHECK. If
#   successfull execute ACTION-IF-SUCCESS otherwise ACTION-IF-FAILURE.
#   PROLOGUE and BODY are optional and should be used as in AC_LANG_PROGRAM
#   macro.
#
#   This code is inspired from KDE_CHECK_COMPILER_FLAG macro. Thanks to
#   Bogdan Drozdowski <bogdandr@op.pl> for testing and bug fixes.
#
# LAST MODIFICATION
#
#   2009-02-09
#
# COPYLEFT
#
#   Copyright (c) 2008 Francesco Salvestrini <salvestrini@users.sourceforge.net>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([SC_C_CHECK_FLAG],[
  AC_PREREQ([2.61])
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_SED])

  flag=`echo "$1" | $SED 'y% .=/+-(){}<>:*,%_______________%'`

  AC_CACHE_CHECK([whether the C compiler accepts the $1 flag],
    [sc_cv_c_check_flag_$flag],[

    AC_LANG_PUSH([C])

    save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $1"
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([$2],[$3])
    ],[
      eval "sc_cv_c_check_flag_$flag=yes"
    ],[
      eval "sc_cv_c_check_flag_$flag=no"
    ])

    CFLAGS="$save_CFLAGS"

    AC_LANG_POP

  ])

  AS_IF([eval "test \"`echo '$sc_cv_c_check_flag_'$flag`\" = yes"],[
    :
    $4
  ],[
    :
    $5
  ])
])

# ===========================================================================
#            From: http://autoconf-archive.cryp.to/ax_gcc_version.html
# and renamed by Carsten Burstedde <carsten@ices.utexas.edu>
# ===========================================================================
#
# SYNOPSIS
#
#   SC_C_VERSION  (Extension of AX_GCC_VERSION to more C compilers)
#
# DESCRIPTION
#
#   This macro retrieves the cc version and returns it in the C_VERSION
#   variable if available, an empty string otherwise.
#
# LAST MODIFICATION
#
#   2009-02-09
#
# COPYLEFT
#
#   Copyright (c) 2008 Lucas Wilcox <lucasw@ices.utexas.edu>
#   Copyright (c) 2008 Francesco Salvestrini <salvestrini@users.sourceforge.net>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([SC_C_VERSION], [
  C_VERSION=""
  AS_IF([test "x$C_VERSION" = "x"],[
    SC_C_CHECK_FLAG([-V],[],[],[
      sc_pgcc_version_option=yes
    ],[
      sc_pgcc_version_option=no
    ])
    AS_IF([test "x$sc_pgcc_version_option" != "xno"],[
      AC_CACHE_CHECK([pgcc version],[sc_cv_pgcc_version],[
        # The sed part removes all new lines
        sc_cv_pgcc_version="`$CC -V 2>/dev/null | sed -e :a -e '$!N; s/\n/ /; ta'`"
        AS_IF([test "x$sc_cv_pgcc_version" = "x"],[
          sc_cv_pgcc_version=""
          ])
        ])
      C_VERSION=$sc_cv_pgcc_version
    ])
  ])

  AS_IF([test "x$C_VERSION" = "x"],[
    SC_C_CHECK_FLAG([-dumpversion],[],[],[
      sc_gcc_version_option=yes
    ],[
      sc_gcc_version_option=no
    ])
    AS_IF([test "x$sc_gcc_version_option" != "xno"],[
      AC_CACHE_CHECK([gcc version],[sc_cv_gcc_version],[
        # The sed part removes all new lines
        sc_cv_gcc_version="`$CC -dumpversion | sed -e :a -e '$!N; s/\n/ /; ta'`"
        AS_IF([test "x$sc_cv_gcc_version" = "x"],[
          sc_cv_gcc_version=""
        ])
      ])
      C_VERSION=$sc_cv_gcc_version
    ])
  ])

  AC_SUBST([C_VERSION])
])


dnl sc_include.m4 - general custom macros
dnl
dnl This file is part of the SC Library.
dnl The SC library provides support for parallel scientific applications.
dnl
dnl Copyright (C) 2008,2009 Carsten Burstedde, Lucas Wilcox.

dnl Documentation for macro names: brackets indicate optional arguments

dnl SC_VERSION(PREFIX)
dnl Expose major, minor, and point version numbers as CPP defines.
dnl Also creates a makefile variable PACKAGE_PREFIX with value PREFIX.
dnl
AC_DEFUN([SC_VERSION],
[
  AX_SPLIT_VERSION
  AC_DEFINE_UNQUOTED([VERSION_MAJOR],[$AX_MAJOR_VERSION],[Package major version])
  AC_DEFINE_UNQUOTED([VERSION_MINOR],[$AX_MINOR_VERSION],[Package minor version])
  AC_DEFINE_UNQUOTED([VERSION_POINT],[$AX_POINT_VERSION],[Package point version])
  AC_SUBST([PACKAGE_PREFIX], [$1])
])

dnl SC_ARG_ENABLE_PREFIX(NAME, COMMENT, TOKEN, PREFIX, HELPEXTRA)
dnl Check for --enable/disable-NAME using shell variable PREFIX_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional PREFIX_ENABLE_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([SC_ARG_ENABLE_PREFIX],
[
AC_ARG_ENABLE([$1],
              [AS_HELP_STRING([--enable-$1$5], [$2])],,
              [enableval=no])
if test "x$enableval" != xno ; then
  AC_DEFINE([$3], 1, [DEPRECATED (use $4_ENABLE_$3 instead)])
  AC_DEFINE([ENABLE_$3], 1, [$2])
fi
AM_CONDITIONAL([$4_ENABLE_$3], [test "x$enableval" != xno])
$4_ENABLE_$3="$enableval"
])
AC_DEFUN([SC_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [SC], [$4])])

dnl SC_ARG_DISABLE_PREFIX(NAME, COMMENT, TOKEN, PREFIX, HELPEXTRA)
dnl Check for --enable/disable-NAME using shell variable PREFIX_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional PREFIX_ENABLE_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([SC_ARG_DISABLE_PREFIX],
[
AC_ARG_ENABLE([$1],
              [AS_HELP_STRING([--disable-$1$5], [$2])],,
              [enableval=yes])
if test "x$enableval" != xno ; then
  AC_DEFINE([$3], 1, [DEPRECATED (use $4_ENABLE_$3 instead)])
  AC_DEFINE([ENABLE_$3], 1, [Undefine if: $2])
fi
AM_CONDITIONAL([$4_ENABLE_$3], [test "x$enableval" != xno])
$4_ENABLE_$3="$enableval"
])
AC_DEFUN([SC_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [SC], [$4])])

dnl SC_ARG_WITH_PREFIX(NAME, COMMENT, TOKEN, PREFIX, HELPEXTRA)
dnl Check for --with/without-NAME using shell variable PREFIX_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional PREFIX_WITH_TOKEN
dnl Default is without
dnl
AC_DEFUN([SC_ARG_WITH_PREFIX],
[
AC_ARG_WITH([$1],
            [AS_HELP_STRING([--with-$1$5], [$2])],,
            [withval=no])
if test "x$withval" != xno ; then
  AC_DEFINE([$3], 1, [DEPRECATED (use $4_WITH_$3 instead)])
  AC_DEFINE([WITH_$3], 1, [$2])
fi
AM_CONDITIONAL([$4_WITH_$3], [test "x$withval" != xno])
$4_WITH_$3="$withval"
])
AC_DEFUN([SC_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [SC], [$4])])

dnl SC_ARG_WITHOUT_PREFIX(NAME, COMMENT, TOKEN, PREFIX, HELPEXTRA)
dnl Check for --with/without-NAME using shell variable PREFIX_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional PREFIX_WITH_TOKEN
dnl Default is with
dnl
AC_DEFUN([SC_ARG_WITHOUT_PREFIX],
[
AC_ARG_WITH([$1],
            [AS_HELP_STRING([--without-$1$5], [$2])],,
            [withval=yes])
if test "x$withval" != xno ; then
  AC_DEFINE([$3], 1, [DEPRECATED (use $4_WITH_$3 instead)])
  AC_DEFINE([WITH_$3], 1, [Undefine if: $2])
fi
AM_CONDITIONAL([$4_WITH_$3], [test "x$withval" != xno])
$4_WITH_$3="$withval"
])
AC_DEFUN([SC_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [SC], [$4])])

dnl SC_REQUIRE_LIB(LIBRARY LIST, FUNCTION)
dnl Check for FUNCTION in LIBRARY, exit with error if not found
dnl
AC_DEFUN([SC_REQUIRE_LIB],
    [AC_SEARCH_LIBS([$2], [$1],,
      [AC_MSG_ERROR([Could not find function $2 in $1])])])

dnl SC_CHECK_LIB(LIBRARY LIST, FUNCTION, TOKEN, PREFIX)
dnl Check for FUNCTION first as is, then in each of the libraries.
dnl Set shell variable PREFIX_HAVE_TOKEN to nonempty if found.
dnl Call AM_CONDITIONAL with PREFIX_HAVE_TOKEN:
dnl   This means it must be called outside of any shell if block.
dnl Call AC_DEFINE with HAVE_TOKEN if found.
AC_DEFUN([SC_CHECK_LIB], [
AC_SEARCH_LIBS([$2], [$1])
AM_CONDITIONAL([$4_HAVE_$3], [test "x$ac_cv_search_$2" != xno])
$4_HAVE_$3=
if test "x$ac_cv_search_$2" != xno ; then
AC_DEFINE([HAVE_$3], [1], [Have we found function $2.])
$4_HAVE_$3=yes
fi
])

dnl SC_CHECK_LIB_NOCOND(LIBRARY LIST, FUNCTION, TOKEN, PREFIX)
dnl Check for FUNCTION first as is, then in each of the libraries.
dnl Set shell variable PREFIX_HAVE_TOKEN to nonempty if found.
dnl Does not establish an automake conditional;
dnl   safe to call from within a shell if block.
dnl Call AC_DEFINE with HAVE_TOKEN if found.
AC_DEFUN([SC_CHECK_LIB_NOCOND], [
AC_SEARCH_LIBS([$2], [$1])
$4_HAVE_$3=
if test "x$ac_cv_search_$2" != xno ; then
AC_DEFINE([HAVE_$3], [1], [Have we found function $2.])
$4_HAVE_$3=yes
fi
])

dnl SC_REQUIRE_FUNCS(FUNCTION LIST)
dnl Check for all functions in FUNCTION LIST, exit with error if not found
dnl
AC_DEFUN([SC_REQUIRE_FUNCS],
[
m4_foreach_w([sc_thefunc], [$1],
             [AC_CHECK_FUNC([sc_thefunc], ,
                            [AC_MSG_ERROR([\
Could not find function sc_thefunc])])])
])

dnl SC_DETERMINE_INSTALL(PREFIX)
dnl This function throws an error if the variable PREFIX_DIR does not exist.
dnl Looks for PREFIX_DIR/{include,lib,bin} to determine installation status.
dnl Set the shell variable PREFIX_INSTALL to "yes" or "no".
dnl
AC_DEFUN([SC_DETERMINE_INSTALL],
[
if test ! -d "$$1_DIR" ; then
  AC_MSG_ERROR([Directory "$$1_DIR" does not exist])
fi
if test -d "$$1_DIR/include" || test -d "$$1_DIR/lib" || \
   test -d "$$1_DIR/bin" || test -d "$$1_DIR/share/aclocal" ; then
  $1_INSTALL=yes
else
  $1_INSTALL=no
fi
])

dnl SC_DETERMINE_INCLUDE_PATH(PREFIX, CPPFLAGS)
dnl This function expects the variable PREFIX_DIR to exist.
dnl Looks for PREFIX_DIR/include and then PREFIX_DIR/src.
dnl If neither is found, throws an error.
dnl Otherwise, set the shell variable PREFIX_CPPFLAGS to -I<dir> CPPFLAGS.
dnl
AC_DEFUN([SC_DETERMINE_INCLUDE_PATH],
[
$1_INC="$$1_DIR/include"
if test ! -d "$$1_INC" ; then
  $1_INC="$$1_DIR/src"
fi
if test ! -d "$$1_INC" ; then
  AC_MSG_ERROR([Include directories based on $$1_DIR not found])
fi
$1_CPPFLAGS="-I$$1_INC $2"
])

dnl SC_DETERMINE_LIBRARY_PATH(PREFIX, LIBS)
dnl This function expects the variable PREFIX_DIR to exist.
dnl Looks for PREFIX_DIR/lib and then PREFIX_DIR/src.
dnl If neither is found, throws an error.
dnl Otherwise, set the shell variable PREFIX_LDADD to -L<dir> LIBS.
dnl
AC_DEFUN([SC_DETERMINE_LIBRARY_PATH],
[
$1_LIB="$$1_DIR/lib"
if test ! -d "$$1_LIB" ; then
  $1_LIB="$$1_DIR/src"
fi
if test ! -d "$$1_LIB" ; then
  AC_MSG_ERROR([Library directories based on $$1_DIR not found])
fi
$1_LDADD="-L$$1_LIB $2"
])

dnl SC_DETERMINE_CONFIG_PATH(PREFIX)
dnl This function expects the variable PREFIX_DIR to exist.
dnl Looks for PREFIX_DIR/share/aclocal and then PREFIX_DIR/src.
dnl If neither is found, throws an error.
dnl Sets shell variables PREFIX_CONFIG and PREFIX_AMFLAGS.
dnl
AC_DEFUN([SC_DETERMINE_CONFIG_PATH],
[
$1_CONFIG="$$1_DIR/share/aclocal"
if test ! -d "$$1_CONFIG" ; then
  $1_CONFIG="$$1_DIR/config"
fi
if test ! -d "$$1_CONFIG" ; then
  AC_MSG_ERROR([Config directories based on $$1_DIR not found])
fi
$1_AMFLAGS="-I $$1_CONFIG"
])

dnl SC_CHECK_BLAS_LAPACK(PREFIX)
dnl This function uses the macros SC_BLAS and SC_LAPACK.
dnl It requires previous configure macros for F77 support,
dnl which are called by SC_MPI_CONFIG/SC_MPI_ENGAGE.
dnl
AC_DEFUN([SC_CHECK_BLAS_LAPACK],
[

dgemm=;AC_F77_FUNC(dgemm)
if test "x$dgemm" = xunknown ; then dgemm=dgemm_ ; fi

AC_MSG_NOTICE([Checking BLAS])
SC_BLAS([$1], [$dgemm],
        [AC_DEFINE([WITH_BLAS], 1, [Define to 1 if BLAS is used])],
        [AC_MSG_ERROR([[\
Cannot find BLAS library, specify a path using LIBS=-L<DIR> (ex.\
 LIBS=-L/usr/path/lib) or a specific library using BLAS_LIBS=DIR/LIB\
 (for example BLAS_LIBS=/usr/path/lib/libcxml.a)]])])

# at this point $sc_blas_ok is either of: yes disable
if test "x$sc_blas_ok" = xdisable ; then
        AC_MSG_NOTICE([Not using BLAS])
fi
AM_CONDITIONAL([$1_WITH_BLAS], [test "x$sc_blas_ok" = xyes])

dgecon=;AC_F77_FUNC(dgecon)
if test "x$dgecon" = xunknown ; then dgecon=dgecon_ ; fi

AC_MSG_NOTICE([Checking LAPACK])
SC_LAPACK([$1], [$dgecon],
          [AC_DEFINE([WITH_LAPACK], 1, [Define to 1 if LAPACK is used])],
          [AC_MSG_ERROR([[\
Cannot find LAPACK library, specify a path using LIBS=-L<DIR> (ex.\
 LIBS=-L/usr/path/lib) or a specific library using LAPACK_LIBS=DIR/LIB\
 (for example LAPACK_LIBS=/usr/path/lib/libcxml.a)]])])

# at this point $sc_lapack_ok is either of: yes disable
if test "x$sc_lapack_ok" = xdisable ; then
        AC_MSG_NOTICE([Not using LAPACK])
fi
AM_CONDITIONAL([$1_WITH_LAPACK], [test "x$sc_lapack_ok" = xyes])

# Append the necessary blas/lapack and fortran libraries to LIBS
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $LAPACK_FLIBS $BLAS_FLIBS"
])

dnl SC_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by libsc.  It can be used by other packages that
dnl link to libsc to add appropriate options to LIBS.
dnl
AC_DEFUN([SC_CHECK_LIBRARIES],
[
SC_REQUIRE_LIB([m], [fabs])
SC_CHECK_LIB([z], [adler32_combine], [ZLIB], [$1])
SC_CHECK_LIB([lua52 lua5.2 lua51 lua5.1 lua lua5], [lua_createtable],
	     [LUA], [$1])
SC_CHECK_BLAS_LAPACK([$1])
SC_BUILTIN_ALL_PREFIX([$1])
SC_CHECK_PTHREAD([$1])
SC_CHECK_OPENMP([$1])
SC_CHECK_MEMALIGN([$1])
dnl SC_CUDA([$1])
])

dnl SC_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using libsc as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if sc is make install'd.
dnl
AC_DEFUN([SC_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [SC], [sc])])

dnl SC_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([SC_FINAL_MESSAGES],
[
if test "x$$1_HAVE_ZLIB" = x ; then
AC_MSG_NOTICE([- $1 -------------------------------------------------
We did not find a recent zlib containing the function adler32_combine.
This is OK if the following does not matter to you:
Calling any sc functions that rely on zlib will abort your program.
These functions include sc_array_checksum and sc_vtk_write_compressed.
You can fix this by compiling a working zlib and pointing LIBS to it.
])
fi
if test "x$$1_HAVE_LUA" = x ; then
AC_MSG_NOTICE([- $1 -------------------------------------------------
We did not find a recent lua containing the function lua_createtable.
This is OK if the following does not matter to you:
Including sc_lua.h in your code will abort the compilation.
You can fix this by compiling a working lua and pointing LIBS to it.
])
fi
])

dnl This is a modified version of the Teuchos config dir from Trilinos
dnl with the following license.
dnl
dnl ***********************************************************************
dnl
dnl                    Teuchos: Common Tools Package
dnl                 Copyright (2004) Sandia Corporation
dnl
dnl Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
dnl license for use of this work by or on behalf of the U.S. Government.
dnl
dnl This library is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU Lesser General Public License as
dnl published by the Free Software Foundation; either version 2.1 of the
dnl License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
dnl USA
dnl Questions? Contact Michael A. Heroux (maherou@sandia.gov)
dnl
dnl ***********************************************************************
dnl
dnl @synopsis SC_LAPACK(PREFIX, DGECON_FUNCTION,
dnl                     [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the SC_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by SC_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id: acx_lapack.m4,v 1.3 2006/04/21 02:29:27 jmwille Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl edited by Jim Willenbring <jmwille@sandia.gov> to check for sgecon
dnl rather than cheev because by default (as of 8-13-2002) Trilinos
dnl does not build the complex portions of the lapack library.  Edited
dnl again on 5-13-2004 to check for dgecon instead of sgecon.
dnl Edited by Jim Willenbring on 4-17-2006 to stop looking for LAPACK if
dnl a specific LAPACK library specified by a user cannot be used.

dnl Edited by Carsten Burstedde <carsten@ices.utexas.edu>
dnl Expect the F77_ autoconf macros to be called outside of this file.
dnl Take as argument a mangled DGECON function to check for.
dnl This way the SC_LAPACK macro can be called multiple times
dnl with different Fortran environments to minimize F77 dependencies.
dnl Replaced obsolete AC_TRY_LINK_FUNC macro.

dnl Subroutine to link a program using lapack
dnl SC_LAPACK_LINK (<added to CHECKING message>)
AC_DEFUN([SC_LAPACK_LINK], [
        AC_MSG_CHECKING([for LAPACK by linking$1])
        AC_LINK_IFELSE([AC_LANG_PROGRAM(dnl
[[
#ifdef __cplusplus
extern "C"
void $sc_lapack_func (char *, int *, double *, int *, double *,
                      double *, double *, int *, int *);
#endif
]], [[
int     i = 1, info = 0, iwork[1];
double  anorm = 1., rcond;
double  A = 1., work[4];
$sc_lapack_func ("1", &i, &A, &i, &anorm, &rcond, work, iwork, &info);
]])],
[AC_MSG_RESULT([successful])],
[AC_MSG_RESULT([failed]); sc_lapack_ok=no])
])

dnl The first argument of this macro should be the package prefix.
dnl The second argument of this macro should be a mangled DGECON function.
AC_DEFUN([SC_LAPACK], [
AC_REQUIRE([SC_BLAS])
sc_lapack_ok=no
user_spec_lapack_failed=no

AC_ARG_WITH([lapack], [AS_HELP_STRING([--with-lapack=<lib>],
            [change default LAPACK library to <lib>
             or specify --without-lapack to use no LAPACK at all])],,
	     [withval=yes])
case $withval in
        yes | "") ;;
        no) sc_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$withval" ;;
        *) LAPACK_LIBS="-l$withval" ;;
esac

dnl Expect the mangled DGECON function name to be in $2.
sc_lapack_func="$2"

# We cannot use LAPACK if BLAS is not found
if test "x$sc_blas_ok" = xdisable ; then
        sc_lapack_ok=disable
elif test "x$sc_blas_ok" != xyes; then
        sc_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$sc_lapack_ok" = xno; then
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $sc_lapack_func in $LAPACK_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$sc_lapack_func])],
                       [sc_lapack_ok=yes], [user_spec_lapack_failed=yes])
        AC_MSG_RESULT($sc_lapack_ok)
        LIBS="$save_LIBS"
        if test "x$sc_lapack_ok" = xno; then
                LAPACK_LIBS=""
        fi
fi
fi

# If the user specified a LAPACK library that could not be used we will
# halt the search process rather than risk finding a LAPACK library that
# the user did not specify.

if test "x$user_spec_lapack_failed" != xyes; then

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test "x$sc_lapack_ok" = xno; then
        save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
        AC_CHECK_FUNC($sc_lapack_func, [sc_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test "x$sc_lapack_ok" = xno; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $sc_lapack_func,
                    [sc_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

dnl AC_SUBST(LAPACK_LIBS)

fi # If the user specified library wasn't found, we skipped the remaining
   # checks.

LAPACK_FLIBS=

# Test link and run a LAPACK program
if test "x$sc_lapack_ok" = xyes ; then
    dnl Link without FLIBS, or with FLIBS if required by BLAS
    sc_lapack_save_run_LIBS="$LIBS"
    LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $BLAS_FLIBS"
    SC_LAPACK_LINK([ w/ BLAS_FLIBS but w/o FLIBS])
    LIBS="$sc_lapack_save_run_LIBS"

    if test "x$sc_lapack_ok" = xno && test "x$BLAS_FLIBS" = x ; then
        dnl Link with FLIBS it didn't work without
        sc_lapack_save_run_LIBS="$LIBS"
        LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        sc_lapack_ok=yes
        SC_LAPACK_LINK([ with FLIBS])
        LIBS="$sc_lapack_save_run_LIBS"
        LAPACK_FLIBS="$FLIBS"
    fi
fi
dnl Now at most one of BLAS_FLIBS and LAPACK_FLIBS may be set, but not both

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$sc_lapack_ok" = xyes; then
        ifelse([$3],,
               [AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.])],[$3])
        :
elif test "x$sc_lapack_ok" != xdisable ; then
        sc_lapack_ok=no
        $4
fi

])

dnl
dnl Copyright 2005-2006 Sun Microsystems, Inc.  All rights reserved.
dnl
dnl Permission is hereby granted, free of charge, to any person obtaining a
dnl copy of this software and associated documentation files (the
dnl "Software"), to deal in the Software without restriction, including
dnl without limitation the rights to use, copy, modify, merge, publish,
dnl distribute, and/or sell copies of the Software, and to permit persons
dnl to whom the Software is furnished to do so, provided that the above
dnl copyright notice(s) and this permission notice appear in all copies of
dnl the Software and that both the above copyright notice(s) and this
dnl permission notice appear in supporting documentation.
dnl
dnl THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
dnl OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
dnl MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT
dnl OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
dnl HOLDERS INCLUDED IN THIS NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL
dnl INDIRECT OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING
dnl FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
dnl NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
dnl WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
dnl
dnl Except as contained in this notice, the name of a copyright holder
dnl shall not be used in advertising or otherwise to promote the sale, use
dnl or other dealings in this Software without prior written authorization
dnl of the copyright holder.
dnl
dnl Renamed by Carsten Burstedde <carsten@ices.utexas.edu>

# SC_PROG_LINT()
# ----------------
# Minimum version: 1.1.0
#
# Sets up flags for source checkers such as lint and sparse if --with-lint
# is specified.   (Use --with-lint=sparse for sparse.)
# Sets $LINT to name of source checker passed with --with-lint (default: splint)
# Sets $LINT_FLAGS to flags to pass to source checker
# Sets LINT automake conditional if enabled (default: disabled)
#
# Note that MPI_INCLUDE_PATH should be defined before this function is called.
#
AC_DEFUN([SC_PROG_LINT],[

# Allow checking code with lint, sparse, etc.
AC_ARG_WITH([lint], [AS_HELP_STRING([--with-lint],
            [use static source code checker (default: splint)])],
            [use_lint=$withval], [use_lint=yes])
if test "$use_lint" = yes ; then
  use_lint="splint"
fi
if test "$use_lint" != no ; then
  AC_PATH_PROG([LINT], [$use_lint], [no])
  if test "$LINT" = no ; then
    AC_MSG_WARN([Static source code checker $use_lint not found])
    use_lint="no"
  fi
fi

if test "$use_lint" != no ; then

if test "$LINT_FLAGS" = "" ; then
    case $LINT in
      lint|*/lint)
        case $host_os in
          solaris*)
            LINT_FLAGS="-u -b -h -erroff=E_INDISTING_FROM_TRUNC2"
            ;;
        esac
        ;;
      splint|*/splint)
        LINT_FLAGS="-weak -fixedformalarray -badflag -preproc -unixlib"
        ;;
    esac
fi

case $LINT in
  splint|*/splint)
    LINT_FLAGS="$LINT_FLAGS -DSC_SPLINT \
                -systemdirs /usr/include:$MPI_INCLUDE_PATH"
    ;;
esac

fi

AC_SUBST(LINT)
AC_SUBST(LINT_FLAGS)
AM_CONDITIONAL(LINT, [test "$use_lint" != no])

])


dnl SC_CHECK_MEMALIGN(PREFIX)
dnl Let the user specify --enable-memalign[=X] or --disable-memalign.
dnl The alignment argument X must be a multiple of sizeof (void *).
dnl
dnl The default is --enable-memalign, which sets X to (sizeof (void *)).
dnl
dnl This macro also searches for the aligned allocation functions aligned_alloc
dnl (C11 / glibc >= 2.16) and posix_memalign (POSIX / glibc >= 2.1.91) and
dnl defines SC_HAVE_ALIGNED_ALLOC and SC_HAVE_POSIX_MEMALIGN, respectively.
dnl If found and alignment is enabled, this macro runs the link tests.
dnl
dnl If memory alignment is selected, the sc_malloc calls and friends will
dnl use the aligned version, relying on posix_memalign if it exists.
dnl
AC_DEFUN([SC_CHECK_MEMALIGN], [

dnl check for size of types
AC_CHECK_SIZEOF([void *])

dnl check for presence of functions
AC_CHECK_FUNCS([aligned_alloc posix_memalign])

dnl custom memory alignment option
AC_MSG_CHECKING([for memory alignment option])
SC_ARG_DISABLE_PREFIX([memalign],
  [use aligned malloc (optionally use --enable-memalign=<bytes>)],
  [MEMALIGN], [$1])

dnl read the value of the configuration argument
if test "x$$1_ENABLE_MEMALIGN" != xno ; then
  if test "x$$1_ENABLE_MEMALIGN" != xyes ; then

    dnl make sure the alignment is a number 
    $1_MEMALIGN_BYTES=`echo "$$1_ENABLE_MEMALIGN" | tr -c -d '[[:digit:]]'`
    $1_MEMALIGN_BYTES_LINK="$$1_MEMALIGN_BYTES"
    if test "x$$1_MEMALIGN_BYTES" = x ; then
      AC_MSG_ERROR([Please provide --enable-memalign with a numeric value or nothing])
    fi
  else
    $1_MEMALIGN_BYTES_LINK="SIZEOF_VOID_P"
    $1_MEMALIGN_BYTES="$1_$$1_MEMALIGN_BYTES_LINK"
  fi
  AC_DEFINE_UNQUOTED([MEMALIGN_BYTES], [($$1_MEMALIGN_BYTES)],
                     [desired alignment of allocations in bytes])
  AC_MSG_RESULT([$$1_MEMALIGN_BYTES])

dnl verify that aligned_alloc can be linked against
  if test "x$ac_cv_func_aligned_alloc" = xyes ; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <stdlib.h>]],
[[
int *a = (int *) aligned_alloc ($$1_MEMALIGN_BYTES_LINK, 3 * sizeof(*a));
free(a);
]])],
                   [], [AC_MSG_ERROR([Linking aligned_alloc failed])])
  fi

dnl verify that posix_memalign can be linked against
  if test "x$ac_cv_func_posix_memalign" = xyes ; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#include<stdlib.h>
#include<errno.h>
]],
[[
int *a;
int err = posix_memalign((void **) &a, $$1_MEMALIGN_BYTES_LINK, 3 * sizeof(*a));
free(a);
]])],
                   [], [AC_MSG_ERROR([Linking posix_memalign failed])])
  fi

dnl the function memalign is obsolete and not used
dnl   $1_WITH_MEMALIGN="yes"
dnl   AC_SEARCH_LIBS([memalign])
dnl   if test "x$ac_cv_search_memalign" != "xnone required" ; then
dnl     $1_WITH_MEMALIGN=no
dnl   fi
dnl   if test "x$$1_WITH_MEMALIGN" = xyes ; then
dnl     AC_LINK_IFELSE([AC_LANG_PROGRAM(
dnl [[
dnl #include<stdlib.h>
dnl ]],
dnl [[
dnl int *a = (int *) memalign($$1_MEMALIGN_BYTES,3*sizeof(*a));
dnl free(a);
dnl ]])],
dnl                    [],[$1_WITH_MEMALIGN=no])
dnl   fi
dnl   if test "x$$1_WITH_MEMALIGN" = xyes ; then
dnl     AC_DEFINE([WITH_MEMALIGN],[1],[define to 1 if memalign() found])
dnl   fi
else
  AC_MSG_RESULT([not used])
fi
])

dnl
dnl SC_MPI_CONFIG(PREFIX, , )
dnl
dnl If the second argument is nonempty, also includes configuration for F77 and FC.
dnl If the third argument is nonempty, also includes configuration for CXX.
dnl
dnl Checks the configure options
dnl --enable-mpi      If enabled and CC is not set, do export CC=mpicc.
dnl                   This may be "too late" if AC_PROG_CC was called earlier.
dnl                   In that case you need to set CC=mpicc (or other compiler)
dnl                   on the configure command line.
dnl                   Likewise for F77, FC and CXX if enabled in SC_MPI_CONFIG.
dnl --disable-mpiio   Only effective if --enable-mpi is given.  In this case,
dnl                   do not use MPI I/O in sc and skip the compile-and-link test.
dnl --disable-mpithread Only effective if --enable-mpi is given.  In this case,
dnl                   do not use MPI_Init_thread () and skip compile-and-link test.
dnl
dnl If MPI is enabled, set AC_DEFINE and AC_CONDITIONAL for PREFIX_ENABLE_MPI.
dnl If MPI I/O is not disabled, set these for PREFIX_ENABLE_MPIIO.
dnl If MPI_Init_thread is not disabled, set these for PREFIX_ENABLE_MPITHREAD.
dnl
dnl SC_MPI_ENGAGE(PREFIX)
dnl
dnl Relies on SC_MPI_CONFIG to be called before.
dnl Calls AC_PROG_CC and other macros related to the C compiler.
dnl Calls AC_PROG_F77 and others if F77 is enabled in SC_MPI_CONFIG.
dnl Calls AC_PROG_FC and others if FC is enabled in SC_MPI_CONFIG.
dnl Calls AC_PROG_CXX and others if CXX is enabled in SC_MPI_CONFIG.
dnl
dnl If MPI is enabled, a compile-and-link test is performed.  It aborts
dnl configuration on failure.
dnl If MPI is enabled and I/O is not disabled, a compile-and-link test
dnl for MPI I/O is performed.  It aborts configuration on failure.
dnl If MPI is enabled and MPITHREAD is not disabled, a compile-and-link test
dnl for MPI_Init_thread is performed.  It aborts configuration on failure.
dnl
dnl These macros are separate because of the AC_REQUIRE logic inside autoconf.

AC_DEFUN([SC_MPI_CONFIG],
[
HAVE_PKG_MPI=no
HAVE_PKG_MPIIO=no
HAVE_PKG_MPITHREAD=no
m4_ifval([$2], [m4_define([SC_CHECK_MPI_F77], [yes])])
m4_ifval([$2], [m4_define([SC_CHECK_MPI_FC], [yes])])
m4_ifval([$3], [m4_define([SC_CHECK_MPI_CXX], [yes])])

dnl The shell variable SC_ENABLE_MPI is set if --enable-mpi is given.
dnl Therefore all further checking uses the HAVE_PKG_MPI shell variable
dnl and neither AC_DEFINE nor AM_CONDITIONAL are invoked at this point.
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
               [enable MPI (force serial code otherwise)])],,
              [enableval=no])
if test "x$enableval" = xyes ; then
  HAVE_PKG_MPI=yes
elif test "x$enableval" != xno ; then
  AC_MSG_ERROR([Please use --enable-mpi without an argument])
fi
AC_MSG_CHECKING([whether we are using MPI])
AC_MSG_RESULT([$HAVE_PKG_MPI])

dnl The shell variable SC_ENABLE_MPIIO is set if --disable-mpiio is not given.
dnl If not disabled, MPI I/O will be verified by a compile/link test below.
AC_ARG_ENABLE([mpiio],
              [AS_HELP_STRING([--disable-mpiio],
               [do not use MPI I/O (even if MPI is enabled)])],,
              [enableval=yes])
if test "x$enableval" = xyes ; then
  if test "x$HAVE_PKG_MPI" = xyes ; then
    HAVE_PKG_MPIIO=yes
  fi
elif test "x$enableval" != xno ; then
  AC_MSG_ERROR([Please don't use --enable-mpiio; it's the default now])
fi
AC_MSG_CHECKING([whether we are using MPI I/O])
AC_MSG_RESULT([$HAVE_PKG_MPIIO])

dnl The variable SC_ENABLE_MPITHREAD is set if --disable-mpithread not given.
dnl If not disabled, MPI_Init_thread will be verified by a compile/link test.
AC_ARG_ENABLE([mpithread],
              [AS_HELP_STRING([--disable-mpithread],
               [do not use MPI_Init_thread (even if MPI is enabled)])],,
              [enableval=yes])
if test "x$enableval" = xyes ; then
  if test "x$HAVE_PKG_MPI" = xyes ; then
    HAVE_PKG_MPITHREAD=yes
  fi
elif test "x$enableval" != xno ; then
  AC_MSG_ERROR([Please don't use --enable-mpithread; it's the default now])
fi
AC_MSG_CHECKING([whether we are using MPI_Init_thread])
AC_MSG_RESULT([$HAVE_PKG_MPITHREAD])

dnl Establish the MPI test environment
$1_MPIRUN=
$1_MPI_TEST_FLAGS=
if test "x$HAVE_PKG_MPI" = xyes ; then
AC_CHECK_PROGS([$1_MPIRUN], [mpiexec mpirun])
if test "x$$1_MPIRUN" = xmpiexec ; then
  # $1_MPIRUN=mpiexec
  $1_MPI_TEST_FLAGS="-n 2"
elif test "x$$1_MPIRUN" = xmpirun ; then
  # $1_MPIRUN=mpirun
  $1_MPI_TEST_FLAGS="-np 2"
else
  $1_MPIRUN=
fi
AC_SUBST([$1_MPIRUN])
AC_SUBST([$1_MPI_TEST_FLAGS])
fi
AM_CONDITIONAL([$1_MPIRUN], [test "x$$1_MPIRUN" != x])

dnl Set compilers if not already set and set define and conditionals
if test "x$HAVE_PKG_MPI" = xyes ; then
m4_ifset([SC_CHECK_MPI_F77], [
  if test "x$F77" = x ; then
    export F77=mpif77
  fi
  AC_MSG_NOTICE([                            F77 set to $F77])
])
m4_ifset([SC_CHECK_MPI_FC], [
  if test "x$FC" = x ; then
    export FC=mpif90
  fi
  AC_MSG_NOTICE([                             FC set to $FC])
])
  if test "x$CC" = x ; then
    export CC=mpicc
  fi
  AC_MSG_NOTICE([                             CC set to $CC])
m4_ifset([SC_CHECK_MPI_CXX], [
  if test "x$CXX" = x ; then
    export CXX=mpicxx
  fi
  AC_MSG_NOTICE([                            CXX set to $CXX])
])
  AC_DEFINE([MPI], 1, [DEPRECATED (use $1_ENABLE_MPI instead)])
  AC_DEFINE([ENABLE_MPI], 1, [Define to 1 if we are using MPI])
  if test "x$HAVE_PKG_MPIIO" = xyes ; then
    AC_DEFINE([MPIIO], 1, [DEPRECATED (use $1_ENABLE_MPIIO instead)])
    AC_DEFINE([ENABLE_MPIIO], 1, [Define to 1 if we are using MPI I/O])
  fi
  if test "x$HAVE_PKG_MPITHREAD" = xyes ; then
    AC_DEFINE([ENABLE_MPITHREAD], 1, [Define to 1 if we are using MPI_Init_thread])
  fi
else
m4_ifset([SC_CHECK_MPI_F77], [
  if test "x$F77" = x ; then
    AC_CHECK_PROGS([$1_F77_COMPILER], [gfortran g77 f77 ifort])
    if test "x$$1_F77_COMPILER" != x ; then
      F77="$$1_F77_COMPILER"
    fi
  fi
], [:])
m4_ifset([SC_CHECK_MPI_FC], [
  if test "x$FC" = x ; then
    AC_CHECK_PROGS([$1_FC_COMPILER], [gfortran ifort])
    if test "x$$1_FC_COMPILER" != x ; then
      FC="$$1_FC_COMPILER"
    fi
  fi
], [:])
fi
AM_CONDITIONAL([$1_ENABLE_MPI], [test "x$HAVE_PKG_MPI" = xyes])
AM_CONDITIONAL([$1_ENABLE_MPIIO], [test "x$HAVE_PKG_MPIIO" = xyes])
AM_CONDITIONAL([$1_ENABLE_MPITHREAD], [test "x$HAVE_PKG_MPITHREAD" = xyes])
])

dnl SC_MPI_F77_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI F77 test program
dnl
dnl DEACTIVATED since it triggers a bug in autoconf:
dnl AC_LANG_PROGRAM(Fortran 77): ignoring PROLOGUE: [
dnl
dnl AC_DEFUN([SC_MPI_F77_COMPILE_AND_LINK],
dnl [
dnl AC_MSG_CHECKING([compile/link for MPI F77 program])
dnl AC_LINK_IFELSE([AC_LANG_PROGRAM(
dnl [[
dnl       include "mpif.h"
dnl ]], [[
dnl       call MPI_INIT (ierror)
dnl       call MPI_COMM_SIZE (MPI_COMM_WORLD, isize, ierror)
dnl       call MPI_COMM_RANK (MPI_COMM_WORLD, irank, ierror)
dnl       print*, isize, irank, ': Hello world'
dnl       call MPI_FINALIZE (ierror)
dnl ]])],
dnl [AC_MSG_RESULT([successful])
dnl  $1],
dnl [AC_MSG_RESULT([failed])
dnl  $2])
dnl ])

dnl SC_MPI_FC_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI FC test program
dnl
dnl DEACTIVATED since it triggers a bug in autoconf:
dnl AC_LANG_PROGRAM(Fortran): ignoring PROLOGUE: [
dnl
dnl AC_DEFUN([SC_MPI_FC_COMPILE_AND_LINK],
dnl [
dnl AC_MSG_CHECKING([compile/link for MPI FC program])
dnl AC_LINK_IFELSE([AC_LANG_PROGRAM(
dnl [[
dnl       include "mpif90.h"
dnl ]], [[
dnl       call MPI_INIT (ierror)
dnl       call MPI_COMM_SIZE (MPI_COMM_WORLD, isize, ierror)
dnl       call MPI_COMM_RANK (MPI_COMM_WORLD, irank, ierror)
dnl       print*, isize, irank, ': Hello world'
dnl       call MPI_FINALIZE (ierror)
dnl ]])],
dnl [AC_MSG_RESULT([successful])
dnl  $1],
dnl [AC_MSG_RESULT([failed])
dnl  $2])
dnl ])

dnl SC_MPI_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI C test program
dnl
AC_DEFUN([SC_MPI_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
]], [[
MPI_Init ((int *) 0, (char ***) 0);
MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPI_CXX_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI CXX test program
dnl
AC_DEFUN([SC_MPI_CXX_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI CXX program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
#include <iostream>
]], [[
std::cout << "Hello C++ MPI" << std::endl;
MPI_Init ((int *) 0, (char ***) 0);
MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPIIO_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI I/O test program
dnl
AC_DEFUN([SC_MPIIO_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI I/O C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
]], [[
MPI_File fh;
MPI_Init ((int *) 0, (char ***) 0);
MPI_File_open (MPI_COMM_WORLD, "filename",
               MPI_MODE_WRONLY | MPI_MODE_APPEND,
               MPI_INFO_NULL, &fh);
MPI_File_close (&fh);
MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPITHREAD_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI_Init_thread test program
dnl
AC_DEFUN([SC_MPITHREAD_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI_Init_thread C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
]], [[
int mpiret;
int mpithr;
mpiret = MPI_Init_thread ((int *) 0, (char ***) 0,
                          MPI_THREAD_MULTIPLE, &mpithr);
mpiret = MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPIWINSHARED_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI_Win_allocate_shared test program
dnl
AC_DEFUN([SC_MPIWINSHARED_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI_Win_allocate_shared C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
]], [[
int mpiret;
int mpithr;
int disp_unit=0;
char *baseptr;
MPI_Win win;
MPI_Init ((int *) 0, (char ***) 0);
mpiret = MPI_Win_allocate_shared(0,disp_unit,MPI_INFO_NULL,MPI_COMM_WORLD,(void *) &baseptr,&win);
mpiret = MPI_Win_shared_query(win,0,0,&disp_unit,(void *) &baseptr);
mpiret = MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,MPI_MODE_NOCHECK,win);
mpiret = MPI_Win_unlock(0,win);
mpiret = MPI_Win_free(&win);
mpiret = MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPICOMMSHARED_C_COMPILE_AND_LINK([action-if-successful], [action-if-failed])
dnl Compile and link an MPI_COMM_TYPE_SHARED test program
dnl
AC_DEFUN([SC_MPICOMMSHARED_C_COMPILE_AND_LINK],
[
AC_MSG_CHECKING([compile/link for MPI_COMM_TYPE_SHARED C program])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#undef MPI
#include <mpi.h>
]], [[
int mpiret;
MPI_Comm subcomm;
MPI_Init ((int *) 0, (char ***) 0);
mpiret = MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&subcomm);
mpiret = MPI_Finalize ();
]])],
[AC_MSG_RESULT([successful])
 $1],
[AC_MSG_RESULT([failed])
 $2])
])

dnl SC_MPI_INCLUDES
dnl Call the compiler with various --show* options
dnl to figure out the MPI_INCLUDES and MPI_INCLUDE_PATH varables
dnl
AC_DEFUN([SC_MPI_INCLUDES],
[
MPI_INCLUDES=
MPI_INCLUDE_PATH=
if test "x$HAVE_PKG_MPI" = xyes ; then
  AC_MSG_NOTICE([Trying to determine MPI_INCLUDES])
  for SHOW in -showme:compile -showme:incdirs -showme -show ; do
    if test "x$MPI_INCLUDES" = x ; then
      AC_MSG_CHECKING([$SHOW])
      if MPI_CC_RESULT=`$CC $SHOW 2> /dev/null` ; then
        AC_MSG_RESULT([Successful])
        for CCARG in $MPI_CC_RESULT ; do
          MPI_INCLUDES="$MPI_INCLUDES `echo $CCARG | grep '^-I'`"
        done
      else
        AC_MSG_RESULT([Failed])
      fi
    fi
  done
  if test "x$MPI_INCLUDES" != x; then
    MPI_INCLUDES=`echo $MPI_INCLUDES | sed -e 's/^ *//' -e 's/  */ /g'`
    AC_MSG_NOTICE([   Found MPI_INCLUDES $MPI_INCLUDES])
  fi
  if test "x$MPI_INCLUDES" != x ; then
    MPI_INCLUDE_PATH=`echo $MPI_INCLUDES | sed -e 's/^-I//'`
    MPI_INCLUDE_PATH=`echo $MPI_INCLUDE_PATH | sed -e 's/-I/:/g'`
    AC_MSG_NOTICE([   Found MPI_INCLUDE_PATH $MPI_INCLUDE_PATH])
  fi
fi
AC_SUBST([MPI_INCLUDES])
AC_SUBST([MPI_INCLUDE_PATH])
])

AC_DEFUN([SC_MPI_ENGAGE],
[
dnl determine compilers
m4_ifset([SC_CHECK_MPI_F77], [
AC_REQUIRE([AC_PROG_F77])
AC_PROG_F77_C_O
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_F77_WRAPPERS
])
m4_ifset([SC_CHECK_MPI_FC], [
AC_REQUIRE([AC_PROG_FC])
AC_PROG_FC_C_O
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
AC_FC_WRAPPERS
])
AC_REQUIRE([AC_PROG_CC])
AC_PROG_CC_C_O
AM_PROG_CC_C_O
m4_ifset([SC_CHECK_MPI_CXX], [
AC_REQUIRE([AC_PROG_CXX])
AC_PROG_CXX_C_O
])

dnl compile and link tests must be done after the AC_PROC_CC lines
if test "x$HAVE_PKG_MPI" = xyes ; then
dnl  m4_ifset([SC_CHECK_MPI_F77], [
dnl    AC_LANG_PUSH([Fortran 77])
dnl    SC_MPI_F77_COMPILE_AND_LINK(, [AC_MSG_ERROR([MPI F77 test failed])])
dnl    AC_LANG_POP([Fortran 77])
dnl  ])
dnl  m4_ifset([SC_CHECK_MPI_FC], [
dnl    AC_LANG_PUSH([Fortran])
dnl    SC_MPI_FC_COMPILE_AND_LINK(, [AC_MSG_ERROR([MPI FC test failed])])
dnl    AC_LANG_POP([Fortran])
dnl  ])
  SC_MPI_C_COMPILE_AND_LINK(, [AC_MSG_ERROR([MPI C test failed])])
  m4_ifset([SC_CHECK_MPI_CXX], [
    AC_LANG_PUSH([C++])
    SC_MPI_CXX_COMPILE_AND_LINK(, [AC_MSG_ERROR([MPI CXX test failed])])
    AC_LANG_POP([C++])
  ])
  if test "x$HAVE_PKG_MPIIO" = xyes ; then
    SC_MPIIO_C_COMPILE_AND_LINK(,
      [AC_MSG_ERROR([MPI I/O not found; you may try --disable-mpiio])])
  fi
  if test "x$HAVE_PKG_MPITHREAD" = xyes ; then
    SC_MPITHREAD_C_COMPILE_AND_LINK(,
      [AC_MSG_ERROR([MPI_Init_thread not found; you may try --disable-mpithread])])
  fi
  $1_ENABLE_MPIWINSHARED=yes
  SC_MPIWINSHARED_C_COMPILE_AND_LINK(,[$1_ENABLE_MPIWINSHARED=no])
  if test "x$$1_ENABLE_MPIWINSHARED" = xyes ; then
    AC_DEFINE([ENABLE_MPIWINSHARED], 1, [Define to 1 if we can use MPI_Win_allocate_shared])
  fi
  $1_ENABLE_MPICOMMSHARED=yes
  SC_MPICOMMSHARED_C_COMPILE_AND_LINK(,[$1_ENABLE_MPICOMMSHARED=no])
  if test "x$$1_ENABLE_MPICOMMSHARED" = xyes ; then
    AC_DEFINE([ENABLE_MPICOMMSHARED], 1, [Define to 1 if we can use MPI_COMM_TYPE_SHARED])
  fi
fi

dnl figure out the MPI include directories
SC_MPI_INCLUDES
])


dnl SC_CHECK_OPENMP(PREFIX)
dnl Check for OpenMP support and link a test program
dnl
dnl This macro tries to link to omp_get_thread_num both as is and with -lgomp.
dnl If neither of this works, we throw an error.
dnl Use the LIBS variable on the configure line to specify a different library.
dnl
dnl Using --enable-openmp without any argument defaults to -fopenmp.
dnl For different CFLAGS use --enable-openmp="-my-openmp-cflags" or similar.
dnl
AC_DEFUN([SC_CHECK_OPENMP], [

dnl This link test changes the LIBS variable in place for posterity
dnl ... which is bad if openmp is NOT used
dnl SAVE_LIBS="$LIBS"
dnl SC_CHECK_LIB([gomp], [omp_get_thread_num], [OPENMP], [$1])
dnl LIBS="$SAVE_LIBS"
AC_MSG_CHECKING([for OpenMP])

SC_ARG_ENABLE_PREFIX([openmp],
  [enable OpenMP:
     Using --enable-openmp without arguments does not specify any CFLAGS;
       to supply CFLAGS use --enable-openmp=<OPENMP_CFLAGS>.
     We check first for linking without any libraries and then with -lgomp;
       to avoid the latter, specify LIBS=<OPENMP_LIBS> on configure line
  ],
  [OPENMP], [$1])
if test "x$$1_ENABLE_OPENMP" != xno ; then
  dnl This is not portable: depends on compiler/architecture
  dnl $1_OPENMP_CFLAGS="-fopenmp"
  $1_OPENMP_CFLAGS=
  if test "x$$1_ENABLE_OPENMP" != xyes ; then
    $1_OPENMP_CFLAGS="$$1_ENABLE_OPENMP"
  fi
  PRE_OPENMP_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $$1_OPENMP_CFLAGS"
  SC_CHECK_LIB_NOCOND([gomp], [omp_get_thread_num], [OPENMP], [$1])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#include <omp.h>
]],[[
  omp_set_num_threads (2);
  #pragma omp parallel
  {
    int id = omp_get_thread_num ();
  }
]])],,
                 [AC_MSG_ERROR([Unable to link with OpenMP])])
dnl Keep the variables changed as done above
dnl CFLAGS="$PRE_OPENMP_CFLAGS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi
])


dnl sc_package.m4 - general custom macros
dnl
dnl This file is part of the SC Library.
dnl The SC library provides support for parallel scientific applications.
dnl
dnl Copyright (C) 2008,2009,2014 Carsten Burstedde, Lucas C. Wilcox.

dnl See SC_AS_SUBPACKAGE below for documentation of the mechanism

dnl Documentation for macro names: brackets indicate optional arguments

dnl SC_CHECK_INSTALL(PREFIX, REQUIRE_INCLUDE, REQUIRE_LDADD,
dnl                  REQUIRE_CONFIG, REQUIRE_ETC)
dnl The REQUIRE_* arguments can be either "true" or "false" (without quotes).
dnl This function throws an error if the variable PREFIX_DIR does not exist.
dnl The package must have been make install'd in that directory.
dnl Optionally require include, lib, config, and etc subdirectories.
dnl Set the shell variable PREFIX_INSTALL to "yes."
dnl
AC_DEFUN([SC_CHECK_INSTALL],
[
if test "x$$1_DIR" = xyes ; then
  AC_MSG_ERROR([Please provide an argument as in --with-PACKAGE=<directory>])
fi
if test ! -d "$$1_DIR" ; then
  AC_MSG_ERROR([Directory "$$1_DIR" does not exist])
fi
$1_INSTALL=yes
$1_INC="$$1_DIR/include"
$1_LIB="$$1_DIR/lib"
$1_CFG="$$1_DIR/share/aclocal"
$1_ETC="$$1_DIR/etc"
if $2 && test ! -d "$$1_INC" ; then
  AC_MSG_ERROR([Specified installation path $$1_INC not found])
fi
if $3 && test ! -d "$$1_LIB" ; then
  AC_MSG_ERROR([Specified installation path $$1_LIB not found])
fi
if $4 && test ! -d "$$1_CFG" ; then
  AC_MSG_ERROR([Specified installation path $$1_CFG not found])
fi
if $5 && test ! -d "$$1_ETC" ; then
  AC_MSG_ERROR([Specified installation path $$1_ETC not found])
fi
])

dnl SC_CHECK_PACKAGE(PREFIX, REQUIRE_INCLUDE, REQUIRE_LDADD,
dnl                  REQUIRE_CONFIG, REQUIRE_ETC)
dnl The REQUIRE_* arguments can be either "true" or "false" (without quotes).
dnl This function throws an error if the variable PREFIX_DIR does not exist.
dnl Looks for PREFIX_DIR/src to identify a source distribution.
dnl If not found, package must have been `make install`ed, in this case
dnl optionally require include, lib, config, and etc directories.
dnl Set the shell variable PREFIX_INSTALL to "yes" or "no."
dnl
AC_DEFUN([SC_CHECK_PACKAGE],
[
if test ! -d "$$1_DIR" ; then
  AC_MSG_ERROR([Directory "$$1_DIR" does not exist])
fi
if test -d "$$1_DIR/src" ; then
  $1_INSTALL=no
  $1_INC="$$1_DIR/src"
  $1_LIB="$$1_DIR/src"
  $1_CFG="$$1_DIR/config"
  $1_ETC=
  if $4 && test ! -d "$$1_CFG" ; then
    AC_MSG_ERROR([Specified source path $$1_CFG not found])
  fi
else
  SC_CHECK_INSTALL([$1], [$2], [$3], [$4], [$5])
fi
])

dnl---------------------- HOW SUBPACKAGES WORK ---------------------------
dnl
dnl A program PROG relies on libsc (or another package called ME, which
dnl could itself be using libsc).  In a build from source, me usually
dnl resides in PROG's subdirectory me.  This location can be overridden by
dnl the environment variable PROG_ME_SOURCE; this situation can arise when
dnl both PROG and me are subpackages to yet another software.  The
dnl path in PROG_ME_SOURCE must be relative to PROG's toplevel directory.
dnl All of this works without specifying a configure command line option.
dnl However, if me is already make install'd in the system and should
dnl be used from there, use --with-me=<path to me's install directory>.
dnl In this case, PROG expects subdirectories etc, include, lib, and
dnl share/aclocal, which are routinely created by me's make install.
dnl
dnl SC_ME_AS_SUBPACKAGE(PREFIX, prefix, ME, me)
dnl Call from a package that is using this package ME as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if me is make install'd.
dnl
AC_DEFUN([SC_ME_AS_SUBPACKAGE],
[
$1_$3_SUBDIR=
$1_$3_MK_USE=
$1_DISTCLEAN="$$1_DISTCLEAN $1_$3_SOURCE.log"

SC_ARG_WITH_PREFIX([$4], [path to installed package $4 (optional)], [$3], [$1])

if test "x$$1_WITH_$3" != xno ; then
  AC_MSG_NOTICE([Using make installed package $4])

  # Verify that we are using a me installation
  $1_DIST_DENY=yes
  $1_$3_DIR="$$1_WITH_$3"
  SC_CHECK_INSTALL([$1_$3], [true], [true], [true], [true])

  # Set variables for using the subpackage
  $1_$3_AMFLAGS="-I $$1_$3_CFG"
  $1_$3_MK_USE=yes
  $1_$3_MK_INCLUDE="include $$1_$3_ETC/Makefile.$4.mk"
  $1_$3_CPPFLAGS="\$($3_CPPFLAGS)"
  $1_$3_LDADD="$$1_$3_DIR/lib/lib$4.la"
else
  AC_MSG_NOTICE([Building with source of package $4])

  # Prepare for a build using me sources
  if test "x$$1_$3_SOURCE" = x ; then
    if test -f "$1_$3_SOURCE.log" ; then
      $1_$3_SOURCE=`cat $1_$3_SOURCE.log`
    else
      $1_$3_SOURCE="$4"
      $1_$3_SUBDIR="$4"
      AC_CONFIG_SUBDIRS([$4])
    fi
  else
    AC_CONFIG_COMMANDS([$1_$3_SOURCE.log],
                       [echo "$$1_$3_SOURCE" >$1_$3_SOURCE.log])
  fi
  $1_$3_AMFLAGS="-I \$(top_srcdir)/$$1_$3_SOURCE/config"
  $1_$3_MK_INCLUDE="include \${$2_sysconfdir}/Makefile.$4.mk"
  $1_$3_CPPFLAGS="-I$$1_$3_SOURCE/src -I\$(top_srcdir)/$$1_$3_SOURCE/src"
  $1_$3_LDADD="$$1_$3_SOURCE/src/lib$4.la"
fi

dnl Make sure we find the m4 macros provided by me
AC_SUBST([$1_$3_AMFLAGS])

dnl We call make in this subdirectory if not empty
AC_SUBST([$1_$3_SUBDIR])

dnl We will need these variables to compile and link with me
AM_CONDITIONAL([$1_$3_MK_USE], [test "x$$1_$3_MK_USE" != x])
AC_SUBST([$1_$3_MK_INCLUDE])
AC_SUBST([$1_$3_CPPFLAGS])
AC_SUBST([$1_$3_LDADD])
])


dnl SC_CHECK_PTHREAD(PREFIX)
dnl Check for POSIX thread support and link a test program
dnl
dnl This macro tries to link to pthread_create both as is and with -lpthread.
dnl If neither of this works, we throw an error.
dnl Use the LIBS variable on the configure line to specify a different library.
dnl
AC_DEFUN([SC_CHECK_PTHREAD], [

dnl This link test changes the LIBS variable in place for posterity
dnl ... which is bad if pthreads are NOT used
dnl SAVE_LIBS="$LIBS"
dnl SC_CHECK_LIB([pthread], [pthread_create], [LPTHREAD], [$1])
dnl LIBS="$SAVE_LIBS"
AC_MSG_CHECKING([for POSIX threads])

SC_ARG_ENABLE_PREFIX([pthread],
  [enable POSIX threads:
     Using --enable-pthread without arguments does not specify any CFLAGS;
       to supply CFLAGS use --enable-pthread=<PTHREAD_CFLAGS>.
     We check first for linking without any libraries and then with -lpthread;
       to avoid the latter, specify LIBS=<PTHREAD_LIBS> on configure line
  ],
  [PTHREAD], [$1])
if test "x$$1_ENABLE_PTHREAD" != xno ; then
  $1_PTHREAD_CFLAGS=
  if test "x$$1_ENABLE_PTHREAD" != xyes ; then
    $1_PTHREAD_CFLAGS="$$1_ENABLE_PTHREAD"
  fi
  PRE_PTHREAD_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $$1_PTHREAD_CFLAGS"
  SC_CHECK_LIB_NOCOND([pthread], [pthread_create], [LPTHREAD], [$1])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
#include <pthread.h>
void *start_routine (void *v)
{
  return NULL;
}
]],[[
  pthread_t thread;
  pthread_create (&thread, NULL, &start_routine, NULL);
  pthread_join (thread, NULL);
  pthread_exit (NULL);
]])],,
                 [AC_MSG_ERROR([Unable to link with POSIX threads])])
dnl Keep the variables changed as done above
dnl CFLAGS="$PRE_PTHREAD_CFLAGS"

  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi
])

m4_include([config/libtool.m4])
m4_include([config/ltoptions.m4])
m4_include([config/ltsugar.m4])
m4_include([config/ltversion.m4])
m4_include([config/lt~obsolete.m4])
m4_include([config/p4est_include.m4])
m4_include([config/p4est_metis.m4])
m4_include([config/p4est_petsc.m4])
