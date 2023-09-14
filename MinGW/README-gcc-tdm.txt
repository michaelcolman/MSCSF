                    ________________________________________
                  _/_                                      _\_
               __/__/  TDM-GCC Compiler Suite for Windows  \__\__
              | « « |             GCC 5 Series             | » » |
               ¯¯\¯¯\         MinGW 32-bit Edition         /¯¯/¯¯
                  ¯\¯                                      ¯/¯
                    ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯


This edition of TDM-GCC is an unofficial replacement for the official GCC
binaries distributed by the MinGW project; please note the following caveats:
 * TDM-GCC is not formally affiliated with or endorsed by the MinGW project.
 * No level of support for TDM-GCC is in any way guaranteed, although a best
     effort is made to fix bugs as they are found or forward them to GCC
     Bugzilla.


BUGS:
If you encounter a problem while using a TDM-GCC build that isn't present in a
previous MinGW or TDM release, you are encouraged to submit a helpful bug
report. Please see <http://tdm-gcc.tdragon.net/bugs> for further instructions.


/=-=- ------------ -=-=\
|=|   INSTALLATION   |=|
\=-=- ------------ -=-=/

-=- TDM/MINGW INSTALLER -=-

Using the TDM/MinGW installer is highly recommended; it can automatically
install TDM-GCC (or the official MinGW GCC) as well as all supplementary base
system packages. The installer uses a standard wizard interface with reasonable
defaults.

-=- MANUAL INSTALLATION -=-

Do not install TDM-GCC packages on top of a previous GCC installation of any
kind.

You will need to download and unpack a set of archives. A minimal base set of
archives is required; there are also some additional components that are
optional, adding support for additional programming languages or GCC features.

TDM-GCC provides a ZIP-compressed version and a TAR.LZMA-compressed version of
each archive. Use whichever is easiest.

REQUIRED BASE:
 * gcc-core (gcc-5.1.0-tdm-1-core)
 * binutils (binutils-2.24-1-mingw32-bin)
 * mingwrt (mingwrt-3.20-2-mingw32-dev, mingwrt-3.20-2-mingw32-dll)
 * w32api (w32api-3.17-2-mingw32-dev)

OPTIONAL:
 * gcc-c++ (gcc-5.1.0-tdm-1-c++) - C++ support
 * gcc-ada (gcc-5.1.0-tdm-1-ada) - Ada support
 * gcc-fortran (gcc-5.1.0-tdm-1-fortran) - Fortran support
 * gcc-objc (gcc-5.1.0-tdm-1-objc) - Objective-C/C++ support
 * gcc-openmp (gcc-5.1.0-tdm-1-openmp) - OpenMP support
 * mingw32-make (make-3.82.90-2-mingw32-cvs-20120902-bin,
     libintl-0.17-1-mingw32-dll-8, libiconv-1.13.1-1-mingw32-dll-2) - GNU make
     for *-mingw32 GCC
 * gdb (gdb32-7.9-tdm-1) - GNU source-level debugger, Python enabled, for
     mingw32
You'll need GDB particularly if you want to use an IDE with debugging support.

Unpack all the archives to an empty directory. You may choose any path, though
it is recommended that you avoid a path with any spaces in the folder names.
Finally, consider adding the bin subdirectory to your Windows PATH environment
variable.


/=-=- ----------- -=-=\
|=|   USAGE NOTES   |=|
\=-=- ----------- -=-=/

-=- GDB32: A TDM BUILD OF THE GNU DEBUGGER -=-

The "gdb32" package is a slightly-modified build of GDB for 32-bit Windows that
includes Python support, enables libstdc++ Python pretty printing by default,
and uses wrapper executables to allow 64-bit and 32-bit cooperation in the
"bin" directory. For more details, see the separate README file
(README-gdb32-tdm.txt).

-=- POSIX THREADS, C++11 STD::THREAD AND OTHER SYNCHRONIZATION FEATURES -=-

As of the 4.8.1 TDM-GCC release, a significant change has been introduced which
allows you to use C++11 features such as std::thread that rely on the POSIX
threading library. TDM-GCC now includes a new pthreads compatibility layer
called "winpthreads" instead of the old "pthreads-w32".

"winpthreads" is one of the libraries distributed by the MinGW-w64 project, and
it allows GCC to be compiled with full pthreads compatibility, which is
necessary to enable std::thread and other threading related functions in the
C++ runtime.

Because of TDM-GCC's continuing goal of minimizing extra DLLs, winpthreads has
been compiled statically. It will be statically linked with every program you
compile, which will increase your baseline executable size.

The same mechanism used in libgcc and libstdc++ to allow EXEs and DLLs to share
state for handling exceptions has also been patched into winpthreads to allow
your EXEs and DLLs to share C++11 thread handles via the underlying pthread
handles.

Because every program you compile will now rely on winpthreads, you should make
sure to read and comply with its MIT-style license, included in the file
"COPYING.winpthreads.txt".

-=- INLINE MEMBER FUNCTIONS AND DLLS -=-

[[[ IMPORTANT NOTE:
[[[ You will probably need to use "-fno-keep-inline-dllexport" when building
[[[ large DLLs with lots of inline member functions, such as the "wxWidgets"
[[[ library's monolithic DLL.

As of the 4.5 series, GCC conforms more closely to the behavior of MSVC on
Windows platforms by always emitting inline functions that are class members
when creating a DLL. This behavior doesn't necessarily conform to the
expectations of libraries that are used to the old behavior, however, and in
some cases can cause the linker to run out of memory when creating the DLL. If
necessary, use the "-fno-keep-inline-dllexport" flag to avoid emitting these
functions.

-=- LTO (LINK-TIME OPTIMIZATION) -=-

Every TDM-GCC release since 4.5.1 includes support for GCC's Link-Time
Optimizer. As long as GCC's own drivers (gcc, g++, etc.) are used at both
compile-time and link-time, and the "-flto" option is specified at both compile-
time and link-time, link-time optimization will be applied. See
<http://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html> for further details.

-=- DW2 VS. SJLJ UNWINDING -=-

GCC currently supports two methods of stack frame unwinding: Dwarf-2 (DW2) or
SJLJ (setjmp/longjmp). Until recently, only SJLJ has been available for the
Windows platform. This affects you, the end user, primarily in programs that
throw and catch exceptions. Programs which utilize the DW2 unwind method
generally execute more quickly than programs which utilize the SJLJ method,
because the DW2 method incurs no runtime overhead until an exception is thrown.
However, the DW2 method does incur a size penalty on code that must handle
exceptions, and more importantly the DW2 method cannot yet unwind (pass
exceptions) through "foreign" stack frames: stack frames compiled by another
non-DW2-enabled compiler, such as OS DLLs in a Windows callback.

This means that you should in general choose the SJLJ version of the TDM-GCC
builds unless you know you need faster exception-aware programs and can be
certain you will never throw an exception through a foreign stack area.

As distributed, the SJLJ and DW2 packages of TDM-GCC can coexist peacefully
extracted to the same directory (i.e. any files in common are for all intents
and purposes identical), because the driver executables (the ones in the "bin"
directory) are suffixed with "-dw2" for the DW2 build, and the libraries and
other executables hide in another "-dw2" directory in "lib(exec)/gcc/mingw32".
This allows you to use the same single addition to your PATH, and use DW2
exceptions only when you need them by calling "gcc-dw2", etc. If you truly want
DW2 exceptions as the default when calling "gcc" (from Makefiles or configury
systems, for example), you can rename or copy the suffixed executables to their
original names.

-=- EXCEPTIONS AND DLLS -=-

[[[ IMPORTANT NOTE:
[[[ TDM-GCC uses a statically-linked libstdc++ by default! To use the libstdc++
[[[ DLL, specify "-shared-libstdc++" on the command line.

With the advent of the GCC 4.5 release series, the mingw32 port finally supports
fully the same method every other platform uses to allow exceptions to propagate
out of shared libraries (DLLs): gcc library DLLs. For any GCC language that
supports exceptions (and DLLs), this method requires the runtime presence of two
additional DLLs: (1) libgcc_s*.dll, which contains common core data, and (2) a
language-specific DLL.

However, TDM-GCC also continues to integrate a versioned shared memory region
for the *static* (non-DLL) runtime libraries, which will still allow you to
throw exceptions between any DLLs or executables that are built with TDM-GCC.
This method incurs a very small execution overhead as compared to the shared
library method, but has the very important benefit of not requiring you to
redistribute extra DLLs with your program.

By default, TDM-GCC will continue to create executables and DLLs that use the
static libraries and do not require you to redistribute further DLLs. If you
would like to use the shared libraries, you should add "-shared-libgcc" to the
command line to use a shared version of libgcc, and additionally ensure that the
shared version of your language-specific runtime library is being used. For C++,
add "-shared-libstdc++".

You cannot use a shared version of libgcc with a static version of a language-
specific runtime. The reverse -- static libgcc with shared language-specific
runtime -- should work fine.

IMPORTANT NOTE:
There has been an update to the license exception clause that permits you to
distribute programs that make use of the GCC runtime libraries without requiring
you to license your programs under the GPLv3. As always, please be familiar with
the terms of GCC's GPLv3 license and exception clauses, and do not redistribute
any portion of GCC, including its runtime DLLs, in any way except as granted by
the license. If you are unclear about which permissions are granted by the
license, please consult a lawyer and/or the Free Software Foundation
(<http://www.fsf.org/>).

A copy of the GPLv3 may be found in the file
COPYING-gcc-tdm.txt, and a copy of the runtime library exception clause may be
found in COPYING.RUNTIME-gcc-tdm.txt. In general, the runtime library exception
clause probably applies to any file found in the "lib" directory or its
subdirectories, and any DLL found in the "bin" directory -- but you should
consult the sources, available for download from the TDM-GCC project site on
SourceForge, if you are unsure.

-=- OPENMP AND WINPTHREADS -=-

TDM-GCC has been built to allow the use of GCC's "-fopenmp" option for
generating parallel code as specified by the OpenMP API. (See
<http://gcc.gnu.org/onlinedocs/libgomp/> for details.) If you want to use
OpenMP in your programs, be sure to install the "openmp" optional package.

The OpenMP support in the TDM-GCC builds has received very little testing; if
you find build or packaging problems, please send a bug report (see BUGS above).

LibGOMP, GCC's implementation of OpenMP, currently only supports the use of the
POSIX Threads (pthreads) api for implementing its threading model. Because the
MinGW project itself doesn't distribute a pthreads implementation, the
"winpthreads" library, available as part of the MinGW-w64 project libraries, is
included in this distribution. The "winpthreads" library is distributed under
the terms of an MIT-style license; see "COPYING.winpthreads.txt" for details.

In order to correctly compile code that utilizes OpenMP/libGOMP, you need to add
the "-fopenmp" option at compile time AND link time. By default, this will link
the static version of winpthreads to your program, and you should not need to
distribute any additional DLLs with your program. If you plan to distribute a
program that relies on OpenMP and winpthreads, be sure to understand and comply
with the terms of winpthreads' license (see COPYING.winpthreads.txt).

"libpthread.a" is included in the "lib" subdirectory of the openmp package along
with three other pthreads library files:
 - "libpthread_s.a" and "libwinpthread.dll.a" link to a DLL version of
     winpthreads - libwinpthread-1.dll.
 - "libwinpthread.a" is the same as "libpthread.a".

-=- WARNINGS AND ERRORS	-=-

Recent GCC releases make significant strides in optimization capabilities, error
detection, and standards compliance. For you, the end user, this often means
that code which compiled and ran without problems on previous GCC releases will
exhibit some warnings and maybe even a few errors.

These meaningful warnings and errors are a very good thing, as they help the
programmer to write safer and more correct code. Unfortunately, there's also a
chance you might encounter incorrect warnings or errors, ICE's (internal
compiler errors, where the compiler makes a mistake and has to bail out), or
even miscompilations (where your code is incorrectly compiled and produces the
wrong result).

If you encounter an ICE while using a TDM-GCC build, feel free to file a bug
report (see BUGS above). With any other unexpected problem, you are urged to
work from the assumption that it stems from user error, and ensure that your
code is correct and standards-compliant.


/=-=- --------------------- -=-=\
|=|   BUGS AND KNOWN ISSUES   |=|
\=-=- --------------------- -=-=/

As these builds are provided on the same basis as the source releases, and the
mingw32 target in GCC tends to receive somewhat less-than-average attention,
some bugs are expected. If you encounter a bug that you are certain is in the
GCC sources (such as an ICE), or that is due to an issue in the building or
packaging process, you are encouraged to report it. Please visit the TDM-GCC
Bugs Page at <http://tdm-gcc.tdragon.net/bugs> for bug reporting instructions.


/=-=- ----------------------- -=-=\
|=|   LOCAL FIXES AND CHANGES   |=|
\=-=- ----------------------- -=-=/

 - [make-rel-pref.patch] Includes a patch to keep GCC from erroneously using the
     CWD as the installation directory.
 - [lfs.patch] Includes a patch to re-enable large file support for C++ fstreams
     (LFS detection fails because there is no definition for struct stat64 in
     mingw-runtime).
 - [libgomp.patch] Includes a patch to allow libgomp to interoperate correctly
     with user-generated pthreads. See
     <https://sourceforge.net/p/tdm-gcc/bugs/76/>.
 - [libgcceh.patch] Includes a patch which reintegrates the code from
     libgcc_eh.a into libgcc.a and the libgcc DLL. As long as the shared memory
     region is used to handle exceptions in the static runtimes, this library is
     unnecessary, and it causes multiple definition errors for the symbols in it
     because it hasn't been added to binutils' exception libraries yet.
 - [defstatic.patch] Includes a patch to make libgcc and libstdc++ link
     statically by default. The "-static-libgcc" and "-static-libstdc++"
     commands thereby do not have any effect; use "-shared-libgcc" and
     "-shared-libstdc++" to link your program to the corresponding DLLs.
 - [ada-lfs.patch] Includes a patch to allow Ada to build for older versions of
     the MSVCRT without a stat64 equivalent.
 - [ssp-mingw32.patch] Includes a patch to allow libssp to build for older
     versions of the MSVCRT without CryptAcquireContext.
 - [relocate.patch] Includes a patch to make all search paths for headers,
     libraries and helper executables relative to the installation directory of
     the driver executables -- in other words, TDM-GCC is fully relocatable and
     does not search any absolute system paths.
 - [eh_shmem.patch] Includes a patch to propagate C++ exceptions out of DLLs
     without the need for shared runtime libraries.
 - [threads.patch] Includes a patch to fix GCC & winpthreads interop for x86
     mingw32 and allow use of a fully static pthreads library.
 - [more-gnattools.patch] Includes a patch to enable extra tools in the Ada
     toolchain for mingw*.
 - [windows-lrealpath.patch] Includes a patch to allow forward slashes in
     libiberty as path separators on Windows, allowing build system integrations
     like map files to work better.
 - [mutex-leak.patch] Includes a patch to prevent pthreads mutexes from leaking
     in libstdc++.
 - [lto-binary.patch] Fix LTO reader to open files in binary mode.
 - [xmmintrin.patch] Add C++ include guards to xmmintrin.h.
 - [crtbegin.patch] Remove static modifier from __EH_FRAME_BEGIN__ to allow
     visibility across compilation units.
 - [gnat-implibs.patch] Create import libraries for the DLL versions of libgnat
     and libgnarl.
 - [mcrtdll.patch] Allow specifying newer MSVCRT versions with -mcrtdll=.
 - [dw2-reg-frame.patch] Prevent DW2 frame register/unregister from getting
     mistakenly stripped.
 - Configured with "--enable-fully-dynamic-string", which fixes a bug when
     passing empty std::string objects between DLLs and EXEs.


/=-=- ----------- -=-=\
|=|   SOURCE CODE   |=|
\=-=- ----------- -=-=/

The source code for the TDM-GCC binary releases is available from the TDM-GCC
download page on SourceForge: <http://sourceforge.net/projects/tdm-gcc/files/>.
(The most up-to-date link to the download site will always be available at
<http://tdm-gcc.tdragon.net/>.)

The source is distributed in the form of the original ("vanilla") separate
source packages as downloaded, plus an additional "TDM Sources" package. The TDM
Sources package includes unified diffs of any changes made to the vanilla
sources, as well as the set of scripts used to build the binary releases.


/=-=- ------------------ -=-=\
|=|   COMPONENT LICENSES   |=|
\=-=- ------------------ -=-=/

The GCC proper packages in TDM-GCC contain binary distributions constituting a
work based on GCC, ISL, MPC, libiconv, GMP, MPFR and winpthreads.

* GCC itself is licensed under the GPLv3; for further details, see
  "COPYING3-gcc-tdm.txt". GCC's runtime libraries are licensed with an
  additional exception clause; see "COPYING.RUNTIME-gcc-tdm.txt".

* MPC, libiconv, GMP, and MPFR are each licensed under the LGPLv3, a somewhat
  more permissive version of the GPL; see "COPYING3.LIB-gcc-tdm.txt".

* ISL and winpthreads use an "MIT-style" license, reproduced in
  "COPYING.ISL.txt" and "COPYING.winpthreads.txt".
