#ifndef ROOT_RConfigOptions
#define ROOT_RConfigOptions

#define R__CONFIGUREOPTIONS   "CMAKE_CXX_STANDARD_LIBRARIES=kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib CMAKE_C_STANDARD_LIBRARIES=kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib GLEW_INCLUDE_DIR=C:/ROOT-CI/src/builtins/glew/inc GLEW_INCLUDE_DIRS=C:/ROOT-CI/src/builtins/glew/inc GLEW_LIBRARIES=GLEW::GLEW GLEW_LIBRARY=$<TARGET_FILE:GLEW> GSL_CBLAS_LIBRARY=C:/libs/x64/GSL/2.7/lib/Release/gslcblas.lib GSL_CBLAS_LIBRARY_DEBUG=C:/libs/x64/GSL/2.7/lib/Debug/gslcblas.lib GSL_INCLUDE_DIR=C:/libs/x64/GSL/2.7/include GSL_LIBRARY=C:/libs/x64/GSL/2.7/lib/Release/gsl.lib GSL_LIBRARY_DEBUG=C:/libs/x64/GSL/2.7/lib/Debug/gsl.lib GSL_VERSION=2.7 LZ4_INCLUDE_DIR=C:/ROOT-CI/src/builtins/lz4 LZ4_INCLUDE_DIRS=C:/ROOT-CI/src/builtins/lz4 LZ4_LIBRARIES=LZ4::LZ4 LZ4_LIBRARY=$<TARGET_FILE:lz4> LZ4_VERSION=1.9.3 LZ4_VERSION_STRING=1.9.3 ODBC_INCLUDE_DIR=C:/Program Files (x86)/Windows Kits/10/Include/10.0.22621.0/um ODBC_LIBRARY=odbc32.lib OPENGL_gl_LIBRARY=opengl32 OPENGL_glu_LIBRARY=glu32 PCRE_INCLUDE_DIR=C:/ROOT-CI/build/builtins/pcre/PCRE-prefix/src/PCRE-build PCRE_LIBRARIES=C:/ROOT-CI/build/builtins/pcre/PCRE-prefix/src/PCRE-build/Release/pcre.lib PCRE_PCRE_LIBRARY=C:/ROOT-CI/build/builtins/pcre/PCRE-prefix/src/PCRE-build/Release/pcre.lib PCRE_VERSION=8.43 ZLIB_INCLUDE_DIR=C:/ROOT-CI/src/builtins/zlib ZLIB_INCLUDE_DIRS=C:/ROOT-CI/src/builtins/zlib ZLIB_LIBRARIES=ZLIB::ZLIB ZLIB_VERSION=1.2.8 ZLIB_VERSION_STRING=1.2.8 ZSTD_INCLUDE_DIR=C:/ROOT-CI/src/builtins/zstd ZSTD_INCLUDE_DIRS=C:/ROOT-CI/src/builtins/zstd ZSTD_LIBRARIES=ZSTD::ZSTD ZSTD_LIBRARY=$<TARGET_FILE:ZSTD> ZSTD_VERSION=1.4.8 ZSTD_VERSION_STRING=1.4.8 xxHash_INCLUDE_DIR=C:/ROOT-CI/src/builtins/xxhash xxHash_INCLUDE_DIRS=C:/ROOT-CI/src/builtins/xxhash xxHash_LIBRARIES=xxHash::xxHash xxHash_LIBRARY=$<TARGET_FILE:xxhash> xxHash_VERSION=0.8.0 xxHash_VERSION_STRING=0.8.0 "
#define R__CONFIGUREFEATURES  "cxx17  asimage builtin_afterimage builtin_cfitsio builtin_clang builtin_cling builtin_cppzmq builtin_freetype builtin_ftgl builtin_gl2ps builtin_glew builtin_gtest builtin_llvm builtin_lz4 builtin_lzma builtin_nlohmannjson builtin_openui5 builtin_pcre builtin_tbb builtin_unuran builtin_xxhash builtin_zeromq builtin_zlib builtin_zstd clad dataframe fitsio fortran gdml http imt mathmore odbc opengl pyroot roofit webgui root7 rpath shared tmva tmva-cpu tmva-pymva spectrum unuran"

#endif
