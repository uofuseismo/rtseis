
# Already in cache, be silent
if (IPP_INCLUDE_DIR AND IPP_IPPS_LIBRARY AND IPP_VM_LIBRARY AND IPP_CORE_LIBRARY)
   set (IPP_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
#   set(IPPS "libipps.a")
#   set(VM "libippvm.a")
#   set(CORE "libippcore.a")
#else()
   set(IPPS "ipps")
   set(VM "ippvm")
   set(CORE "ippcore")
#endif()

find_path(IPP_INCLUDE_DIR
          NAMES ipps.h
          HINTS $ENV{IPP_INC_DIR}
                $ENV{IPP_ROOT}/include
                /opt/intel/ipp/include
                /opt/intel/oneapi/ipp/latest/include)
find_library(IPP_IPPS_LIBRARY
             NAMES ${IPPS}
             PATHS $ENV{IPP_LIB_DIR}
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   /opt/intel/ipp/lib/intel64
                   /opt/intel/ipp/lib
                   /opt/intel/oneapi/ipp/latest/lib/intel64)
find_library(IPP_VM_LIBRARY
             NAMES ${VM}
             PATHS $ENV{IPP_LIB_DIR}
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   /opt/intel/ipp/lib/intel64
                   /opt/intel/ipp/lib
                   /opt/intel/oneapi/ipp/latest/lib/intel64)
find_library(IPP_CORE_LIBRARY
             NAMES ${CORE}
             PATHS $ENV{IPP_LIB_DIR}
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   /opt/intel/ipp/lib/intel64
                   /opt/intel/ipp/lib
                   /opt/intel/oneapi/ipp/latest/lib/intel64)

set(IPP_LIBRARY ${IPP_IPPS_LIBRARY} ${IPP_VM_LIBRARY} ${IPP_CORE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FindIPP DEFAULT_MSG IPP_LIBRARY IPP_INCLUDE_DIR IPP_IPPS_LIBRARY IPP_VM_LIBRARY IPP_CORE_LIBRARY)
mark_as_advanced(IPP_INCLUDE_DIR IPP_LIBRARY IPP_IPPS_LIBRARY IPP_VM_LIBRARY IPP_CORE_LIBRARY)
