cmake_minimum_required(VERSION 2.8)
project(Bembel CXX)

MACRO(MACRO_ENSURE_OUT_OF_SOURCE_BUILD MSG)
     STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}"
     "${CMAKE_BINARY_DIR}" insource)
     GET_FILENAME_COMPONENT(PARENTDIR ${CMAKE_SOURCE_DIR} PATH)
     STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}"
     "${PARENTDIR}" insourcesubdir)
    IF(insource OR insourcesubdir)
        MESSAGE(FATAL_ERROR "${MSG}")
    ENDIF(insource OR insourcesubdir)
ENDMACRO(MACRO_ENSURE_OUT_OF_SOURCE_BUILD)


MACRO_ENSURE_OUT_OF_SOURCE_BUILD(
    "${CMAKE_PROJECT_NAME} requires an out of source build."
)

message(STATUS "
          ____  ______ __  __ ____  ______ _         _____           
         |  _ \\|  ____|  \\/  |  _ \\|  ____| |       / ____|_     _   
         | |_) | |__  | \\  / | |_) | |__  | |      | |   _| |_ _| |_ 
         |  _ <|  __| | |\\/| |  _ <|  __| | |      | |  |_   _|_   _|
         | |_) | |____| |  | | |_) | |____| |____  | |____|_|   |_|  
         |____/|______|_|  |_|____/|______|______|  \\_____|          
                                                                                                                      

                         ,;:cclloooooooodddddddoolc;, 
                      .cc::::::;;;;;;;,,,,,,;;;;::c:,.             
                  ...'lkc',;;::::::::::::::::;;,'..,cd;            
              ';:::::::dOdllccc::::::::;;;:::ccloddoll'            
            ,ll,..     'oc.',,,,;;;;;;;;;;;;,,''''lk:              
          .:o,      ...'okdooollllllllllllllllllookd.              
         .co.    .;cc::cxx'                      .oo.              
         ,d;    .ol.    ;d,                       co.              
         :o.    :o.    .ll.  .:lo;.               'll'             
         co.    co.  .cll;..ckOxdc'....            .,lc.           
         ;d,    ,d;.cOx:od:dkollodO00Ok:.            .:o;.         
         .ll.   .:xkOOkkxloxoxkl,'lk0XX0:              'oc.        
          'oc.  .:kl..'dkkkdOk;'coxO0XXKl.              .ol.       
           .ll'.;d:...',l00O0c.,kXXXXXKk,                'o:       
            .;ooxo.'d00kooOKO;  'lxkko;.   .              :d'      
              .ckc.:kkxO0OkKKc .;okOx;.;ldxxdl.           'd:      
               'd:.,ll::coOKKx';O0o:;ckKKK00K0o.          .oc      
               .oc;kKXK0xl:lkOll0x. ,k0o;,,;:c;.          .oc       
               .lock0xoodxxocllcxd'.lOc..,oxOOOd:.        'd;      
                ,d:,coxxdddddl:':o;.lo..l0XXXXXXKd.       co.      
                .lo'c00occcloddl:cc;::.:0KklcxKXXOc;'    .o:       
                 .oll0d. ':,',:lol::;'.o0o.  :OOkdOKl.  .co.       
                  'odxk:.oK0c..:cc;,...lO:   .,';xOl.  .:d,        
                   'oddkk0KXk,;dxdol:'.'dd'  ,xlok:.;,.;d;         
                    'oo:lxOkc..;ldxkkkxlldd;.'xddd;lkc;d:.         
                     .ol.... .dKKkdxOKXK0OOOkdkkkxdkkxx:.          
                      .ll.   .ckk:;dOKOdc,,;:ldOKKK0Okc.           
                       .ll.    ...'cokkdccdxdlloxxOOOl.            
                        .ol.         ...'cl;':dl;codl.             
                         ;ko;'....         ..,:;;lxx,              
                         cxlclooolllllllllllllllccdx'              
                         ,lc:;,,''''''''''''',,;;:l:.              
                          ..,;:::::::::::::::::;,'.                


Hi. I am the CMake-File of Bembel. Note that this software requires Eigen 3
to be installed. Bembel itself was written by Jürgen Dölz, Helmut Harbrecht, 
Michael Multerer (born Peters) and Felix Wolf. This software was written as
a basis for boundary element codes of higher order. It relies heavily on 
template libraries, and thus, compiler optimisation. Please, only use this 
if compiled with -O2 or -O3. See <http://www.bembel.eu/> for further infor-
mation.

This program is part of free software, software which you can redistribute 
and/or modify under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, or (at your 
option) any later version. The software is distributed in the hope that it 
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
Public License for more details. You should have received a copy of the GNU 
General Public License along with this software. If not, you will find it at
<http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------
")

## SETTING COMPILE FLAGS

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" MARCHNATIVE)

## -D_BEMBEL_OLD_ is required to enable support of the legacy C vector3 struct in the spline module
SET ( CMAKE_CXX_FLAGS "-std=c++11 -O3 -D_BEMBEL_OLD_" CACHE STRING "compile flags" FORCE)

if (MARCHNATIVE)
SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native" CACHE STRING "compile flags" FORCE)
endif()

link_directories(${CMAKE_BINARY_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
endif()

message( "Compiling everything with ${CMAKE_CXX_FLAGS}" )

##

set( PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
include_directories()
include_directories("${PROJECT_SOURCE_DIR}/include/includeC")
include_directories("${PROJECT_SOURCE_DIR}/include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include/")
include_directories("${PROJECT_SOURCE_DIR}")

add_library(bem STATIC
src/bemlibC/pdeproblem
src/bemlibC/pdeproblem_LaplaceSingle
src/bemlibC/pdeproblem_HelmholtzSingle
src/bemlibC/pdeproblem_MaxwellSingle
src/bemlibC/discretization
src/bemlibC/discretization_NB
src/bemlibC/projector
src/bemlibC/glue_maxwell
src/bemlibC/trafos
src/bemlibC/gauss_legendre
src/bemlibC/gauss_square
src/bemlibC/topology2DK
src/bemlibC/H_print_geometry
src/bemlibC/meshdata
src/bemlibC/BEMRHS_LaplaceDirichlet_can
src/bemlibC/BEMRHS_HelmholtzDirichlet_can
src/bemlibC/BEMRHS_MaxwellDirichlet_can
src/bemlibC/error_can
src/bemlibC/pot_can
src/bemlibC/fuev
src/bemlibC/gram
src/bemlibC/gram_can
src/bemlibC/hmatrixsettings
src/bemlibC/element_tree
src/bemlibC/cluster_tree
src/bemlibC/int_core
src/bemlibC/Hl2_FtimesV
src/bemlibC/Hl2_getHdiag
src/bemlibC/multipoleH2_can
src/bemlibC/tree_leafs_H2
src/bemlibC/H2transforms
src/bemlibC/H2l2_HtimesV
src/bemlibC/H2l2_RtimesV
src/bemlibC/constants
src/bemlibC/glue_laplace
src/bemlibC/Visualization
src/bemlibC/basis
src/bemlibC/patch
src/bemlibC/deformation
src/bemlibC/shape
)

target_compile_definitions(bem PRIVATE _DK _BEMBEL_OLD_ _BEMBEL_TEST_WITH_ASSERT_)

add_library(bemtest STATIC 
  src/test/test_phi
  src/test/test_phiphi
  src/test/test_projector
  src/test/test_glue
  src/test/test_HMatrixWrapper
  src/test/test_testdivconf
  src/test/test_divdiv
  src/test/test_geom_sphere
  src/test/test_geom_screen
)

add_library(spltest STATIC 
  src/test/spline/test_bernstein
  src/test/spline/test_deBoor
  src/test/spline/test_import
  src/test/spline/test_knots
  src/test/spline/test_localise
  src/test/spline/test_patch
)

target_link_libraries(bemtest bem)
target_link_libraries(spltest bem)

add_executable(bem_testsuite.out src/test/testsuite)
target_link_libraries(bem_testsuite.out bem bemtest)

add_executable(spline_testsuite.out src/test/spline/testsuite)
target_link_libraries(spline_testsuite.out bem spltest)

add_executable(LaplaceSingle.out src/examples/example_LaplaceSingle)
target_link_libraries(LaplaceSingle.out bem)

add_executable(HelmholtzSingle.out src/examples/example_HelmholtzSingle)
target_link_libraries(HelmholtzSingle.out bem)

add_executable(MaxwellSingle.out src/examples/example_MaxwellSingle)
target_link_libraries(MaxwellSingle.out bem)

add_executable(Visualization.out src/examples/example_Visualization)
target_link_libraries(Visualization.out bem)

add_executable(VisualizeDatfile.out src/bemlibC/vtkMaker)
target_link_libraries(VisualizeDatfile.out bem)

