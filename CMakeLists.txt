cmake_minimum_required(VERSION 2.8)
project(Bembel CXX)

message(STATUS "

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

	
This is Bembel, the higher order C++ boundary element library. It was 
written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz, 
M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
source code is subject to the GNU General Public License version 3 and
provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
information.
")


#   Instructions for VisualStudio users (tested with VS 2019 Community)
#
#   Assuming you do not have some preixisting setup, you can compile 
#   tests and examples as follows: Uncomment the options
#
#   add_compile_options(-bigobj)
#   include_directories(eigen3)
#
#   in this file, and comment out the line  

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#   below. Then place the eigen3 directory obtained from downloading 
#   Eigen in Bembels root directory, where this CMake file sits. 
#   Afterwards, comment or remove the 
#
#   target_link_libraries(example_Quadrature.out Eigen3::Eigen)
#
#   lines in the CMake-Files of "\examples" and "\tests".


set (BEMBEL_PATH "${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 11)
include_directories(${BEMBEL_PATH})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" MARCHNATIVE)
if (MARCHNATIVE)
message(STATUS "Found -march=native, adding option to builds of type >Release< only")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" CACHE STRING "compile flags" FORCE)
endif()

CHECK_CXX_COMPILER_FLAG("-flto" FLTO)
if (FLTO)
message(STATUS "Found -flto, adding option to builds of type >Release< only")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto" CACHE STRING "compile flags" FORCE)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
message(STATUS "Found -fopenmp, adding option to builds of type >Release< only")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")
endif()

set(CMAKE_CXX_FLAGS_COVERAGE "-g --coverage -fprofile-abs-path" CACHE STRING
    "Flags used by the coverage test."
    FORCE)

add_subdirectory(tests)
add_subdirectory(examples)
enable_testing()
