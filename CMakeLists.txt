cmake_minimum_required(VERSION 3.5)

project(CryoSimulation)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/CryoSimulation")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Pedro (BWH)")
set(EXTENSION_DESCRIPTION "MonteCarlos simulation totest the impact of placement error in prostate cryoablation")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/CryoSimulation.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/CryoSimulation/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CryoMC)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
