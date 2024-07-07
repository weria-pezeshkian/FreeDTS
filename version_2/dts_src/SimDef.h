// This file contains a series of preprocessor directives defining various constants, flags,
// and configurations used throughout the software system.
// These constants are utilized to configure the behavior, settings, and options of the software,
// allowing for easy customization and adaptation to different requirements.
//
// Mathematical constants:
// - PI: Represents the mathematical constant pi.
// - S60: Represents the square root of 3 divided by 2.
// - SQ3: Represents the square root of 3.
//
// Constants related to system configuration:
// - FFType: Defines the forcefield type used in the system.
// - Inclusion_Type_Number: Specifies the number of inclusion types in the system.
// - SoftWareVersion: Indicates the version of the software.
// - Precision: Specifies the precision of calculations, such as the number of decimal places.
//
// Enabled/Disabled flags:
// - Enabled: Indicates the enabled state or option.
// - Disabled: Indicates the disabled state or option.
//
// Random number generator types:
// - UNIFROMTYPE1: Specifies the type of random number generator (not yet completed, avoid using).
// - UNIFROMTYPE0: Specifies an alternative type of random number generator.
//
// File extensions:
// - InExt: Specifies the file extension for input files.
// - TSIExt: Specifies the file extension for TSI files.
// - TopExt: Specifies the file extension for topology files.
// - TSExt: Specifies the file extension for TS files.
// - BTSExt: Specifies the file extension for BTS files.
// - RestartExt: Specifies the file extension for restart files.
//
// Names and identifiers:
// - S_Name: Represents the name of the software.
// - EXE_NAME: Represents the name of the executable file.
//
// Library configurations:
// - Libarmadillo: Specifies the status of the Armadillo library.
// - EigenLib: Specifies the status of the Eigen library.
// - NoLib: Specifies the absence of any external library.
//
// Testing and debugging modes:
// - TEST_MODE: Specifies the status of the test mode.
// - DEBUG_MODE: Specifies the status of the debug mode.
//
// Random number generator configuration:
// - RNGTYPE: Specifies the type of random number generator used.
// - BACKMAP: Specifies the status of the backmapping feature.
//
// These constants are essential for configuring the software system and can be adjusted as needed to customize its behavior.
#include <iostream>
#include <algorithm>
#include <iomanip> // Include this header for std::setprecision
#include "Nfunction.h"

// Mathematical constants
#define PI 3.14159265359  // Value of pi
#define S60 0.8660254037844  // Square root of 3 divided by 2
#define SQ3 1.73205080757  // Square root of 3


#define FreeDTS_BUFFERSIZE   4096// Adjust buffer size as needed

// Constants related to system configuration
#define FFType 2  // Forcefield type
#define Inclusion_Type_Number 10  // Number of inclusion types
#define SoftWareVersion 2.0  // Software version
#define Precision 8  // Precision of writing in files 

// Enabled/Disabled flags
#define Enabled 1  // Flag indicating enabled state
#define Disabled 2  // Flag indicating disabled state

// Random number generator types
#define UNIFROMTYPE1 2  // Random number generator type 1 (not yet completed, avoid using)
#define UNIFROMTYPE0 1  // Random number generator type 0

// File extensions
#define InExt "dts"  // Input file extension
#define TSIExt "tsi"  // TSI file extension
#define TopExt "top"  // Topology file extension
#define TSExt "q"  // TS file extension
#define BTSExt "bts"  // BTS file extension
#define RestartExt "res"  // Restart file extension
#define TimeSeriDataExt "-en.xvg"  // 

// Names and identifiers
#define S_Name "DTS"  // Software name
#define EXE_NAME "DTS"  // Executable name

// Library configurations
#define Libarmadillo Disabled  // Armadillo library status
#define EigenLib Disabled  // Eigen library status
#define NoLib Enabled  // No library status

// Testing and debugging modes
#define TEST_MODE Disabled  // Test mode status
#define DEBUG_MODE Disabled  // Debug mode status
#define DEVELOPMENT_MODE Disabled  // Test mode status

// Random number generator configuration
#define RNGTYPE UNIFROMTYPE0  // Random number generator type


//Precisions

#define TSI_Precisions  "18.10"
