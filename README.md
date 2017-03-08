#ASL_Sensor_Suite

###Purpose

This program is used to analyze various aspects of seismic sensor data in order to determine information about their configuration, such as gain and orientation. It is meant to be a simple program that can be used to generate data on a wide range of sensor tests, a one-stop-shop for sensor diagnostics designed to replace several disparate programs across multiple languages. The program is designed around an interface meant to be simple and intuitive.

###Requirements
#####Software
This program is designed to be used with Java 1.8, but should also be compatible with Java 1.7
Apache Ant is used to generate builds. Jar compilation was successful using Ant 1.9.7. 
Future releases plan to use Gradle as primary build tool; because Gradle is backward-compatible with Ant configurations, Gradle can be installed instead of Ant in order to create jar builds.

#####Hardware
Because SEED files must be decompressed and stored in memory, the footprint of this application has the potential to be rather large. Running this program on OSX v. 10.11.15, the total memory usage was 76.6 MB compressed but 1.82 GB uncompressed. As a result, this program's performance may be dependent on the memory management systems used by the OS of the system it runs on.

###Compilation
A jar file can be built from command line inside the code subfolder by using the command "ant jar". The intermediate files created during compilation can be removed with "ant clean". Running the program can be done by either opening the jar through a filebrowser or running the command "java -jar SensorTestSuite.jar".

###File Selection
The program will default to looking for SEED files in the "data" subdirectory in the same folder the jar is, if it exists. It is possible to choose a different folder, which will become the new default folder when loading in additional files. It will not, however, persist after the program is closed.
SEED files must have an intersecting time range, which will be displayed when multiple files are loaded in. The bars below the input plots can be used to select a narrower range of data to zoom in on. Loading in a SEED file that does not have any common time range with the other data will produce an error; to reset the loaded SEED files, they can either be unloaded individually with a corresponding remove button for that data, or all data can be cleared out with the 'clear all' button.

The program also comes with a number of response files (in RESP format) embedded in the jar, selectable from a drop-down menu, that correspond to several common sensor configurations. It is also possible to load in other response files from the 'load custom response' option in a manner similar to SEED files described above, the default directory being the 'responses' folder in the same folder as the jar file. Note that clearing a SEED file also clears out the corresponding response.

###Output
Plots of the input files used in a sensor test can be output in PNG, as can the output of a given sensor test, using the corresponding buttons in each panel. Both can be compiled into a single-page PDF of all currently displayed data using the button at the bottom of the program.

###Usage

For more information on the specifics of certain tests, consult the javadoc. 

####Self-noise

Self-noise requires three components and an appropriate response file for each. The test computes the cross-power (power-spectral density) of each pair of files, and use that data to extract estimations of the inherent noise of each sensor. Plots of the seismic NLNM and NHNM are also included.

The input files do not need to be in any particular order. They all must have responses specified.

####Relative Gain

Relative gain computes the mean of the PSD of each of two sensors, and estimates the gain from the mean of the ratio of the values over a selected range of the output.

The input files, again, do not need to be in any order (the panel allows for choosing which sensor to be used as reference), but they must both have responses specified.

####Step Calibration

Step calibration takes in a step input signal and the response to that signal from a sensor. It attempts to find response parameters that best fit the application of that response to the input, or rather, the response that when deconvolved with the input's signal produces a function closest to the step.

The input files have a specific order: the step input signal must be placed first, though it does not use a response. The second input, then, is the output from the sensor of interest, and a response should be chosen to help guide the initial guess of the solver.

####Orthogonality

Orthogonality takes in four inputs, two each from sensors known or assumed to be orthogonal, and finds the true (interior) angle between the second two sensors if the first two sensors are truly orthogonal.

The input files have a specific order: the first and third inputs are for north-facing sensors, and the second and fourth are for east-facing sensors. As noted above, the first two sensors are assumed to be 90 degrees apart for the purpose of the test; the second two sensors' orientation is what is solved for.

### Further Work / Known Issues

Currently the application does its best to show the complete range among all data, there are some issues in doing so. If there are three SEED files loaded and the first two SEED files have more data than the third, then when switching to a test using only two inputs, the entire range of the first two sensors should be visible. However, if there are loaded inputs not included in a test and a new file is loaded in one of the input slots, it must still have a time range in common with the unused inputs. While not ideal behavior, it prevents additional bugs from handling non-matching time ranges if a test using the non-active data is selected again.

There are plans to include additional tests in this project. Azimuth and the 9-input self noise calculation are intended to be included.