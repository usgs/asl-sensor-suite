# ASL_Sensor_Suite

### Purpose

This program is used to analyze various aspects of seismic sensor data in order to determine information about their configuration, such as gain and orientation. It is meant to be a simple program that can be used to generate data on a wide range of sensor tests, a one-stop-shop for sensor diagnostics designed to replace several disparate programs across multiple languages. The program is designed around an interface meant to be simple and intuitive.

### Requirements
##### Software
This program is designed to be used with Java 1.8, but should also be compatible with Java 1.7
This program also requires Gradle in order to be run from source. For instructions on installing Gradle, see https://gradle.org/install
NOTE: because Gradle requires access to the maven repositories to download dependencies, there may be build issues when running the program on DOI-networked computers. The instructions for using DOI certificates for maven authentication can be found under the Java subheader at https://github.com/usgs/best-practices/blob/master/WorkingWithinSSLIntercept.md. You will also need the file "DOIRootCA.crt" from the linked page near the top detailing how to configure git.
For those using Mac computers, the last step for trusting the certificate in Java may be slightly different. If installing the certificate using the instructions above fails, try using this command instead (NOTE: requires administrative access due to the use of the sudo command) https://blog.alwold.com/2011/06/30/how-to-trust-a-certificate-in-java-on-mac-os-x/

##### Hardware
Because SEED files must be decompressed and stored in memory, the footprint of this application has the potential to be rather large. Running this program on OSX v. 10.11.15, the total memory usage was 76.6 MB compressed but 1.82 GB uncompressed. As a result, this program's performance may be dependent on the memory management systems used by the OS of the system it runs on.

### Compilation

While releases are regularly updated with major feature changes included in the releases tab on this project's Github and approximately-daily snapshots included among the projects files, it may be desired to build the project direct from source. There are a few ways of doing this, explained below.

##### Command Line
The program can be compiled by using the commands `gradle compileJava` which will compile the source code, or `gradle build` which will also run the unit tests.
Running the program can be done by either opening the jar through a filebrowser or running either `gradle run`, which launches the jar file, or `java -jar build/libs/SensorTestSuite$version_number$.jar` after the program has been built, with $version_number$ replaced with the current version, such as 0.9.0. The gradle build script also allows the built jar file to be placed in the root directory; if `gradle compileJava` was previously run, then `gradle copyJar` will move it there. Note that `gradle build` includes this step by default.

##### Eclipse
For those who wish to compile and run this program with Eclipse, run the command `gradle eclipse` and then, inside eclipse, go to File>"Open projects from file system..." and direct Eclipse to the root folder of the test suite. Now the code will be available as an Eclipse project. For more information on using Eclipse, consult the Eclipse documentation.

### File Selection
The program will default to looking for SEED files in the "data" subdirectory in the same folder the jar is, if it exists. It is possible to choose a different folder, which will become the new default folder when loading in additional files. It will not, however, persist after the program is closed.
SEED files must have an intersecting time range, which will be displayed when multiple files are loaded in. The bars below the input plots can be used to select a narrower range of data to zoom in on. Loading in a SEED file that does not have any common time range with the other data will produce an error; to reset the loaded SEED files, they can either be unloaded individually with a corresponding remove button for that data, or all data can be cleared out with the 'clear all' button.

The program also comes with a number of response files (in RESP format) embedded in the jar, selectable from a drop-down menu, that correspond to several common sensor configurations. It is also possible to load in other response files from the 'load custom response' option in a manner similar to SEED files described above, the default directory being the 'responses' folder in the same folder as the jar file. Note that clearing a SEED file also clears out the corresponding response.

### Output
Plots of the input files used in a sensor test can be output in PNG, as can the output of a given sensor test, using the corresponding buttons in each panel. Both can be compiled into a single-page PDF of all currently displayed data using the button at the bottom of the program.

### Usage

For more information on the specifics of certain tests, consult the javadoc. 

#### Self-noise

Self-noise requires three components and an appropriate response file for each. The test computes the cross-power (power-spectral density) of each pair of files, and uses that data to extract estimations of the inherent noise of each sensor. Plots of the seismic NLNM and NHNM are also included. Units of frequency (Hz) or period (seconds, default) can be selected using the checkmark in the bottom-left of the panel.

The input files do not need to be in any particular order. They all must have responses specified. For three-component self-noise, they should all be pointing in the same direction (i.e., all facing north).

There is also a nine-component self-noise test that takes in horizontal north, east, and vertical sensor data for each of the three components, finds the best angle to rotate the horizontal components to maximize coherence, and then performs the same test on the 3 sensors in each direction. 

#### Relative Gain

Relative gain computes the mean of the PSD of each of two sensors, and estimates the gain from the mean of the ratio of the values over a selected range of the output.

The input files, again, do not need to be in any order (the panel allows for choosing which sensor to be used as reference by way of the selection menus below the chart), but they must both have responses specified.

The gain is initially calculated using the octave around the peak frequency, but a custom range can be specified using the sliders, same as the 

#### Step Calibration

Step calibration takes in a step input signal and the response to that signal from a sensor. It attempts to find response parameters that best fit the application of that response to the input, or rather, the response that when deconvolved with the input's signal produces a function closest to the step.

The input files have a specific order: the step input signal must be placed first, though it does not use a response. The second input, then, is the output from the sensor of interest, and a response should be chosen to help guide the initial guess of the solver.

#### Randomized calibration

This function solves for poles to attempt to fit the response curve calculated from deconvolving the given calibration input from the sensor output. Low-frequency (the two lowest poles) and high-frequency (all other poles) are fitted to minimize the difference between the estimated response, based on the response specified for the sensor. The inputs follow the same structure as step calculation, though what response parameters are solved for is dependent on whether a high or low frequency calculation is chosen. Both the magnitude and argument (angle of the response curve along the real axis) of the response curve are displayed in plots, and saving the plot to an image will include both such plots.

When using embedded response files, it is strongly recommended to use an appropriate response file with "nocoil" in the name, as these remove the calibration coil's response from the file and thus generate more accurate results of calculations.

Note that plots have been scaled in order to produce more representative fits of response curves. For high-frequency calibrations, the curves are all set to be equal to zero at 1 Hz; for low-frequency calibrations, this point occurs at 0.2 Hz.

This function is still work-in-progress but has been tested with good results on data from a KS54000 sensor. Other sensors may not produce as good results (see known issues, below).

Older high-frequency cals may produce lots of noise on the high-frequency end confounding the solver, especially depending on how the calibration was produced.
This program includes a second checker tab which does not run the solver for 
response parameters, but can be used to determine whether or not the calculated response from the sensor output is good enough to be used for the solver. 
Noisy calibrations or ones whose output otherwise varies significantly from the given nominal response may take a long time to solve or produce errors that lead to the solver being unable to converge on any solution. 
IT IS STRONGLY ENCOURAGED TO RUN THE NO-SOLVER RANDOM CAL PANEL ON DATA THAT PRODUCES A CONVERGENCE ERROR DURING SOLVING IN ORDER TO DIAGNOSE POTENTIAL ISSUES WITH THE CALIBRATION ITSELF.

#### Azimuth

Azimuth takes in 3 inputs. The first two are orthogonal sensors assumed to be respectively facing near north and east. The third is a reference sensor assumed to point north, though the offset angle field can be used to specify a clockwise offset from north. The code will try to find a clockwise rotation angle that maximizes the coherence estimation between the rotated unknown-angle sensor data and the reference angle. This angle is added to the offset to produce the (clockwise) azimuth estimation. A value of the coherence estimations per-frequency for the found angle is also given as a separate plot. 

#### Orthogonality

Orthogonality takes in four inputs, two each from sensors known or assumed to be orthogonal, and finds the true (interior) angle between the second two sensors if the first two sensors are truly orthogonal.

The input files have a specific order: the first and third inputs are for north-facing sensors, and the second and fourth are for east-facing sensors. As noted above, the first two sensors are assumed to be 90 degrees apart for the purpose of the test; the second two sensors' orientation is what is solved for.


#### Response

This plots 1-3 different response attenuation and phase curves (Bode plots) for given response files. The image generated from this plot will include both plots, though the program can only display one at a time (selectable with the drop-down menu in the bottom-left of the panel). Units of frequency (Hz, default) or period can be selected by the selection box on the bottom-right, much like with the self-noise plot.

This program also allows for extracting a response file embedded in the program, in order to be edited by hand. This way a field engineer with access to the program will be able to define a custom response file using the text editor of their choice. The response selected will be copied into the "responses" subdirectory with a date-stamped version of the same name as the embedded response file. (Copying the responses is done so that the nominal responses embedded in the program cannot be edited unwittingly. These new response files can be loaded in using the "load custom response" option with the response loader pane).

### Further Work / Known Issues

Currently the application does its best to show the complete range among all data, there are some issues in doing so. If there are three SEED files loaded and the first two SEED files have more data than the third, then when switching to a test using only two inputs, the entire range of the first two sensors should be visible. However, if there are loaded inputs not included in a test and a new file is loaded in one of the input slots, it must still have a time range in common with the unused inputs. While not ideal behavior, it prevents additional bugs from handling non-matching time ranges if a test using the non-active data is selected again.

Most sensors have a specific response related to the calibration signal produced by their calibration coils. These responses are necessary for producing accurate plots of the calculated response from a calibration, but do not yet exist as part of the program. As a result, trying to solve for poles from a random calibration is likely to produce incorrect results for sensors that require calibration response correction. The KS54000 is an example of a sensor that does not require such correction, and results using the calibration test with one is much closer to expectation.

## DISCLAIMER

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.
