#ASL_Sensor_Suite

This program is used to do various testing and calibration on seismic sensors.
Running this program requires SEED files and and corresponding RESP files. It is currently capable of doing self-noise tests using 3 sets of sensor data (SEED and RESP) and relative gain tests using 2 or 3 sets of data. A test for calculating step calibration response is currently in development but not complete. Additional tests are planned but not yet implemented.

If a sensor test requires more data than has currently been loaded, attempting to generate a result for that test will fail; the error message "INSUFFICIENT DATA LOADED" will display instead.

Currently, RESP files much define a response over a single epoch only. RESP files with multiple epochs cannot be parsed correctly; only one epoch will be loaded.

The program currently attempts to trim data so that all 3 displayed time series are truncated to their overlapping time series. If two SEED files do not have an intersecting time range, the second file will not be loaded. In that case, the data will have to be cleared out (press the "Clear Data" button at the bottom of the applet) before the other file can be loaded. 