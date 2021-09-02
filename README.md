# medwatch

A simple set of utilities to interact with the medwatch.  Tested on a mac, but should work on windows and linux with no modification.

Requires numpy and bleak.  After installing python 3x, install these using:

pip install numpy bleak
 
Currently there are three programs (found in the scripts folder).  These are all heavily based on the bleak example code.:

**blelist**:  Lists nearby BLE devices

**findmedwatch**: finds any BLE devices nearby whose device name starts with "MED-WATCH"

**readwatch**: Continuously reads the output of a running medwatch, dumps the output to the terminal.  Currently it will attach to the first medwatch it finds.
