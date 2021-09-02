# medwatch

A simple set of utilities to interact with the medwatch.  Should work on mac, windows, and linux.

Requires numpy and bleak.  After installing python 3x, install these using:

pip install numpy bleak
 
Currently there are three programs (found in the scripts folder):

**blelist**:  Lists nearby BLE devices

**findmedwatch**: finds any BLE devices nearby whose device name starts with "MED-WATCH"

**readwatch**: Continuously reads the output of a running medwatch, dumps the output to the terminal.  Currently it will attach to the first medwatch it finds.
