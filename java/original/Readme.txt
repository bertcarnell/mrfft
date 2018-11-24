cct: Contains utility functions (like Complex) called by my code.  

FFTTestData: Test data Mike ran.  The numbers are the length of the array.  The file with pts is the input
array and the file without pts is the output.   Due to miscommunication Mike used my input array twice, so really what you have here is a complex array for inputs with real and complex component matching.  It should serve just as well for verification as the all real arrays I had intended.

signal_processing: has the code I've been developing off the FORTRAN code.  There is a driver and then the algorithm in MRFFT.  You only need the java files, the class files are created when you build your app

.classpath/.project:  You will recreate when you set up Eclipse

fftsingleton.f:  the FORTRAN implementation

*.xls:  Various informal tests I was playing with.  Your mileage my vary.
