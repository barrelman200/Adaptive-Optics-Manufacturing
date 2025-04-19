# Adaptive-Optics-Manufacturing
MATLAB program that implements the Gerchberg-Saxton algorithm to generate multiple foci, using a phase mask written onto an SLM.
Example shown below:
![Image](https://github.com/user-attachments/assets/91af7e3b-4157-4fdc-9e2b-1bbe4134efc0)
This version contains aberration correction for a tilt (linear phase ramp), but future versions will include more general aberration correction from Zernike polynomials.

The program allows the user to input a wide range of parameters, including the number of foci, real spacing, lens NA, etc. However, this code is specialised for foci generation, and not general-purpose hologram generation (there is no option to convert an image yet).
