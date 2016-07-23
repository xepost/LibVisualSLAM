LibVisualSLAM
=============
LibVisualSLAM is an algorithm library accompany with [CoSLAM](https://github.com/danping/CoSLAM).

System requirements
-------------

* Linux Unbuntu or Linux Mint (64 bit). 

Dependent packages:
-------------
* OpenCV
* libatlas-dev
* libswscale
* libavformat
* libavutil
* libavcodec


Download & Installation
-------------
1. OpenCV

    Download and install OpenCV, notice that there openCV must be compiled with
    ffmpeg and v4l and  -D WITH_OPENCL=OFF 

2. Libatlas
	
    sudo apt-get install libatlas-dev

3. LibVisualSLAM

    git clone https://github.com/danping/LibVisualSLAM
    cd LibVisualSLAM
    cmake .
    make
    sudo make install
    
    
Liciense
-------------
LibVisualSLAM is released under [GPL v2.0](http://www.gnu.org/licenses/gpl-2.0.html).
