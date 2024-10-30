# How to install
There are several important steps you need to follow during the installation of the program. 
Here is a step by step guide with some troubleshooting steps to help you during the installation.

## Install Visual Studio
The project is written in Visual Studio and packaged with a Visual Studio solution file (.sln). It is 
thus highly recommended that you specifically use Visual Studio and not any other editor. Visual Studio
also comes preinstalled with Microsofts MSVC compiler, so there is no need to install a compiler yourself or
cmake. 

You may download Visual Studio with C++ via the following link or search for it on google by yourself:
https://visualstudio.microsoft.com/vs/features/cplusplus/

## Install gnuplot
The project makes use of the open source graphing library gnuplot. You will specifically need the version 
gnuplot Release 5.4.10, newer versions will not work! Here you can find a link to the gnuplot website:
http://www.gnuplot.info/
You may, of course, search for the gnuplot webiste by yourself via google. There you will find this 
website. You need press the button indicated by the red arrow in the image below:

![gnuplotDownload](https://github.com/user-attachments/assets/952cdef6-552f-40cd-a270-ace1086ff95f)

After clicking the link you will arrive at a sourceforge page, which provides the download links for 
the different gnuplot versions. You need to download the file gp5410-win64-mingw.exe as indicated in 
the image.

![sourceforgeDownload](https://github.com/user-attachments/assets/660612fe-6379-469d-b152-1a7d32c90b4d)

The download will start shortly afterwards. After clicking on the executable you may get a windwos defender 
warning as gnuplot is not a known trusted executable. It is however a well known and trusted program by 
the community. You may, of course, choose to inspect the source code and build the program yourself if you 
feel like it. Otherwise, click the "More info" option and press "Run anyway". Accept the license agreement
and follow the steps as indicated. 

You need to add gnuplot to your system path. You may do this manually or check the following box during the 
installation (this is recommended).

![pathvariable](https://github.com/user-attachments/assets/c7a467a2-55da-412a-8b4d-50a7067bf23f)

If you happened to miss this during the installtion, you may just run the installation executable file again
and click the box this time.

## Matplot++
The project relies on the Matplot++ library written by Alan de Freitas for plotting. There is no installation of 
Matplot++ required, as it already comes packaged and set up with the project. You may visit the Matplot++ GitHub
page via the following link: https://github.com/alandefreitas/matplotplusplus.

## Running the project
After cloning the GitHub project, you can open the solution file using Visual Studio. Since the Matplot++ library
only works in Release Mode, you should make sure to change the Visual Studio configuration to Release Mode | x64 before
running the program. Afterwards, you should just need to press CTRL + F5 and the program should run automatically. 
A console window will open itself and display data from the program. Shortly afterwards, the gnuplot figure
windows open automatically and display the data. 

# Troubleshooting
## gnuplot erros
If you get gnuplot errors in the console such as "gnuplot not found", you need to make sure that you have correctly
installed gnuplot and added gnuplot to your system path. It is best to make extra sure and redo the installation
of gnuplot and specifically check the add to system path box as instructed during the installation.

    line 0: warning: Reading from '-' inside a multiplot not supported; use a datablock instead

When this error occurs, you may or may not get an actual gnuplot figure. This errors means however, that you 
installed the wrong version of gnuplot i.e. 6.01 or higher. Make sure you install gnuplot version 5.4.10.

If the program runs without any errors, but there are no gnuplot figures created, this may mean, that you are in 
the incorrect configuration. Check if the Matplot++ code and functions (e.g. plot(x,y), figure(), etc.) are greyed out 
via a preprocessor block. This means, you are in the incorrect configuration and have to switch to Release mode.

## Linker Errors

    Error	LNK2038	mismatch detected for 'RuntimeLibrary': value 'MD_DynamicRelease' doesn't match value 'MDd_DynamicDebug' in main.obj	

This error occurs when the program is not set to release mode, and you are trying to run a Matplot++ command. Usually this does 
not occur, as all Matplot++ commands are excluded in Debug mode via preprocessor statements. If it occurs nontheless, you may 
need to check if the program is actually in Release Mode. 

If other linker or even compiler errors occur, make sure that the libraries are included correctly. To do this manually
go to: Project > Properties then look for the Configuration drop down menu, and select "All configurations". Afterwards, check if 
all of the following items are in the specified directories:

C/C++ > General > Additional Include Directories : `$(SolutionDir)\Dependencies\Matplot++\include` | `$(ProjectDir)\Header Files`
Linker > General > Additional Library Directories : `$(SolutionDir)\Dependencies\Matplot++\lib\Matplot++` | `$(SolutionDir)\Dependencies\Matplot++\lib`
Linker > General > Input > Additional Dependencies : `matplot.lib` | `nodesoup.lib`

    Error	LNK1104	cannot open file 'C:\CPProjects\CompPhys\x64\Release\UE1.exe'	UE1	C:\CPProjects\CompPhys\UE1\LINK	1		

This may occur when you didnt close the last console and output properly before rerunning the program. To fix this
simply make sure, that all console windows and gnuplot windows are closed. Then try rerunning the program.

## Compiler Errors

    Error	C2923	'std::pair': 'string_view' is not a valid template type argument for parameter '_Ty1'	

If an error like this occurs, make sure that your C++ standard is set to C++17 or higher for ALL configurations.
To do this, go to: 
Project > Properties and select "All Configurations" in the Configuration drop down menu then go to:
Configuration Properties > General > C++ Language Standard and select "ISO C++17 Standard (/std:c++17)

# Additional Information
When using the program, you may change the line 39: #define EXERCICE 0 to any of the following values: 0, 11, 12, 13, 2 or 3. This
selects which section of the program is run and which is not. Specifically, 0 means that the entire program is run, 11 means that
only the part for exercice 1a) is run and the rest corresponds to 1b), 1c), 2a) and 2b) or 3a) and 3b).

If you cannot get the gnuplot installation to work as intended, you may also run the program in Debug mode. This way, all
parts of the program related to plotting, are automatically excluded via preprocessor statements. You can then just
view the results as a console ouput. You may also check out the Graphs folder to find images of the graphs.

