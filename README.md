Welcome to the GitHub repository for the hands-on project of the Advanced Methods for Scientific Computing (AMSC) course. This project focuses on the implementation of the Lattice Boltzmann Method (LBM) in C++, providing a powerful tool for simulating complex fluid dynamics and computational physics.

Authors:
- Martina Raffaelli
- Andrea Torti
- Marco Cioci

## Requirements
The project builds using cmake and make and compiles using g++, then requiring ffmpeg for python to render the video.
The project requires python to be installed and a virtual environment to be created. To do so, type the following commands in the terminal:
```
$ python3 -m venv env
```
then activate the virtual environment:
```
$ source env/bin/activate
```
and install the required packages:
```
$ pip install -r requirements.txt
```

## Usage
To run the project, simply type the following command in the terminal:
```
$ ./run.sh
```
