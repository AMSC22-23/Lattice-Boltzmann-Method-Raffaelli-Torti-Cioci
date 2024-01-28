# Lattice Boltzmann Method in C++ and CUDA

Welcome to our cutting-edge implementation of the Lattice Boltzmann Method (LBM), a powerful computational fluid dynamics technique. In this project, we seamlessly blend clarity in CPU code with optimized GPU performance, offering an efficient and versatile solution. Whether you're a developer exploring our object-oriented CPU code for its clarity or delving into the GPU-accelerated calculations for maximum performance, this documentation should provide a comprehensive guide to navigating and understanding our LBM implementation.

## Building
After cloning the repo, the first step is to verify you have all the software needed to run the project. You'll need CMake, a C++ compiler (GCC suggested), Python, ffmpeg and, optionally, CUDA Toolkit.

Create python virtual environment and install requirements

```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

Compile the project

```bash
./compile.sh
```

## Running

**You'll find ALL outputs in the `outputs` folder.**

Test run the project

```bash
./run.sh lid-driven-cavity 100
```

Optionally, try running a more complex example on the gpu

```bash
./run.sh ball 100 --gpu
```

To generate an airfoil, run the following command

```bash
python scripts/airfoil_generator.py
```

You can also generate an obstacle of your own, by loading from a png of only black and white pixels

```bash
python scripts/png_reader.py
```

Circular obstacles can be generated with the following command

```bash
python scripts/ball.py
```

For the usage of `plotting2D.py` and `dragLift.py`, please check out `run.sh`

# What, How, Why and Performance Analysis
In the following sections, we will discuss the main features of our implementation and the rationale behind our design choices. We will also provide a brief overview of the performance of our code, comparing the CPU and GPU implementations.

## *What*

### Lattice Boltzmann Method
We have implemented an advanced Computational Fluid Dynamics (CFD) method known as the Lattice Boltzmann Method (LBM), grounded in the principles of particle dynamics and specific space-time discretization techniques. LBM, an accurate mesoscale model, leverages the Kinetic Theory and the Lattice Gas Model to depict fluid behavior by simulating particle movement on a regular lattice. Despite its inherent computational complexity, this model finds widespread application in large-scale simulations.

### Our Implementation
Our model is characterized by a two-dimensional portrayal of fluid dynamics within a lattice. Depending on the nature of the problem at hand, the fluid may either be confined within the lattice or allowed to pass through it. The foundational structure of our code not only facilitates seamless two-dimensional simulations but also sets the stage for future implementations of three-dimensional simulations. This adaptability is made possible through the integration of the `Structure` class and the `NDimensionalMatrix` template class, allowing for modifications in specific methods while preserving the fundamental structures and overall code organization.

### Technologies used
- Ubuntu 22.04
- CMake 3.22.1
- GCC 11.4.0
- CUDA 12.4
- C++ 17 (with OpenMP)
- Python 3.10.12
- ffmpeg 4.4.2

### How different problems are represented
We have developed a system that loads problems from plain text files. Each problem type is denoted by a number, and additional simulation parameters are also specified in the file. Some problem types support arbitrarily shaped obstacles in the lattice, and the shapes of those are also outlined in the text file.

### Separable compilation
The CMakeLists file is designed to cater to a broad range of systems. Computers with the CUDA toolkit installed will compile the executable with CUDA support, while those without the CUDA toolkit will still compile but lack support for GPU execution.

### Implementation Choices
We have chosen to use the **D2Q9 velocity set** along with their respective weights as the discretization of the velocity field to describe the possible movements of particles within the lattice.

For managing interactions among fluid particles within our discretized environment, we've implemented the **Two-Relaxation-Time (TRT)** model. This model employs two distinct relaxation times, tailored for different components of the particle distribution function. The use of TRT effectively addresses limitations present in alternative models, offering enhanced flexibility, consistency, and accuracy in describing a wide array of fluid dynamic phenomena.

We have selected the **Zou-He boundary conditions** to govern the interaction between fluid particles and solid surfaces in our implementation. These conditions facilitate an efficient and accurate representation of the no-slip condition, assuming that the fluid velocity at the surface contact equals the velocity of the surface itself. The Zou-He conditions play a crucial role in calculating the unknown particle distributions in cells adjacent to the domain boundary. By enforcing the bounce-back rule for the non-equilibrium part of particles, these conditions ensure the correct reflection of particles, contributing to the overall fidelity of our simulations.

The conditions governing the behavior of the fluid in contact with an obstacle employ the **non-interpolated bounce-back** method.

### Lift and Drag
Lift and Drag forces are essential for understanding the behavior of a moving fluid and are particularly relevant in evaluating specific sections of a project, such as vehicles, airplanes, and marine structures.

The **Momentum Exchange Method**, at the core of the implementation, manages microscopic interactions in fluids. After the collision phase of particles, momentum exchange occurs between adjacent cells, allowing the LBM to accurately model complex flows, velocity gradients, and non-stationary flows.

For a parabolic setInlets configuration, we anticipate a monotonic behavior in the lift and drag graphs, barring the presence of significant turbulence. The parabolic setInlets should lead to a predictable and smoothly varying flow, resulting in **gradual changes in lift and drag forces**. Any deviations from monotonic behavior may indicate the influence of turbulent phenomena within the simulated fluid.

## *HOW*

### Code structure
The CPU code is structured as an Array of Structures, following an object-oriented paradigm, facilitating a clear code structure and a seamless implementation. On the other hand, the GPU code utilizes a Structure of Arrays data structure, optimizing for maximum performance through memory coalescence. This approach capitalizes on the GPU memory bus width, reducing the overall number of data requests when multiple threads require neighboring data.

### Input Files
The requisite information within these text files includes:
- The problem type (1 for the lid-driven cavity, 2 for parabolic inlet / poiseuille flow)
- The width and height of the lattice
- Parameters such as the Reynolds number, total number of steps, and maximum forced velocity
- For problem type 2, a list of coordinates representing cells designated as obstacles

Furthermore, our program supports a problem type 3, facilitating the *inaccurate* simulation of a ball oscillating up and down.

### Output Files
The program generates output files within the `outputs` folder, with the output displayed every *total steps / frames* steps, where frames is a user-specified parameter from the command line.

- **velocity_out.txt**: This file contains the velocity field at each time step.
- **lift_drag_out.txt**: It includes the drag and lift forces at each time step.

The program also produces visualization files:
- **movie[...].mp4**: This video, generated by the `plotting2D.py` script, illustrates the simulation's velocity modulus at each time step.
- **lift_drag.png**: A plot of drag and lift forces over time, created by the `dragLift.py` script.

### Parallelization
Parallelization is achieved through both OpenMP and CUDA, utilizing multiple threads on the CPU and the parallel architecture of the GPU, respectively. The OpenMP implementation employs the `#pragma omp parallel for` directive, automatically distributing the workload among available threads. The CUDA implementation utilizes a 2D grid of blocks, where each block contains a grid of 24x24 threads. The number of blocks is determined by the lattice size. Implementing OpenMP involves minimal code changes, allowing seamless execution with or without OpenMP support.

### Atomic Operations for Lift and Drag
The use of `#pragma omp atomic` directive and `atomicAdd` CUDA function is pivotal during force updates to guarantee the consistency and correctness of results in parallel environments. By enforcing atomic operations, this directive prevents write concurrency issues, preserving the integrity of the results.

## *Why*

### Design Choices
We opted for the **regular bounce-back** method for handling obstacles instead of the interpolated version. This choice allows us to achieve a sufficiently accurate simulation with a significantly simpler code structure and higher performance.

For the velocity set in two dimensions, we settled on **D2Q9**. This choice strikes a balance between accuracy and computational efficiency, with the nine velocity directions proving effective in capturing a wide range of fluid dynamic phenomena in two dimensions, including vortical flows and other complex behaviors.

In terms of **conventions for x and y coordinates** representing lattice cell positions, we adopted an approach aligned with matrix indexing. Specifically, an increase in the x coordinate signifies movement to the right, while an increase in the y coordinate corresponds to downward movement.

When organizing data in a memory data structure, we chose **row-major** order. In C++, sequential access to a matrix in row-major order is often more efficient. This decision aligns with how memory is loaded and managed in the processor cache during program execution.

### Object-Oriented Programming
The prevalence of object-oriented programming led us to adopt this approach over other contemporary methodologies, such as functional programming. This decision is grounded in the understanding that the primary performance gains come from the GPU. The primary objective of the CPU code is clarity and conciseness rather than aggressively optimizing for every last bit of available performance.

## Performance Analysis

### Credits
This project takes inspiration from [jviquerat's lbm python code](https://github.com/jviquerat/lbm)
