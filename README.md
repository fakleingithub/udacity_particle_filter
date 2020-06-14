# **Particle Filter** 

---

**Kidnapped Vehicle Project**

The goals / steps of this project are the following:
* Localization of the robot car with the usage of a particle filter.
* Particle filter localizes the vehicle to within desired accuracy.
* Particle filter runs within the specified time of 100 seconds.

[//]: # (Image References)

[gif1]: ./media/particle_filter.gif "Particle Filter"

## Project Introduction
>Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

>In this project you will implement a 2 dimensional particle filter in C++. Your particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step your filter will also get observation and control data.

## Project Result

Here you can see the final simulation result of the particle filter:

![alt text][gif1]

# code structure
The directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data
|   |   
|   |   map_data.txt
|   
|   
|___src
    |   helper_functions.h
    |   main.cpp
    |   map.h
    |   particle_filter.cpp
    |   particle_filter.h
```

#### code compilation

The main program can be built and run by doing the following from the project top directory:

```sh
1. ./clean.sh
2. ./build.sh
3. ./run.sh
```


