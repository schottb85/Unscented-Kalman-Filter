# kalman
Unscented Kalman filter (UKF) C++ implementation based on Udacity Nanodegree Programm "Sensor Fusion Engineer"

# Project content
Kalman filter are used in this project to estimate **positioning** and **velocity** and **yaw angle/rate** using lidar or/and radar measurement data.
- An *unscented Kalman filter* is implemented in ukf.cpp and ukf.h
- tools.cpp/.h controls the ukf filter and its visualization
- highway.cpp/.h define an example scenario for testing


# Installation
`git clone https://github.com/schottb85/kalman.git`

`cd kalman`

`mkdir build`

`cd build`

`cmake ..`

`make`

# Variants
##Visualization options
There are different flags to be selected for different visualization in highway.h
-	`bool visualize_lidar = true;`
-	`bool visualize_radar = true;`
-	`bool visualize_pcd = false;`

##Lidar and/or Radar UKF measurements
- in ukf.cpp two parameters `bool use_laser_` and `bool use_lidar_` can be used activate/deactivate the use of lidar and radar measurements, respectively

##Prediction of path
the path prediction for the different cars can be configured using e.g.
- `double projectedTime = 2;`
- `int projectedSteps = 6;`
defining the prediction/projection time in seconds and the number of discrete ukf projection steps
