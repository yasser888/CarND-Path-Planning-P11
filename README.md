# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

https://youtu.be/lP9r2uGjh_E

[![Watch the video](https://img.youtube.com/vi/lP9r2uGjh_E/0.jpg)](https://youtu.be/lP9r2uGjh_E)

First we have defined constana values to use:

 1- NUM_OF_LANES 3 : Number of lane on each side of the road which is 3 here.

 2- LANE_WIDTH 4.0 : Lane of each lane in width is 4m

 3- MPH_TO_MPS 0.44704 : Main car's speed provided by the simulator is in MPH. constant is used to convert it into Meters per seconds as position of the car and other distance are given in meters

 4- TIME_BETWEEN_WAYPOINTS 0.02 : Simulator travel a distance between 2 consecutive waypoints in 0.02 seconds

 5- MAX_VELOCITY_IN_MPH 47 : Maximum allowed Speed is 50 MPH in simulator , but to have some buffer I have set it to 47 MPH

 6- MAX_ACCELERATION 7 : Max allowed acceleration is 10 mps^2 kept it 7 mps^2 to have some buffer.

 7- FRONT_CLEAR_DISTANCE_FOR_LANE_CHANGE 30 : It defines clear distance (without any vehicle) in front of our car, This distance of 30m is used at the time of making lane change decision.

 8- REAR_CLEAR_DISTANCE_FOR_LANE_CHANGE 10 : It defines clear distance (without any vehicle) in front of our car, This distance of 30m is used at the time of making lane change decision.

 9- BUFFER_DISTANCE_TO_FRONT_VEHICLE 30 : It defines the distance to keep from a vehicle in front.

10- NUM_WAYPOINTS 50 : Total number of waypoints sent to simulator is 50. which is enough to have smooth trajectory and to cop up with changing postion of surrounding vehicles.

### Prediction
Through using the car's state and sensor fusion data, 
possible collisions are predicted. Frenet coordinates of the surrounding vehicles are used to determine their respective
lanes and distance from the car. The outputs of this step are array of flags representing occupancy of each lane 
based on current and predicted trajectory of Ego Vehicle and other target vehicles.  

### Behavior Planning 
Based on the previous predictions, 
corrective behaviour is needed. The method used can be considered a state machine with 3 states 
( Keep lane, change left and change right). 
Cost functions are considered too complex for such a simple application. However, 
a cost function is essential when optimized lane changes to improve forward progress is needed,
 which is a considerable improvement to the current implementation.
 
 If we have a slow moving vehicle in front of us. I am trying to check a free left and right Lane. 
 If there are no vehicle in 30m front and 10m back , I am changing the lane. 
 If there are any I am following the front car with slower speed till we find a lane to change

### Trajectory Generation 
This is where the next target waypoints are calculated and passed to the simulator. 
The first two points in the path are the last two points in the previous path. 
Given the target lane, three evenly spaced points (30 meters apart) are added to the path. 
A spline is generated through these points using the spline.h library. 
Waypoints are extracted as the first 30 meters of the spline and are sent to the simulator. 
A distance of 0.02 * ref_velocity / 2.24 is kept between the path points to maintain a velocity below the speed limit (50 mph).
Division by 2.24 converts from miles per hour to meters per second.

   
### Simulator.
You can download the Term3 Simulator which contains the Path Planning Project from the [releases tab (https://github.com/udacity/self-driving-car-sim/releases/tag/T3_v1.2).

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!


## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).

