#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

//---------------- Defining Constant --------------------------------

#define NUM_OF_LANES 3
#define LANE_WIDTH 4.0
#define MPH_TO_MPS 0.44704
#define TIME_BETWEEN_WAYPOINTS 0.02
#define MAX_VELOCITY_IN_MPH 45
//	#define MIN_VELOCITY_IN_MPH 30
#define MAX_ACCELERATION 7 //in mps^2 kept it 7 to have some before max allowed value in simulator is 10mps^2
#define FRONT_CLEAR_DISTANCE_FOR_LANE_CHANGE 30
#define REAR_CLEAR_DISTANCE_FOR_LANE_CHANGE 10
#define BUFFER_DISTANCE_TO_FRONT_VEHICLE 30
#define NUM_WAYPOINTS 50

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  float ref_vel = 0.0;
  int lane = 1 ;
  float target_velocity = MAX_VELOCITY_IN_MPH;


  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&target_velocity, &ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

			vector<double> anchor_ptsx, anchor_ptsy;
			
			
			int prev_path_size = previous_path_x.size();
			
			
			//reference x,y, yaw states
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);
			bool is_lane_safe_to_change[] = {true,true,true};
			bool is_slow_car_in_front = false;
			
			lane = ((int)car_d/LANE_WIDTH);
			
			if(prev_path_size>0){
				car_s = end_path_s;
			}
			
			// Deciding Safe lanes to change
			for(int i=0 ; i<sensor_fusion.size(); i++){
				float d_other = sensor_fusion[i][6];
				float vx = sensor_fusion[i][3];
				float vy = sensor_fusion[i][4];
				
				double speed_other = sqrt(vx*vx+vy*vy);
				double car_s_other = sensor_fusion[i][5];
				
				//predicting the future position of car based on last point in previous path
				car_s_other += (prev_path_size*TIME_BETWEEN_WAYPOINTS*speed_other);
				
				
				//checking if car is in our lane
				if(d_other >=(LANE_WIDTH*lane) && d_other <= (LANE_WIDTH*(lane+1))){
					//cheking if slow moving car is in 30 m front range
					if(car_s_other>car_s && car_s_other<(car_s+BUFFER_DISTANCE_TO_FRONT_VEHICLE)){
						// cout<<"car detected"<<endl;
						if(speed_other < MAX_VELOCITY_IN_MPH){
							is_slow_car_in_front = true;
							target_velocity = speed_other;
							// cout<<"slow car detected"<<endl;
						}
					}
				}
				
				int lane_of_other_car = ((int)d_other/LANE_WIDTH);
				//checking if car is there in front 
				if((car_s_other > car_s - REAR_CLEAR_DISTANCE_FOR_LANE_CHANGE) && (car_s_other < car_s + FRONT_CLEAR_DISTANCE_FOR_LANE_CHANGE)){
					is_lane_safe_to_change[lane_of_other_car] = false;
				}
			}
			
			if(is_slow_car_in_front){
				int left_lane = lane-1;
				int right_lane = lane+1;
				if(left_lane >= 0 && is_lane_safe_to_change[left_lane]){
					lane = left_lane;
				}else if(right_lane < NUM_OF_LANES && is_lane_safe_to_change[right_lane]){
					lane = right_lane;
				}
				else if(ref_vel > target_velocity){
					ref_vel -= MAX_ACCELERATION*TIME_BETWEEN_WAYPOINTS/MPH_TO_MPS;
				}
			}
			else {
				target_velocity = MAX_VELOCITY_IN_MPH;
				if(ref_vel < target_velocity){
					ref_vel += MAX_ACCELERATION*TIME_BETWEEN_WAYPOINTS/MPH_TO_MPS;
				}
			}
			
			float spline_ref_dist = 30.0;
			
			
			if(prev_path_size < 2){
				
				anchor_ptsx.push_back(car_x-cos(ref_yaw));
				anchor_ptsy.push_back(car_y-sin(ref_yaw));
				
				//adding car current postion in global Map coordinates
				anchor_ptsx.push_back(ref_x);
				anchor_ptsy.push_back(ref_y);
			}else{
				float last_1_point_x =  previous_path_x[prev_path_size-1];
				float last_1_point_y =  previous_path_y[prev_path_size-1];
				
				ref_x = last_1_point_x;
				ref_y = last_1_point_y;
				
				float last_2_point_x =  previous_path_x[prev_path_size-2];
				float last_2_point_y =  previous_path_y[prev_path_size-2];
				
				ref_yaw = atan2(last_1_point_y - last_2_point_y, last_1_point_x - last_2_point_x);
				
				anchor_ptsx.push_back(last_2_point_x);
				anchor_ptsy.push_back(last_2_point_y);
				
				anchor_ptsx.push_back(last_1_point_x);
				anchor_ptsy.push_back(last_1_point_y);
			}
			vector<double> next_wp0 = getXY( car_s+spline_ref_dist,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY( car_s+2.0*spline_ref_dist,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY( car_s+3.0*spline_ref_dist,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

			anchor_ptsx.push_back(next_wp0[0]);
			anchor_ptsx.push_back(next_wp1[0]);
			anchor_ptsx.push_back(next_wp2[0]);

			anchor_ptsy.push_back(next_wp0[1]);
			anchor_ptsy.push_back(next_wp1[1]);
			anchor_ptsy.push_back(next_wp2[1]);
			
			for( int i=0; i<anchor_ptsx.size(); i++){
				//shifting anchor points to car's reference frame
				double shift_x = anchor_ptsx[i] - ref_x;
				double shift_y = anchor_ptsy[i] - ref_y;

				anchor_ptsx[i] = (shift_x * cos(- ref_yaw) - shift_y*sin(- ref_yaw));
				anchor_ptsy[i] = (shift_x * sin(- ref_yaw) + shift_y*cos(- ref_yaw));
			}
			
				
			for(int i=0; i<prev_path_size; i++){
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}
			
			double dist_inc = ref_vel*MPH_TO_MPS*TIME_BETWEEN_WAYPOINTS; //in meters
			
			
			tk::spline s;
			s.set_points(anchor_ptsx,anchor_ptsy);
			//---------------------------- Generating Waypoints -----------------------
			for(int i = 0; i < NUM_WAYPOINTS-prev_path_size; i++)
			{	  
				    float way_x = (i+1)*dist_inc;
				    float way_y = s(way_x);
					 
					// cout<<"way_x : "<<way_x<<" way_y : "<<way_y<<endl;
					 
					//shifting to global map cordinate system
					float new_x = way_x*cos(ref_yaw)-way_y*sin(ref_yaw)+ref_x;
					float new_y = way_x*sin(ref_yaw)+way_y*cos(ref_yaw)+ref_y;
					
					// cout<<"new_x : "<<new_x<<" new_y : "<<new_y<<endl;
					 
				    next_x_vals.push_back(new_x);
					next_y_vals.push_back(new_y);
	
			}
			
	

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}