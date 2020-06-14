/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::default_random_engine;
using std::normal_distribution;
using std::discrete_distribution;
using std::numeric_limits;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 70;  // Set the number of particles (can't be too high, so alogrithm does not take too long to compute)
  
  weights = std::vector<double>(static_cast<unsigned long>(num_particles), 1.0);
  
  // Set standard deviations for x, y and theta
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  default_random_engine gen;
  
  // Create normal distributions for x, y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  // initialize particles and add gaussian noise
  
  particles = std::vector<Particle>(static_cast<unsigned long>(num_particles));
  for (auto i = 0; i < num_particles; ++i){   
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = weights[i];
  }
  
  // mark initialization as finished
  is_initialized = true;
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   */
  
  default_random_engine gen;
  // preparing gaussian noise
  normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);

  

  for (auto i = 0; i < num_particles; ++i){   
    // predict using the process model and adding gaussian noise
    // handle case of very small yaw rate respectively constant yaw angle
    if (fabs(yaw_rate) < 0.000001) {
    	particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
    	particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
    } else {
	particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)) + dist_x(gen);
	particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)) + dist_y(gen);
	    particles[i].theta = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);
    }
  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   */
  

  for (unsigned int i=0; i<observations.size(); ++i) {
    double min_distance = numeric_limits<double>::max();
    int landmark_id = -1;
    LandmarkObs obs = observations[i];
    for (unsigned int j=0; j< predicted.size(); ++j) {
      LandmarkObs pred = predicted[j];
      double act_distance = dist(pred.x, pred.y, obs.x, obs.y);
      // find Landmark that is closest to measurement
      if (act_distance <= min_distance) {
         min_distance = act_distance;
         landmark_id = pred.id;
      }
      }
    observations[i].id = landmark_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. 
   */
    
  for (unsigned int i = 0; i < particles.size(); ++i) {
     Particle const &particle = particles[i];
       
     vector<LandmarkObs> obs_map(observations.size());
     for (unsigned int j = 0; j < observations.size(); ++j) {
            LandmarkObs obs = observations[j];
            // transform to map coordinates
            obs_map[j].x = particle.x + (cos(particle.theta) * obs.x) - (sin(particle.theta) * obs.y);
            obs_map[j].y = particle.y + (sin(particle.theta) * obs.x) + (cos(particle.theta) * obs.y);
            obs_map[j].id = -1;
          }
    
    vector<LandmarkObs> landmarks_sensorrange;
    for (auto const &landmark : map_landmarks.landmark_list) {
       // get only landmarks within sensor range
       if (dist(particle.x, particle.y, landmark.x_f, landmark.y_f) <= sensor_range) {
           landmarks_sensorrange.push_back(LandmarkObs{ landmark.id_i, landmark.x_f, landmark.y_f});
       }
    }
    
    dataAssociation(landmarks_sensorrange, obs_map);
    
    vector<double> particle_weights(obs_map.size());
    particles[i].weight = 1.0;
    for (unsigned int l = 0; l < observations.size(); ++l) {
      
      LandmarkObs trans_obs = obs_map[l];
      
      LandmarkObs closest_landmark = {
        .id = -1,
        .x = static_cast<double>(map_landmarks.landmark_list[trans_obs.id - 1].x_f),
        .y = static_cast<double>(map_landmarks.landmark_list[trans_obs.id - 1].y_f),
      };
      
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double x_obs = trans_obs.x;
      double y_obs = trans_obs.y;
      double mu_x = closest_landmark.x;
      double mu_y = closest_landmark.y;
      
      // calculate normalization term
      double gauss_norm;
      gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

      // calculate exponent
      double exponent;
      exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
                   + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

      // calculate weight using normalization terms and exponent
      double weight;
      weight = gauss_norm * exp(-exponent);
      
      particles[i].weight *= weight;
    }
    
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional 
   *   to their weight. 
   */
  default_random_engine gen;
  discrete_distribution<size_t> dist_particle(weights.begin(), weights.end());
  
  vector<Particle> particles_resampled(particles.size());
  
  for (unsigned int i = 0; i < particles.size(); i++) {
  	particles_resampled[i] = particles[dist_particle(gen)];
  }
  
  particles = particles_resampled;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
