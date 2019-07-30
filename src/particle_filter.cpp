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
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  // Set the number of particles
  num_particles = 105;
  // resize weights vector for resample use
  weights.resize(num_particles);
  // start a random engine
  std::default_random_engine gen;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; ++i){

    // sample x,y,theta from normal distribution

    // particle data structure
    Particle instanceP;
    instanceP.id = i;
    instanceP.x = dist_x(gen);
    instanceP.y = dist_y(gen);
    instanceP.theta = dist_theta(gen);
    instanceP.weight = 1.0f;
    
    // adding p to particle vector
    particles.push_back(instanceP);

  }
  // initial is true
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
  normal_distribution<double> x_noise(0.0, std_pos[0]);
  normal_distribution<double> y_noise(0.0, std_pos[1]);
  normal_distribution<double> t_noise(0.0, std_pos[2]);

  for (int i=0; i < particles.size(); ++i){

    // prediction
    // keep straight if yaw_rate is too small, keep stable
    if (fabs(yaw_rate) < 0.00001){
      particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
    } 
    else{
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + x_noise(gen);
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + y_noise(gen);
      particles[i].theta += yaw_rate * delta_t + t_noise(gen);
    }
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // for (int i=0; i < observations.size(); ++i){

  //   double d_min = 100000.0;
  //   int closest_id;

  //   for (int j=0; j < predicted.size(); ++j){

  //     d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y)

  //     if (d < d_min){
  //       d_min = d;
  //       closest_id = predicted.id;
  //     }
  //   }
  //   observations[i].id = closest_id;
  // }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double std_landmark_x = std_landmark[0];
  double std_landmark_y = std_landmark[1];
  double gussian_norm =  1.0 / (std_landmark_x*std_landmark_y*2.0*M_PI);

  for (int p=0; p < particles.size(); ++p){

    double px = particles[p].x;
    double py = particles[p].y;
    double ptheta = particles[p].theta;
    // reset weights for each particle
    particles[p].weight = 1.0;

    // create predictions vector
    vector<LandmarkObs> predicteds;
    // filter map landmarks within sensor range
    for (int m=0; m < map_landmarks.landmark_list.size(); ++m){

      int maplm_i = map_landmarks.landmark_list[m].id_i;
      double maplm_x = map_landmarks.landmark_list[m].x_f;
      double maplm_y = map_landmarks.landmark_list[m].y_f;

      if (fabs(px-maplm_x)<=sensor_range && fabs(py-maplm_y)<=sensor_range){
        predicteds.push_back(LandmarkObs{maplm_i, maplm_x, maplm_y});
      }

    }
    // create converted map coordinate observation vector
    // vector<LandmarkObs> t_observations;
    long double w = 1.0;

    for (int j=0; j < observations.size(); ++j){
      // sensor sense obs in vehicle coordinate, need to convert to map coordinate
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      double trans_x = px + cos(ptheta) * x_obs - sin(ptheta) * y_obs;
      double trans_y = py + sin(ptheta) * x_obs + cos(ptheta) * y_obs;

      double d_min = 10000.0;
      double pred_x_min, pred_y_min;

      // find the closest landmarks
      for (int k=0; k < predicteds.size(); ++k){

        double pred_x = predicteds[k].x;
        double pred_y = predicteds[k].y;

        double d = dist(pred_x, pred_y, trans_x, trans_y);

        if (d < d_min){
          d_min = d;
          pred_x_min = pred_x;
          pred_y_min = pred_y;
        }
      }

      double current_weight = gussian_norm * exp(-(0.5*pow((trans_x-pred_x_min)/std_landmark_x, 2.0)+0.5*pow((trans_y-pred_y_min)/std_landmark_y, 2.0)));
      w *= current_weight;
    }

    particles[p].weight = w;
    weights[p] = w;

  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;

  std::discrete_distribution<> discrete_dist(weights.begin(), weights.end());
  std::vector<Particle> particles_sampled;

  for (int i=0; i < particles.size(); ++i){
    // sample more index with higher prob
    int id_sampled = discrete_dist(gen);
    particles_sampled.push_back(particles[id_sampled]);
    
  }
  // update the particles set
  particles = particles_sampled;
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