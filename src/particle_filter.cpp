/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 20;

	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_p(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		particles.emplace_back(
			i,
			dist_x(gen),
			dist_y(gen),
			dist_p(gen),
			1.0
		);

		weights.push_back(1.0);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	double w = velocity / yaw_rate;
	double hdt = yaw_rate * delta_t;

	std::normal_distribution<double> dist_x(0, std_pos[0]);
	std::normal_distribution<double> dist_y(0, std_pos[1]);
	std::normal_distribution<double> dist_h(0, std_pos[2]);

	for (auto i = particles.begin(), n = particles.end(); i != n; ++i) {
		Particle &p = *i;
		double h0 = p.theta;
		p.x += w * (sin(h0 + hdt) - sin(h0)) + dist_x(gen);
		p.y += w * (cos(h0) - cos(h0 + hdt)) + dist_y(gen);
		p.theta += hdt + dist_h(gen);
	}
}

std::vector<LandmarkObs> ParticleFilter::toMapCoordinates(const Particle &p, const std::vector<LandmarkObs> &observations) {
	std::vector<LandmarkObs> mapped(observations.size());
	double x0 = p.x;
	double y0 = p.y;
	double h = p.theta;
	double cos_h = cos(h);
	double sin_h = sin(h);
	auto j = mapped.begin();
	for (auto i = observations.begin(), n = observations.end(); i != n; ++i, ++j) {
		double xi = i->x;
		double yi = i->y;
		j->x = x0 + xi * cos_h - yi * sin_h;
		j->y = y0 + xi * sin_h + yi * cos_h;
		j->id = i->id;
	}

	return mapped;
}

std::vector<LandmarkObs> ParticleFilter::predict(const std::vector<LandmarkObs> &mapped, const Map &map_landmarks) {
	std::vector<LandmarkObs> predicted(mapped.size());

	auto j = predicted.begin();
	const auto &landmarks = map_landmarks.landmark_list;
	for (auto i = mapped.begin(), m = mapped.end(); i != m; ++i, ++j) {
		double d = std::numeric_limits<double>::max();
		for (auto k = landmarks.begin(), n = landmarks.end(); k != n; ++k) {
			double d_ik = dist(*i, *k);
			if (d > d_ik) {
				d = d_ik;
				j->id = k->id;
				j->x = k->x;
				j->y = k->y;
			}
		}
	}

	return predicted;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; ++i) {
		Particle &p = particles[i];

		// Transform observations from vehicle to map coordinate system.
		std::vector<LandmarkObs> mapped = toMapCoordinates(p, observations);

		// Return a predicted observation for each mapped observation.
		std::vector<LandmarkObs> predicted = predict(mapped, map_landmarks);

		// Update the particle's weight.
		double w = pdf(mapped, predicted, std_landmark);
		p.weight = w;
		weights[i] = w;
	}
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::discrete_distribution<> sampler(weights.begin(), weights.end());

	std::vector<Particle> resampled;
	for (int i = 0; i < num_particles; ++i) {
		Particle &p = particles[sampler(gen)];
		resampled.push_back(p);
	}

	particles = resampled;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
