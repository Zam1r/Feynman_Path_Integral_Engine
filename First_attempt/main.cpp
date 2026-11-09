#include <vector> 
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream> 

struct QuantumPath {
    std::vector<double> positions; 
    double delta_t;
};

double calculate_action(const QuantumPath& path, double mass, double omega) {
    double action = 0.0;
    int num_points = path.positions.size();
  
    for (int i = 0; i < num_points - 1; i = i + 1) {
        double dx = path.positions[i+1] - path.positions[i];
        double kinetic = 0.5 * mass * (dx * dx) / path.delta_t;
        double current_x = path.positions[i];
        double potential = 0.5 * mass * (omega * omega) * (current_x * current_x) * path.delta_t;
        action += (kinetic + potential);
    }

    return action; 
}

void mutate_path(QuantumPath& path, double max_shift) {
    int num_points = path.positions.size();
    int j = 1 + (rand() % (num_points - 2));
    double random_fraction = (2.0 * rand() / RAND_MAX) - 1.0; 
    path.positions[j] = path.positions[j] + (random_fraction * max_shift);
}

void monte_carlo_step(QuantumPath& path, double mass, double omega, double max_shift, double hbar) {
    double old_action = calculate_action(path, mass, omega);
    QuantumPath new_path = path; 
    mutate_path(new_path, max_shift);
    double new_action = calculate_action(new_path, mass, omega);

    if (new_action < old_action) {
        path = new_path; 
    } else {
        double P = std::exp(-(new_action - old_action) / hbar); 
        double random_decimal = (double)rand() / RAND_MAX; 
        if (random_decimal < P) {
            path = new_path; 
        }
    }
}

int main() {
    QuantumPath path;
    path.positions = std::vector<double>(100, 0.0);
    path.delta_t = 0.1;
    
    double mass = 1.0;
    double omega = 1.0;
    double max_shift = 0.5;
    double hbar = 1.0;
    
    std::cout << "Starting thermalization..." << std::endl;
    int thermal_steps = 100000;
    
    for (int i = 0; i < thermal_steps; i = i + 1) {
        monte_carlo_step(path, mass, omega, max_shift, hbar);
    }
    
    std::ofstream outfile("histogram_data.txt");
    
    std::cout << "Starting data collection..." << std::endl;
    int production_steps = 900000;
    int mid_index = path.positions.size() / 2; // This equals 50
    
    for (int i = 0; i < production_steps; i = i + 1) {
        monte_carlo_step(path, mass, omega, max_shift, hbar);
        
        outfile << path.positions[mid_index] << "\n";
    }
    
    outfile.close();
    std::cout << "Simulation complete! Data saved to histogram_data.txt" << std::endl;
    
    return 0;
}
