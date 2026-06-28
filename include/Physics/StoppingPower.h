#ifndef STOPPING_POWER_H
#define STOPPING_POWER_H

#include "Transport/Particle.h"
#include <cmath>

namespace PET_MC {

class StoppingPower {
public:
    // Constructor
    StoppingPower();
    
    // Main interface: compute dE/dx in MeV/mm
    // particle: the traveling particle
    // material: target material (water, tissue, etc.)
    double computeEnergyLoss(const Particle& particle, const std::string& material) const;
    
    // Bethe-Bloch formula components
    double betheBlochFormula(double kinetic_energy, 
                            int Z,      // atomic number of projectile
                            double mass, // MeV/c^2
                            int Z_target, // atomic number of target material
                            double A_target, // mass number of target
                            double rho) const; // density g/cm^3
    
    // Mean ionization potential (approximation)
    double getMeanIonizationEnergy(int Z) const;
    
    // Density effect correction
    double getDensityCorrection(double beta, int Z_target) const;
    
    // Multiple scattering angle (Highland formula)
    double getMultipleScatteringAngle(double distance, 
                                     const Particle& particle,
                                     int Z_target) const;
    
    // Range-energy relation (Geisler-Bethe)
    double getRangeFromEnergy(double energy, const Particle& particle, 
                             const std::string& material) const;
    
    // Energy from range (inverse relation)
    double getEnergyFromRange(double range, const Particle& particle,
                             const std::string& material) const;
    
private:
    // Physical constants
    static constexpr double K = 0.307075;      // MeV·cm²/g (Bethe-Bloch constant)
    static constexpr double M_ELECTRON = 0.511; // MeV/c²
    static constexpr double C = 2.998e8;       // m/s (speed of light)
    
    // Water parameters (for reference tissue)
    static constexpr double WATER_Z = 7.42;
    static constexpr double WATER_A = 18.02;
    static constexpr double WATER_RHO = 1.0;  // g/cm³
    static constexpr double WATER_I = 0.0755; // MeV (mean ionization)
    
    // Logarithmic term computation
    double computeLogTerm(double beta, double gamma, int Z) const;
    
    // Shell correction
    double getShellCorrection(double beta, int Z_target) const;
};

} // namespace PET_MC

#endif // STOPPING_POWER_H
