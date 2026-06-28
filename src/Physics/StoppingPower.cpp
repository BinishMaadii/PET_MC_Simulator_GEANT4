#include "Physics/StoppingPower.h"
#include <cmath>
#include <stdexcept>

namespace PET_MC {

StoppingPower::StoppingPower() {}

double StoppingPower::computeEnergyLoss(const Particle& particle, 
                                       const std::string& material) const {
    // Get particle properties
    double Z = particle.getAtomicNumber();
    double mass = particle.getMass();
    double energy = particle.getEnergy();
    double beta = particle.getBeta();
    double gamma = particle.getGamma();
    
    if (energy <= 0) return 0;
    
    // For now, use water as reference material
    // In full implementation, would have material database
    int Z_target = 7;  // Oxygen in water
    double A_target = 16;
    double rho = 1.0;  // g/cm³
    
    if (material != "water") {
        // Could add other materials here
        Z_target = 7;
        A_target = 16;
        rho = 1.0;
    }
    
    // Bethe-Bloch formula: -dE/dx = K * z² * Z_target / A_target * 
    //                               1/β² * [ln(2*m*β²*γ²*T_max/I²) - 2*β² - δ]
    
    double z_squared = Z * Z;
    double beta_squared = beta * beta;
    
    // Electron rest mass
    double m_e = 0.511;  // MeV/c²
    
    // Maximum energy transfer to electron
    double T_max = (2.0 * m_e * beta_squared * gamma * gamma) / 
                   (1.0 + 2.0 * gamma * m_e / mass);
    
    // Mean ionization potential (approximation)
    double I = getMeanIonizationEnergy(Z_target);
    
    // Logarithmic term
    double log_term = std::log(2.0 * m_e * beta_squared * gamma * gamma * T_max / 
                              (I * I));
    
    // Density correction (simplified)
    double delta = getDensityCorrection(beta, Z_target);
    
    // Bethe-Bloch formula (MeV/cm)
    double dE_dx = K * z_squared * Z_target / A_target * 
                   (1.0 / beta_squared) * 
                   (log_term - 2.0 * beta_squared - delta);
    
    // Convert from MeV/cm to MeV/mm
    return dE_dx / 10.0;
}

double StoppingPower::getMeanIonizationEnergy(int Z) const {
    // Approximation: I ≈ 9.76 * Z + 58.8 * Z^(-0.19) eV
    // Convert to MeV
    double I_eV = 9.76 * Z + 58.8 * std::pow(Z, -0.19);
    return I_eV / 1e6;  // Convert eV to MeV
}

double StoppingPower::getDensityCorrection(double beta, int Z_target) const {
    // Simplified density correction (usually small for water)
    // Full Sternheimer parametrization would go here
    if (beta < 0.1) return 0;
    return 0.0;  // Negligible for protons in water
}

double StoppingPower::getMultipleScatteringAngle(double distance,
                                                const Particle& particle,
                                                int Z_target) const {
    // Highland formula: θ_rms ≈ (13.6 MeV / p*β*c) * sqrt(t/X0) * (1 + 0.038*ln(t/X0))
    // Simplified version here
    
    double beta = particle.getBeta();
    double momentum = particle.getMomentum();  // MeV/c
    
    if (momentum < 1) return 0;
    
    // Radiation length in water ≈ 36 cm
    double X0_water = 360.0;  // mm
    
    double theta_rms = (13.6 / (momentum * beta)) * 
                      std::sqrt(distance / X0_water) * 
                      (1.0 + 0.038 * std::log(distance / X0_water));
    
    return theta_rms;  // radians
}

double StoppingPower::getRangeFromEnergy(double energy,
                                        const Particle& particle,
                                        const std::string& material) const {
    // Simplified: integrate stopping power
    // For protons in water: rough approximation
    
    if (particle.getAtomicNumber() != 1 || energy <= 0) return 0;
    
    // Geisler-Bethe formula for proton range in water (mm)
    if (energy < 2.0) {
        return 0.56 * energy;
    } else if (energy < 15.0) {
        return 0.31 * energy - 0.06;
    } else {
        return 1.8 * energy - 2.2;
    }
}

double StoppingPower::getEnergyFromRange(double range,
                                        const Particle& particle,
                                        const std::string& material) const {
    // Inverse of getRangeFromEnergy
    if (particle.getAtomicNumber() != 1 || range <= 0) return 0;
    
    if (range < 1.12) {
        return range / 0.56;
    } else if (range < 4.44) {
        return (range + 0.06) / 0.31;
    } else {
        return (range + 2.2) / 1.8;
    }
}

double StoppingPower::betheBlochFormula(double kinetic_energy,
                                       int Z, double mass,
                                       int Z_target, double A_target,
                                       double rho) const {
    // Direct Bethe-Bloch evaluation (MeV/mm)
    if (kinetic_energy <= 0) return 0;
    
    double gamma = 1.0 + kinetic_energy / mass;
    double beta_squared = 1.0 - 1.0 / (gamma * gamma);
    double beta = std::sqrt(beta_squared);
    
    double z_squared = Z * Z;
    double m_e = 0.511;
    double T_max = (2.0 * m_e * beta_squared * gamma * gamma) / 
                   (1.0 + 2.0 * gamma * m_e / mass);
    
    double I = getMeanIonizationEnergy(Z_target);
    double log_term = std::log(2.0 * m_e * beta_squared * gamma * gamma * T_max / 
                              (I * I));
    
    double dE_dx = K * z_squared * (rho / 1.0) * Z_target / A_target * 
                   (1.0 / beta_squared) * (log_term - 2.0 * beta_squared);
    
    return dE_dx / 10.0;  // Convert to MeV/mm
}

} // namespace PET_MC
