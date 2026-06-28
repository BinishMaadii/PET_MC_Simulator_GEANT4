#include "Transport/Particle.h"
#include <cmath>
#include <stdexcept>

namespace PET_MC {

// Static mass and charge data
static const double PROTON_MASS = 938.272;      // MeV/c²
static const double ALPHA_MASS = 3727.379;      // MeV/c²
static const double CARBON12_MASS = 11177.929;  // MeV/c²
static const double ELECTRON_MASS = 0.511;      // MeV/c²

Particle::Particle(ParticleType type, double energy, Vector3D position, Vector3D direction)
    : type_(type), energy_(energy), position_(position), direction_(direction.normalized()) {
    
    // Set particle properties based on type
    switch (type) {
        case ParticleType::PROTON:
            mass_ = PROTON_MASS;
            charge_ = 1.0;
            Z_ = 1;
            A_ = 1;
            break;
        case ParticleType::ALPHA:
            mass_ = ALPHA_MASS;
            charge_ = 2.0;
            Z_ = 2;
            A_ = 4;
            break;
        case ParticleType::CARBON_12:
            mass_ = CARBON12_MASS;
            charge_ = 6.0;
            Z_ = 6;
            A_ = 12;
            break;
        case ParticleType::ELECTRON:
        case ParticleType::POSITRON:
            mass_ = ELECTRON_MASS;
            charge_ = -1.0;
            Z_ = 0;
            A_ = 0;
            break;
        default:
            throw std::runtime_error("Unknown particle type");
    }
}

double Particle::getRemainingRange() const {
    // Simplified range approximation (Geisler formula for protons in water)
    if (type_ != ParticleType::PROTON) return 0;
    
    if (energy_ < 0) return 0;
    if (energy_ < 2.0) {
        return 0.56 * energy_;  // mm in water, below 2 MeV
    } else if (energy_ < 15.0) {
        return 0.31 * energy_ - 0.06;  // mm
    } else {
        return 1.8 * energy_ - 2.2;    // mm
    }
}

double Particle::getGamma() const {
    // γ = 1 + (KE / m*c²)
    return 1.0 + energy_ / mass_;
}

double Particle::getBeta() const {
    // β = v/c = sqrt(1 - 1/γ²)
    double gamma = getGamma();
    return std::sqrt(1.0 - 1.0/(gamma*gamma));
}

double Particle::getMomentum() const {
    // p*c = sqrt((E_total)² - (m*c²)²)
    double E_total = energy_ + mass_;
    double pc_squared = E_total*E_total - mass_*mass_;
    if (pc_squared < 0) return 0;
    return std::sqrt(pc_squared);
}

} // namespace PET_MC
