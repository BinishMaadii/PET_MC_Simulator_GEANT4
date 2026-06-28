#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cmath>

namespace PET_MC {

// Particle types
enum class ParticleType {
    PROTON,
    ALPHA,
    CARBON_12,
    ELECTRON,
    POSITRON,
    PHOTON,
    NEUTRON
};

// 3D Vector for position/momentum
struct Vector3D {
    double x, y, z;
    
    Vector3D(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
    
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    
    Vector3D normalized() const {
        double mag = magnitude();
        if (mag < 1e-10) return Vector3D(0, 0, 1);
        return Vector3D(x/mag, y/mag, z/mag);
    }
};

// Main Particle class
class Particle {
public:
    Particle(ParticleType type, double energy, Vector3D position, Vector3D direction);
    
    // Getters
    ParticleType getType() const { return type_; }
    double getEnergy() const { return energy_; }
    double getMass() const { return mass_; }
    double getCharge() const { return charge_; }
    int getAtomicNumber() const { return Z_; }
    int getMassNumber() const { return A_; }
    
    Vector3D getPosition() const { return position_; }
    Vector3D getDirection() const { return direction_; }
    
    // Setters
    void setEnergy(double energy) { energy_ = energy; }
    void setPosition(const Vector3D& pos) { position_ = pos; }
    void setDirection(const Vector3D& dir) { direction_ = dir.normalized(); }
    
    // Track history
    void recordStep(const Vector3D& pos, double dE) {
        trajectory_.push_back(pos);
        energy_losses_.push_back(dE);
    }
    
    const std::vector<Vector3D>& getTrajectory() const { return trajectory_; }
    const std::vector<double>& getEnergyLosses() const { return energy_losses_; }
    
    // Check if particle is alive
    bool isAlive() const { return energy_ > 0; }
    double getRemainingRange() const; // In mm
    
    // Relativistic properties
    double getGamma() const;
    double getBeta() const;
    double getMomentum() const; // MeV/c
    
private:
    ParticleType type_;
    double energy_;           // MeV
    double mass_;             // MeV/c^2
    double charge_;           // Elementary charges
    int Z_;                   // Atomic number
    int A_;                   // Mass number
    
    Vector3D position_;       // mm
    Vector3D direction_;      // Unit vector
    
    std::vector<Vector3D> trajectory_;
    std::vector<double> energy_losses_;
};

} // namespace PET_MC

#endif // PARTICLE_H
