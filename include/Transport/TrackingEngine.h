#ifndef TRACKING_ENGINE_H
#define TRACKING_ENGINE_H

#include "Transport/Particle.h"
#include "Transport/Material.h"
#include "Physics/StoppingPower.h"
#include "Physics/NuclearReaction.h"
#include "Physics/CrossSection.h"
#include "PET/Activity.h"
#include <vector>
#include <memory>

namespace PET_MC {

// Simulation configuration
struct SimulationConfig {
    double stepSize;            // mm (step length for tracking)
    double minEnergy;           // MeV (threshold to stop tracking)
    int maxSteps;               // Maximum steps per particle
    int randomSeed;
    bool recordTrajectory;      // Store full path?
    bool computeSecondaries;    // Track secondary particles?
    std::string outputPrefix;   // Output file prefix
};

// Simulation statistics
struct SimulationStats {
    int particlesTracked;
    int reactionEvents;
    int positronEmittersProduced;
    double totalDoseDeposited;  // Gy
    double computeTime;         // seconds
};

class TrackingEngine {
public:
    // Constructor
    TrackingEngine(const SimulationConfig& config);
    
    // Initialize geometry (material map)
    bool loadGeometry(const std::string& geomFile);
    
    // Define a simple slab geometry
    void createSlabGeometry(int thickness_mm, 
                           const std::string& material = "water");
    
    // Track a single particle
    void trackParticle(const Particle& particle);
    
    // Run batch simulation
    void runSimulation(const std::vector<Particle>& particles);
    
    // Generate a pencil beam (common clinical setup)
    std::vector<Particle> generatePencilBeam(
        ParticleType type,
        double energy_MeV,
        int num_particles,
        double spotSigma_mm = 2.0);  // Gaussian beam sigma
    
    // Get activity distribution
    const Activity& getActivityDistribution() const { return activity_; }
    Activity& getActivityDistribution() { return activity_; }
    
    // Get statistics
    const SimulationStats& getStats() const { return stats_; }
    
    // Export results
    void exportResults() const;
    
    // Set material map directly
    void setMaterialAtPosition(const Vector3D& pos, const std::string& material);
    
private:
    SimulationConfig config_;
    StoppingPower stoppingPower_;
    NuclearReaction nuclearReaction_;
    CrossSection crossSections_;
    Activity activity_;
    SimulationStats stats_;
    
    // Material map
    std::map<Vector3D, std::string> materialMap_;  // Position -> material name
    
    // Track a particle until it stops
    void trackToStop(Particle& particle);
    
    // Single tracking step
    void step(Particle& particle, double stepLength);
    
    // Get material at position (nearest neighbor)
    std::string getMaterialAtPosition(const Vector3D& pos) const;
    
    // Boundary handling
    bool checkBoundary(const Vector3D& pos) const;
    
    // Time tracking (for decay simulation)
    double totalSimulationTime_;
};

} // namespace PET_MC

#endif // TRACKING_ENGINE_H
