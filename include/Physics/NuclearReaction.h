#ifndef NUCLEAR_REACTION_H
#define NUCLEAR_REACTION_H

#include "Transport/Particle.h"
#include "Physics/CrossSection.h"
#include <random>
#include <vector>

namespace PET_MC {

// Result of a nuclear reaction
struct ReactionEvent {
    bool occurred;                      // Did reaction happen?
    std::string productName;            // e.g., "C11"
    bool isPositronEmitter;
    double reactionEnergy;              // Energy at which reaction occurred
    Vector3D reactionPosition;          // Where reaction happened
    std::vector<Particle> secondaries;  // Secondary particles produced
};

class NuclearReaction {
public:
    // Constructor
    NuclearReaction(unsigned int seed = 0);
    
    // Main interface: sample a reaction
    ReactionEvent sampleReaction(const Particle& incident,
                                const std::string& targetMaterial,
                                double stepLength) const;
    
    // Probability of at least one reaction occurring
    double getReactionProbability(const Particle& incident,
                                 const std::string& targetMaterial,
                                 double stepLength) const;
    
    // Sample secondary particle kinematics (simplified isotropic)
    void generateSecondaryParticles(const Particle& incident,
                                   const ReactionProduct& product,
                                   ReactionEvent& event) const;
    
    // Get target nucleus info
    struct TargetNucleus {
        int Z;
        int A;
        double atomicMass; // u
    };
    
    TargetNucleus getTargetNucleus(const std::string& material) const;
    
    // Cross-section utility
    double getCrossSection(const Particle& incident,
                          const std::string& targetMaterial) const;
    
private:
    mutable std::mt19937 rng_;
    mutable std::uniform_real_distribution<> uniform_;
    
    // Recoil energy calculation
    double computeRecoilEnergy(const Particle& incident,
                              const ReactionProduct& product) const;
    
    // Forward angle bias (reactions are forward-peaked)
    void biasSecondaryAngle(Particle& secondary,
                           double incidentEnergy) const;
};

} // namespace PET_MC

#endif // NUCLEAR_REACTION_H
