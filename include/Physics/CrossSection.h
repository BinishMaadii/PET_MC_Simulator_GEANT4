#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <string>
#include <vector>
#include <map>
#include <memory>

namespace PET_MC {

// Reaction product definition
struct ReactionProduct {
    std::string name;           // e.g., "C11", "O15", "F18"
    int Z;                      // Atomic number
    int A;                      // Mass number
    bool isPositronEmitter;     // Relevant for PET
    double branchingRatio;      // Fraction of reactions producing this
    double halfLife;            // seconds (for decay calculation)
};

// Cross-section data point
struct CrossSectionPoint {
    double energy;              // MeV
    double crossSection;        // millibarn (mb)
};

class CrossSection {
public:
    // Constructor
    CrossSection();
    
    // Load cross-section data from file
    bool loadFromFile(const std::string& filename);
    
    // Get cross-section at arbitrary energy (interpolation)
    // Reaction type: "p+C12", "p+O16", "He+C12", etc.
    double getCrossSection(const std::string& reaction, double energy) const;
    
    // Get list of products for a reaction
    std::vector<ReactionProduct> getProducts(const std::string& reaction) const;
    
    // Register a reaction channel with cross-section data
    void registerReaction(const std::string& name,
                         const std::vector<CrossSectionPoint>& data,
                         const std::vector<ReactionProduct>& products);
    
    // Get total cross-section summed over all products
    double getTotalCrossSection(const std::string& reaction, double energy) const;
    
    // Check if reaction is known
    bool hasReaction(const std::string& reaction) const;
    
    // List all registered reactions
    std::vector<std::string> getAvailableReactions() const;
    
private:
    // Data structure for reaction cross-sections
    struct ReactionData {
        std::vector<CrossSectionPoint> crossSectionData;
        std::vector<ReactionProduct> products;
    };
    
    std::map<std::string, ReactionData> reactions_;
    
    // Linear interpolation between energy points
    double interpolate(const std::vector<CrossSectionPoint>& data, double energy) const;
    
    // Add default/built-in cross-sections
    void initializeDefaultCrossSections();
};

// Global cross-section database accessor
extern CrossSection g_crossSectionDB;

} // namespace PET_MC

#endif // CROSS_SECTION_H
