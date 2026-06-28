#ifndef ISOTOPE_H
#define ISOTOPE_H

#include <string>
#include <cmath>

namespace PET_MC {

class Isotope {
public:
    // Constructor with isotope name
    Isotope(const std::string& name);
    
    // Isotope properties
    std::string getName() const { return name_; }
    int getZ() const { return Z_; }
    int getA() const { return A_; }
    double getHalfLife() const { return halfLife_; }        // seconds
    double getDecayConstant() const { return decayConstant_; } // 1/s
    bool isPositronEmitter() const { return isPET_; }
    bool isBetaMinus() const { return isBetaMinus_; }
    double getQValue() const { return Q_value_; } // MeV
    
    // Activity calculation
    // N(t) = N0 * exp(-lambda * t)
    double getActivityAtTime(double initialActivity,  // Bq
                            double timeElapsed) const; // seconds
    
    // Number of decays between times t1 and t2
    double getDecayedCount(double initialCount,
                          double t1, double t2) const;
    
    // Branching ratio to positron emission
    double getPositronBranchingRatio() const { return positronBranching_; }
    
    // Mean positron energy (for Doppler correction)
    double getMeanPositronEnergy() const { return meanPositronEnergy_; } // MeV
    
    // String representation
    std::string toString() const;
    
    // Database methods
    static Isotope* getIsotope(const std::string& name);
    static bool registerIsotope(const std::string& name,
                               int Z, int A, double halfLife_s,
                               bool isPET, double Q_MeV,
                               double positronBranching = 1.0);
    
private:
    std::string name_;
    int Z_;                     // Atomic number
    int A_;                     // Mass number
    double halfLife_;           // seconds
    double decayConstant_;      // 1/s, ln(2)/halfLife
    bool isPET_;                // Positron emitter?
    bool isBetaMinus_;          // Beta-minus decay?
    double Q_value_;            // Q-value in MeV
    double positronBranching_;  // Branching ratio to β+
    double meanPositronEnergy_; // MeV (typical ~0.64 MeV for C-11)
};

// Common PET isotopes
inline void initializeCommonIsotopes() {
    // C-11: ¹¹C → ¹¹B + e⁺ + νₑ, half-life = 20.4 min
    Isotope::registerIsotope("C11", 6, 11, 20.4*60, true, 1.982, 1.0);
    
    // O-15: ¹⁵O → ¹⁵N + e⁺ + νₑ, half-life = 2.03 min
    Isotope::registerIsotope("O15", 8, 15, 2.03*60, true, 1.732, 1.0);
    
    // F-18: ¹⁸F → ¹⁸O + e⁺ + νₑ, half-life = 110 min
    Isotope::registerIsotope("F18", 9, 18, 110*60, true, 0.635, 0.97);
    
    // N-13: ¹³N → ¹³C + e⁺ + νₑ, half-life = 10.0 min
    Isotope::registerIsotope("N13", 7, 13, 10.0*60, true, 1.944, 1.0);
}

} // namespace PET_MC

#endif // ISOTOPE_H
