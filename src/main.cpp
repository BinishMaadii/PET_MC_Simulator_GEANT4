#include "Transport/Particle.h"
#include "Transport/Material.h"
#include "Physics/StoppingPower.h"
#include "Physics/CrossSection.h"
#include "Physics/NuclearReaction.h"
#include "PET/Activity.h"
#include "PET/Isotope.h"
#include "Geometry/DetectorGeometry.h"

// ROOT Includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

using namespace PET_MC;

// ============================================================================
// GLOBAL ROOT OBJECTS
// ============================================================================

TFile* rootFile = nullptr;
TH1F* hEnergyLoss = nullptr;
TH1F* hBraggPeak = nullptr;
TH2F* hActivityMap = nullptr;
TH1F* hDecayC11 = nullptr;
TGraph* gCrossSection = nullptr;

// ============================================================================
// INITIALIZATION FUNCTIONS
// ============================================================================

void initializeROOT() {
    std::cout << "\nInitializing ROOT..." << std::endl;
    
    // Create ROOT file
    rootFile = new TFile("PET_Simulation_Results.root", "RECREATE");
    
    // Create histograms
    hEnergyLoss = new TH1F("hEnergyLoss", "Proton Energy Loss in Water;Distance (mm);Energy (MeV)",
                           300, 0, 300);
    hEnergyLoss->SetLineColor(kBlue);
    hEnergyLoss->SetLineWidth(2);
    
    hBraggPeak = new TH1F("hBraggPeak", "Bragg Peak (dE/dx vs Distance);Distance (mm);-dE/dx (MeV/mm)",
                          300, 0, 300);
    hBraggPeak->SetLineColor(kRed);
    hBraggPeak->SetLineWidth(2);
    
    hActivityMap = new TH2F("hActivityMap", "Activity Distribution (3D projected to 2D);X (mm);Z (mm)",
                            10, 0, 10, 100, 0, 100);
    
    hDecayC11 = new TH1F("hDecayC11", "C-11 Radioactive Decay;Time (min);Activity (Bq)",
                         101, 0, 100);
    hDecayC11->SetLineColor(kGreen);
    hDecayC11->SetLineWidth(2);
    
    std::cout << "ROOT initialized - file: PET_Simulation_Results.root" << std::endl;
}

void closeROOT() {
    if (rootFile) {
        rootFile->Write();
        rootFile->Close();
        std::cout << "\nROOT file closed: PET_Simulation_Results.root" << std::endl;
    }
}

// ============================================================================
// SIMULATION FUNCTIONS
// ============================================================================

void printWelcome() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  PET ACTIVITY PREDICTION MONTE CARLO SIMULATOR" << std::endl;
    std::cout << "  Proton Therapy Dose Verification" << std::endl;
    std::cout << "  With GEANT4 Physics & ROOT Visualization" << std::endl;
    std::cout << std::string(60, '=') << "\n" << std::endl;
}

void printParticleInfo(const Particle& p) {
    std::cout << "Particle Information:" << std::endl;
    std::cout << "  Energy: " << std::fixed << std::setprecision(3) 
              << p.getEnergy() << " MeV" << std::endl;
    std::cout << "  Beta: " << std::setprecision(4) << p.getBeta() << std::endl;
    std::cout << "  Gamma: " << std::setprecision(4) << p.getGamma() << std::endl;
    std::cout << "  Momentum: " << std::setprecision(3) << p.getMomentum() << " MeV/c" << std::endl;
    std::cout << "  Atomic Number: " << p.getAtomicNumber() << std::endl;
    std::cout << "  Remaining Range (in water): " << p.getRemainingRange() << " mm\n" << std::endl;
}

void runProtonSimulation() {
    std::cout << "\n--- PROTON THERAPY SIMULATION (with Bethe-Bloch) ---\n" << std::endl;
    
    Material::initializeCommonMaterials();
    
    double proton_energy = 150.0;  // MeV
    Vector3D position(0, 0, 0);
    Vector3D direction(0, 0, 1);
    
    Particle proton(ParticleType::PROTON, proton_energy, position, direction);
    
    std::cout << "Created " << proton_energy << " MeV proton" << std::endl;
    printParticleInfo(proton);
    
    StoppingPower stoppingPower;
    
    std::cout << "Energy Loss Calculation (Bethe-Bloch):" << std::endl;
    std::cout << "Material: Water" << std::endl;
    
    double step_size = 0.5;
    double current_energy = proton_energy;
    double distance_traveled = 0;
    
    std::cout << "\nStep-by-step energy loss:" << std::endl;
    std::cout << std::setw(10) << "Distance (mm)" 
              << std::setw(15) << "Energy (MeV)"
              << std::setw(15) << "-dE/dx (MeV/mm)" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    // Store data for ROOT histogram
    std::vector<double> distances;
    std::vector<double> energies;
    std::vector<double> stopping_powers;
    
    while (current_energy > 1.0) {
        Particle temp(ParticleType::PROTON, current_energy, position, direction);
        double dE_dx = stoppingPower.computeEnergyLoss(temp, "water");
        
        // Store for plotting
        distances.push_back(distance_traveled);
        energies.push_back(current_energy);
        stopping_powers.push_back(dE_dx);
        
        // Print every 10mm
        if (std::fmod(distance_traveled, 10.0) < step_size) {
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(10) << distance_traveled
                      << std::setw(15) << current_energy
                      << std::setw(15) << dE_dx << std::endl;
        }
        
        // Fill histograms
        hEnergyLoss->Fill(distance_traveled, current_energy);
        hBraggPeak->Fill(distance_traveled, dE_dx);
        
        current_energy -= dE_dx * step_size;
        distance_traveled += step_size;
        
        if (distance_traveled > 300) break;
    }
    
    double bragg_peak_range = distance_traveled;
    std::cout << "\nBragg Peak Range: " << std::fixed << std::setprecision(2) 
              << bragg_peak_range << " mm" << std::endl;
    
    // Cross-section information
    std::cout << "\n--- NUCLEAR REACTION CROSS-SECTIONS ---\n" << std::endl;
    
    CrossSection& xsdb = g_crossSectionDB;
    std::vector<std::string> reactions = xsdb.getAvailableReactions();
    
    std::cout << "Available Reactions:" << std::endl;
    for (const auto& rxn : reactions) {
        std::cout << "  - " << rxn << std::endl;
    }
    
    std::cout << "\nProton + Carbon-12 Cross-section (p+C12 → C11+n):" << std::endl;
    std::cout << std::setw(15) << "Energy (MeV)"
              << std::setw(20) << "σ (millibarn)" << std::endl;
    std::cout << std::string(35, '-') << std::endl;
    
    // Create cross-section graph for ROOT
    std::vector<double> energies_xs;
    std::vector<double> sigmas;
    
    for (double E = 50; E <= 200; E += 10) {
        double sigma = xsdb.getCrossSection("p+C12", E);
        energies_xs.push_back(E);
        sigmas.push_back(sigma);
        
        if (std::fmod(E - 50, 30.0) < 0.1) {  // Print every 30 MeV
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(15) << E
                      << std::setw(20) << sigma << std::endl;
        }
    }
    
    // Create ROOT graph for cross-sections
    gCrossSection = new TGraph(energies_xs.size(), energies_xs.data(), sigmas.data());
    gCrossSection->SetTitle("p+C12 Cross-section;Energy (MeV);σ (millibarn)");
    gCrossSection->SetLineColor(kMagenta);
    gCrossSection->SetLineWidth(2);
    gCrossSection->SetMarkerColor(kMagenta);
    gCrossSection->SetMarkerStyle(20);
}

void runPETActivitySimulation() {
    std::cout << "\n--- PET ACTIVITY PREDICTION SIMULATION ---\n" << std::endl;
    
    initializeCommonIsotopes();
    
    std::cout << "Initialized PET Isotopes:" << std::endl;
    std::vector<std::string> isotope_names = {"C11", "O15", "F18", "N13"};
    
    for (const auto& name : isotope_names) {
        Isotope iso(name);
        std::cout << std::setw(8) << iso.getName()
                  << " | Z=" << iso.getZ() << ", A=" << iso.getA()
                  << " | Half-life: " << std::setprecision(2) 
                  << iso.getHalfLife() / 60.0 << " min"
                  << " | λ = " << iso.getDecayConstant() << " s⁻¹" << std::endl;
    }
    
    // Simulate decay of C-11
    std::cout << "\n--- C-11 Decay Simulation ---" << std::endl;
    Isotope c11("C11");
    double initial_activity = 1000.0;  // Bq
    
    std::cout << "\nActivity vs Time (C-11):" << std::endl;
    std::cout << std::setw(15) << "Time (min)"
              << std::setw(20) << "Activity (Bq)"
              << std::setw(20) << "Activity (% initial)" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (double t_min = 0; t_min <= 100; t_min += 10) {
        double t_sec = t_min * 60;
        double activity = c11.getActivityAtTime(initial_activity, t_sec);
        double percent = (activity / initial_activity) * 100.0;
        
        // Fill ROOT histogram
        hDecayC11->Fill(t_min, activity);
        
        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(15) << t_min
                  << std::setw(20) << activity
                  << std::setw(20) << percent << std::endl;
    }
}

void runVoxelScoringExample() {
    std::cout << "\n--- VOXEL SCORING & ACTIVITY DISTRIBUTION ---\n" << std::endl;
    
    Activity activity_map(10, 10, 10, 1.0);
    
    std::cout << "Created 10×10×10 voxel grid (1 mm³ voxels)" << std::endl;
    std::cout << "\nSimulating activity deposition from particle track..." << std::endl;
    
    // Record C-11 production along a track
    for (double z = 0; z < 100; z += 2) {
        double intensity = 100.0 * std::exp(-std::pow(z - 80, 2) / (2 * 10*10));
        for (int i = 0; i < static_cast<int>(intensity); i++) {
            Vector3D pos(5.0, 5.0, z);
            activity_map.recordIsotope("C11", pos);
            
            // Fill 2D histogram
            hActivityMap->Fill(5.0, z, 1.0);
        }
    }
    
    std::cout << "Total activity recorded: " 
              << activity_map.getTotalActivity("C11") 
              << " (arbitrary units)" << std::endl;
    
    Vector3D max_pos = activity_map.getMaxActivityPosition();
    std::cout << "Max activity position: (" << max_pos.x << ", " 
              << max_pos.y << ", " << max_pos.z << ")" << std::endl;
}

void printInstructions() {
    std::cout << "\n--- NEXT STEPS ---\n" << std::endl;
    std::cout << "1. Visualization (ROOT files):" << std::endl;
    std::cout << "   root PET_Simulation_Results.root" << std::endl;
    std::cout << "   root> hEnergyLoss->Draw()" << std::endl;
    std::cout << "   root> hBraggPeak->Draw()" << std::endl;
    std::cout << "   root> hDecayC11->Draw()" << std::endl;
    std::cout << "\n2. Further extensions:" << std::endl;
    std::cout << "   - Full GEANT4 tracking engine" << std::endl;
    std::cout << "   - Realistic nuclear physics lists" << std::endl;
    std::cout << "   - Multi-particle batch simulation" << std::endl;
    std::cout << "   - Advanced ROOT analysis" << std::endl;
    std::cout << "\n3. Clinical integration:" << std::endl;
    std::cout << "   - DICOM CT import" << std::endl;
    std::cout << "   - Patient-specific simulations" << std::endl;
    std::cout << "   - Treatment plan verification" << std::endl;
    std::cout << std::endl;
}

int main() {
    try {
        printWelcome();
        
        // Initialize ROOT
        initializeROOT();
        
        // Run simulations
        runProtonSimulation();
        runPETActivitySimulation();
        runVoxelScoringExample();
        
        printInstructions();
        
        // Close ROOT file
        closeROOT();
        
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "  Simulation completed successfully!" << std::endl;
        std::cout << "  Results saved to: PET_Simulation_Results.root" << std::endl;
        std::cout << std::string(60, '=') << "\n" << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
