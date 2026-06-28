#ifndef ACTIVITY_H
#define ACTIVITY_H

#include "Transport/Particle.h"
#include "PET/Isotope.h"
#include <map>
#include <vector>
#include <string>

namespace PET_MC {

// Single voxel activity record
struct VoxelActivity {
    int x, y, z;                                // Voxel indices
    std::map<std::string, double> isotopeActivity; // Activity per isotope (Bq)
    double totalActivity;                       // Sum across all isotopes
};

class Activity {
public:
    // Constructor
    Activity(int nx, int ny, int nz,          // Voxel grid dimensions
            double voxel_size_mm = 1.0);      // Voxel edge length in mm
    
    // Record a positron-emitting isotope production
    void recordIsotope(const std::string& isotopeName,
                      const Vector3D& position,
                      double energy = 0.0);  // Optional: energy deposited
    
    // Record energy deposition for dose calculation
    void recordDose(const Vector3D& position, double dE);
    
    // Get activity at a voxel
    double getActivityAtVoxel(int x, int y, int z, 
                            const std::string& isotopeFilter = "") const;
    
    // Get dose at a voxel (in Gy, assuming 1 kg tissue)
    double getDoseAtVoxel(int x, int y, int z) const;
    
    // Apply decay from time t to time t+dt
    void applyDecay(double deltaTime_seconds);
    
    // Apply washout (biological clearance)
    void applyWashout(double washoutFraction, const std::string& tissue = "brain");
    
    // Get total activity (sum all voxels)
    double getTotalActivity(const std::string& isotopeFilter = "") const;
    
    // Export to file (simple ASCII format)
    bool exportToFile(const std::string& filename) const;
    
    // Export to ROOT file (if ROOT available)
    bool exportToROOT(const std::string& filename) const;
    
    // Get statistics
    double getMaxActivity() const;
    Vector3D getMaxActivityPosition() const;
    
    // Get voxel list (for iteration)
    const std::vector<VoxelActivity>& getVoxels() const { return voxels_; }
    
    // Clear all data
    void clear();
    
    // Voxel conversion utilities
    int positionToVoxelX(double x_mm) const;
    int positionToVoxelY(double y_mm) const;
    int positionToVoxelZ(double z_mm) const;
    
    Vector3D voxelToPosition(int x, int y, int z) const;
    bool isInGrid(int x, int y, int z) const;
    
private:
    int nx_, ny_, nz_;                  // Grid dimensions
    double voxelSize_;                   // mm
    std::vector<VoxelActivity> voxels_;
    
    // Find or create voxel
    VoxelActivity* getOrCreateVoxel(int x, int y, int z);
    
    // Tissue-dependent washout rates (1/s)
    double getWashoutRate(const std::string& tissue, 
                         const std::string& isotope) const;
};

} // namespace PET_MC

#endif // ACTIVITY_H
