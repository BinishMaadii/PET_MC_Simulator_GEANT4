#ifndef DETECTOR_GEOMETRY_H
#define DETECTOR_GEOMETRY_H

#include "Geometry/Voxel.h"
#include "Transport/Particle.h"
#include <vector>
#include <memory>

namespace PET_MC {

class DetectorGeometry {
public:
    // Constructor
    DetectorGeometry(int nx, int ny, int nz, double voxel_size_mm);
    
    // Get voxel at position
    Voxel* getVoxelAtPosition(const Vector3D& pos);
    const Voxel* getVoxelAtPosition(const Vector3D& pos) const;
    
    // Get voxel by index
    Voxel* getVoxel(int x, int y, int z);
    
    // Dimensions
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    int getNz() const { return nz_; }
    
    // Size
    double getVoxelSize() const { return voxelSize_; }
    
    // Total volume
    double getTotalVolume() const;
    
    // List all voxels
    const std::vector<std::unique_ptr<Voxel>>& getVoxels() const { return voxels_; }
    
    // Reset all voxels
    void reset();
    
private:
    int nx_, ny_, nz_;
    double voxelSize_;
    std::vector<std::unique_ptr<Voxel>> voxels_;
    
    // Helper for linear indexing
    int getLinearIndex(int x, int y, int z) const {
        return x + y * nx_ + z * nx_ * ny_;
    }
    
    bool isValidIndex(int x, int y, int z) const;
};

} // namespace PET_MC

#endif // DETECTOR_GEOMETRY_H
