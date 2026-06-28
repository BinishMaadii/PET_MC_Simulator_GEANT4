#ifndef VOXEL_H
#define VOXEL_H

#include "Transport/Particle.h"
#include <vector>

namespace PET_MC {

class Voxel {
public:
    Voxel(int x, int y, int z, double size_mm);
    
    // Voxel identification
    int getX() const { return x_; }
    int getY() const { return y_; }
    int getZ() const { return z_; }
    
    // Boundaries
    Vector3D getCenter() const;
    Vector3D getMin() const;
    Vector3D getMax() const;
    
    // Check if point is inside
    bool contains(const Vector3D& pos) const;
    
    // Voxel size
    double getSize() const { return size_; }
    double getVolume() const { return size_ * size_ * size_; } // mm³
    
    // Scoring
    double getDoseDeposited() const { return doseDeposited_; }
    void addDose(double dE) { doseDeposited_ += dE; }
    
    void reset() { doseDeposited_ = 0; }
    
private:
    int x_, y_, z_;
    double size_;
    double doseDeposited_;
};

} // namespace PET_MC

#endif // VOXEL_H
