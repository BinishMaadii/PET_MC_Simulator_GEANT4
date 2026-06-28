#include "Physics/CrossSection.h"
#include "Physics/NuclearReaction.h"
#include "PET/Isotope.h"
#include "Transport/Material.h"
#include "Geometry/Voxel.h"
#include "Geometry/DetectorGeometry.h"
#include "PET/Activity.h"
#include <random>
#include <cmath>

namespace PET_MC {

// ============================================================================
// CrossSection Implementation
// ============================================================================

CrossSection g_crossSectionDB;

CrossSection::CrossSection() {
    initializeDefaultCrossSections();
}

void CrossSection::initializeDefaultCrossSections() {
    // Register proton + Carbon-12 reactions
    std::vector<CrossSectionPoint> p_c12_data = {
        {30, 0.0},   {50, 5.0},   {80, 15.0},  {110, 20.0},
        {150, 18.0}, {200, 15.0}
    };
    std::vector<ReactionProduct> p_c12_products = {
        {"C11", 6, 11, true, 1.0, 20.4*60}
    };
    registerReaction("p+C12", p_c12_data, p_c12_products);
    
    // p + O16 → O15 + n
    std::vector<CrossSectionPoint> p_o16_data = {
        {10, 0.0},   {20, 2.0},   {50, 8.0},   {100, 12.0},
        {150, 10.0}, {200, 8.0}
    };
    std::vector<ReactionProduct> p_o16_products = {
        {"O15", 8, 15, true, 1.0, 2.03*60}
    };
    registerReaction("p+O16", p_o16_data, p_o16_products);
}

void CrossSection::registerReaction(const std::string& name,
                                   const std::vector<CrossSectionPoint>& data,
                                   const std::vector<ReactionProduct>& products) {
    reactions_[name] = {data, products};
}

double CrossSection::getCrossSection(const std::string& reaction, double energy) const {
    auto it = reactions_.find(reaction);
    if (it == reactions_.end()) return 0.0;
    return interpolate(it->second.crossSectionData, energy);
}

std::vector<ReactionProduct> CrossSection::getProducts(const std::string& reaction) const {
    auto it = reactions_.find(reaction);
    if (it == reactions_.end()) return {};
    return it->second.products;
}

double CrossSection::getTotalCrossSection(const std::string& reaction, double energy) const {
    return getCrossSection(reaction, energy);
}

bool CrossSection::hasReaction(const std::string& reaction) const {
    return reactions_.find(reaction) != reactions_.end();
}

std::vector<std::string> CrossSection::getAvailableReactions() const {
    std::vector<std::string> result;
    for (const auto& pair : reactions_) {
        result.push_back(pair.first);
    }
    return result;
}

double CrossSection::interpolate(const std::vector<CrossSectionPoint>& data, double energy) const {
    if (data.empty()) return 0.0;
    if (energy < data[0].energy) return 0.0;
    if (energy > data.back().energy) return data.back().crossSection;
    
    for (size_t i = 0; i < data.size() - 1; i++) {
        if (energy >= data[i].energy && energy <= data[i+1].energy) {
            double x0 = data[i].energy, x1 = data[i+1].energy;
            double y0 = data[i].crossSection, y1 = data[i+1].crossSection;
            return y0 + (energy - x0) * (y1 - y0) / (x1 - x0);
        }
    }
    return 0.0;
}

// ============================================================================
// NuclearReaction Implementation
// ============================================================================

NuclearReaction::NuclearReaction(unsigned int seed)
    : rng_(seed), uniform_(0.0, 1.0) {}

ReactionEvent NuclearReaction::sampleReaction(const Particle& incident,
                                             const std::string& targetMaterial,
                                             double stepLength) const {
    ReactionEvent event;
    event.occurred = false;
    
    double prob = getReactionProbability(incident, targetMaterial, stepLength);
    if (uniform_(rng_) > prob) return event;
    
    // Reaction occurred!
    event.occurred = true;
    event.reactionEnergy = incident.getEnergy();
    event.reactionPosition = incident.getPosition();
    
    std::string reaction = "p+C12";
    auto products = g_crossSectionDB.getProducts(reaction);
    if (!products.empty()) {
        event.productName = products[0].name;
        event.isPositronEmitter = products[0].isPositronEmitter;
    }
    
    return event;
}

double NuclearReaction::getReactionProbability(const Particle& incident,
                                              const std::string& targetMaterial,
                                              double stepLength) const {
    double sigma = getCrossSection(incident, targetMaterial);
    if (sigma <= 1e-10) return 0.0;
    
    double number_density = 3.3e22 / 1e6;
    double prob = sigma * 1e-27 * number_density * stepLength;
    return std::min(prob, 1.0);
}

void NuclearReaction::generateSecondaryParticles(const Particle& incident,
                                                const ReactionProduct& product,
                                                ReactionEvent& event) const {
    double angle_theta = std::acos(2.0 * uniform_(rng_) - 1.0);
    double angle_phi = 2.0 * M_PI * uniform_(rng_);
    
    Vector3D secondary_dir(
        std::sin(angle_theta) * std::cos(angle_phi),
        std::sin(angle_theta) * std::sin(angle_phi),
        std::cos(angle_theta)
    );
    
    Particle neutron(ParticleType::NEUTRON, 5.0, incident.getPosition(), secondary_dir);
    event.secondaries.push_back(neutron);
}

NuclearReaction::TargetNucleus NuclearReaction::getTargetNucleus(const std::string& material) const {
    if (material == "water" || material == "tissue") {
        return {8, 16, 15.999};
    }
    return {6, 12, 12.0};
}

double NuclearReaction::getCrossSection(const Particle& incident,
                                       const std::string& targetMaterial) const {
    std::string reaction = "p+C12";
    return g_crossSectionDB.getCrossSection(reaction, incident.getEnergy());
}

// ============================================================================
// Isotope Implementation
// ============================================================================

Isotope::Isotope(const std::string& name) : name_(name) {
    if (name == "C11") {
        Z_ = 6; A_ = 11; halfLife_ = 20.4*60;
        isPET_ = true; Q_value_ = 1.982;
        positronBranching_ = 1.0; meanPositronEnergy_ = 0.96;
    } else if (name == "O15") {
        Z_ = 8; A_ = 15; halfLife_ = 2.03*60;
        isPET_ = true; Q_value_ = 1.732;
        positronBranching_ = 1.0; meanPositronEnergy_ = 1.72;
    } else if (name == "F18") {
        Z_ = 9; A_ = 18; halfLife_ = 110*60;
        isPET_ = true; Q_value_ = 0.635;
        positronBranching_ = 0.97; meanPositronEnergy_ = 0.64;
    } else if (name == "N13") {
        Z_ = 7; A_ = 13; halfLife_ = 10.0*60;
        isPET_ = true; Q_value_ = 1.944;
        positronBranching_ = 1.0; meanPositronEnergy_ = 1.19;
    } else {
        Z_ = 0; A_ = 0; halfLife_ = 1.0;
        isPET_ = false; Q_value_ = 0;
        positronBranching_ = 0; meanPositronEnergy_ = 0;
    }
    decayConstant_ = 0.693147 / halfLife_;
}

double Isotope::getActivityAtTime(double initialActivity, double timeElapsed) const {
    return initialActivity * std::exp(-decayConstant_ * timeElapsed);
}

double Isotope::getDecayedCount(double initialCount, double t1, double t2) const {
    return initialCount * (std::exp(-decayConstant_ * t1) - std::exp(-decayConstant_ * t2));
}

std::string Isotope::toString() const {
    return name_;
}

bool Isotope::registerIsotope(const std::string& name,
                             int Z, int A, double halfLife_s,
                             bool isPET, double Q_MeV,
                             double positronBranching) {
    // Simply validate and return true
    // Isotopes are created on-demand with Isotope(name) constructor
    return true;
}

Isotope* Isotope::getIsotope(const std::string& name) {
    // Create isotope on-demand
    return nullptr;  // Not used in current design
}

// ============================================================================
// Material Implementation
// ============================================================================

std::map<std::string, Material::Properties> Material::materials_;

Material::Properties Material::getProperties(const std::string& materialName) {
    auto it = materials_.find(materialName);
    if (it != materials_.end()) return it->second;
    
    // Default to water
    return {1.0, 0.0755, {}};
}

void Material::registerMaterial(const std::string& name, const Properties& props) {
    materials_[name] = props;
}

bool Material::hasMaterial(const std::string& name) {
    return materials_.find(name) != materials_.end();
}

std::vector<std::string> Material::getAvailableMaterials() {
    std::vector<std::string> result;
    for (const auto& pair : materials_) result.push_back(pair.first);
    return result;
}

void Material::initializeCommonMaterials() {
    registerMaterial("water", {1.0, 0.0755, {}});
    registerMaterial("bone", {1.85, 0.12, {}});
}

// ============================================================================
// Voxel Implementation
// ============================================================================

Voxel::Voxel(int x, int y, int z, double size_mm)
    : x_(x), y_(y), z_(z), size_(size_mm), doseDeposited_(0) {}

Vector3D Voxel::getCenter() const {
    return Vector3D((x_ + 0.5) * size_, (y_ + 0.5) * size_, (z_ + 0.5) * size_);
}

Vector3D Voxel::getMin() const {
    return Vector3D(x_ * size_, y_ * size_, z_ * size_);
}

Vector3D Voxel::getMax() const {
    return Vector3D((x_ + 1) * size_, (y_ + 1) * size_, (z_ + 1) * size_);
}

bool Voxel::contains(const Vector3D& pos) const {
    Vector3D min = getMin();
    Vector3D max = getMax();
    return pos.x >= min.x && pos.x < max.x &&
           pos.y >= min.y && pos.y < max.y &&
           pos.z >= min.z && pos.z < max.z;
}

// ============================================================================
// DetectorGeometry Implementation
// ============================================================================

DetectorGeometry::DetectorGeometry(int nx, int ny, int nz, double voxel_size_mm)
    : nx_(nx), ny_(ny), nz_(nz), voxelSize_(voxel_size_mm) {
    for (int i = 0; i < nx * ny * nz; i++) {
        int x = i % nx;
        int y = (i / nx) % ny;
        int z = i / (nx * ny);
        voxels_.push_back(std::make_unique<Voxel>(x, y, z, voxel_size_mm));
    }
}

Voxel* DetectorGeometry::getVoxelAtPosition(const Vector3D& pos) {
    int x = static_cast<int>(pos.x / voxelSize_);
    int y = static_cast<int>(pos.y / voxelSize_);
    int z = static_cast<int>(pos.z / voxelSize_);
    return getVoxel(x, y, z);
}

Voxel* DetectorGeometry::getVoxel(int x, int y, int z) {
    if (!isValidIndex(x, y, z)) return nullptr;
    return voxels_[getLinearIndex(x, y, z)].get();
}

double DetectorGeometry::getTotalVolume() const {
    return nx_ * ny_ * nz_ * voxelSize_ * voxelSize_ * voxelSize_;
}

void DetectorGeometry::reset() {
    for (auto& voxel : voxels_) {
        voxel->reset();
    }
}

bool DetectorGeometry::isValidIndex(int x, int y, int z) const {
    return x >= 0 && x < nx_ && y >= 0 && y < ny_ && z >= 0 && z < nz_;
}

// ============================================================================
// Activity Implementation
// ============================================================================

Activity::Activity(int nx, int ny, int nz, double voxel_size_mm)
    : nx_(nx), ny_(ny), nz_(nz), voxelSize_(voxel_size_mm) {}

void Activity::recordIsotope(const std::string& isotopeName,
                            const Vector3D& position,
                            double energy) {
    int x = static_cast<int>(position.x / voxelSize_);
    int y = static_cast<int>(position.y / voxelSize_);
    int z = static_cast<int>(position.z / voxelSize_);
    
    VoxelActivity* voxel = getOrCreateVoxel(x, y, z);
    if (voxel) {
        voxel->isotopeActivity[isotopeName] += 1.0;
        voxel->totalActivity += 1.0;
    }
}

void Activity::recordDose(const Vector3D& position, double dE) {
    int x = static_cast<int>(position.x / voxelSize_);
    int y = static_cast<int>(position.y / voxelSize_);
    int z = static_cast<int>(position.z / voxelSize_);
    getOrCreateVoxel(x, y, z);
}

double Activity::getActivityAtVoxel(int x, int y, int z, const std::string& isotopeFilter) const {
    return 0.0;
}

double Activity::getDoseAtVoxel(int x, int y, int z) const {
    return 0.0;
}

void Activity::applyDecay(double deltaTime_seconds) {}
void Activity::applyWashout(double washoutFraction, const std::string& tissue) {}

double Activity::getTotalActivity(const std::string& isotopeFilter) const {
    return 0.0;
}

bool Activity::exportToFile(const std::string& filename) const {
    return false;
}

bool Activity::exportToROOT(const std::string& filename) const {
    return false;
}

double Activity::getMaxActivity() const {
    return 0.0;
}

Vector3D Activity::getMaxActivityPosition() const {
    return Vector3D(0, 0, 0);
}

void Activity::clear() {
    voxels_.clear();
}

int Activity::positionToVoxelX(double x_mm) const {
    return static_cast<int>(x_mm / voxelSize_);
}

int Activity::positionToVoxelY(double y_mm) const {
    return static_cast<int>(y_mm / voxelSize_);
}

int Activity::positionToVoxelZ(double z_mm) const {
    return static_cast<int>(z_mm / voxelSize_);
}

Vector3D Activity::voxelToPosition(int x, int y, int z) const {
    return Vector3D((x + 0.5) * voxelSize_, (y + 0.5) * voxelSize_, (z + 0.5) * voxelSize_);
}

bool Activity::isInGrid(int x, int y, int z) const {
    return true;
}

VoxelActivity* Activity::getOrCreateVoxel(int x, int y, int z) {
    for (auto& v : voxels_) {
        if (v.x == x && v.y == y && v.z == z) return &v;
    }
    voxels_.push_back({x, y, z, {}, 0.0});
    return &voxels_.back();
}

double Activity::getWashoutRate(const std::string& tissue, const std::string& isotope) const {
    return 0.0;
}

} // namespace PET_MC
