#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <map>

namespace PET_MC {

class Material {
public:
    // Material composition
    struct Composition {
        int Z;              // Atomic number
        int A;              // Mass number
        double fraction;    // Weight fraction
    };
    
    // Material properties
    struct Properties {
        double density;         // g/cm³
        double meanExcitationE; // MeV (mean ionization potential)
        std::vector<Composition> composition;
    };
    
    // Get material properties
    static Properties getProperties(const std::string& materialName);
    
    // Register a material
    static void registerMaterial(const std::string& name, 
                                const Properties& props);
    
    // Check if material exists
    static bool hasMaterial(const std::string& name);
    
    // List all registered materials
    static std::vector<std::string> getAvailableMaterials();
    
    // Common materials (pre-registered)
    static void initializeCommonMaterials();
    
private:
    static std::map<std::string, Properties> materials_;
};

} // namespace PET_MC

#endif // MATERIAL_H
