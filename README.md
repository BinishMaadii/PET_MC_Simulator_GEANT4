# PET Activity Prediction Monte Carlo Simulator

A professional C++ Monte Carlo simulator for predicting Positron Emission Tomography (PET) activity distribution in particle therapy, with GEANT4 physics engine and ROOT visualization.

## 🎯 Overview

This simulator calculates:
- **Proton energy loss** using Bethe-Bloch formula
- **Nuclear reaction cross-sections** for PET isotope production (¹¹C, ¹⁵O, ¹⁸F)
- **Radioactive decay** modeling for activity predictions
- **Voxel-based dose/activity scoring** in 3D geometry
- **Bragg peak distributions** for hadron therapy verification

## ✨ Features

### Phase 1: Core Simulator ✅
- Bethe-Bloch energy loss implementation
- Nuclear reaction cross-section tables
- PET isotope decay chains (C-11, O-15, F-18, N-13)
- 3D voxel-based scoring geometry
- Custom Monte Carlo transport engine

### Phase 2: GEANT4 + ROOT Integration ✅
- **GEANT4** physics engine for realistic cross-sections
- **ROOT** histograms for professional visualization
- Bragg peak curves with publication-quality plots
- Decay simulation analysis
- ROOT file export for advanced analysis

## 📋 Requirements

### System
- macOS (Intel or Apple Silicon)
- CMake 3.16+
- C++17 compiler (Apple clang or GCC)

### Dependencies

## 🔧 Installation

### Option A: Quick Setup (15 minutes) - ROOT Only

```bash
# 1. Install ROOT
brew install root

# 2. Clone this repository
git clone https://github.com/BinishMaadii/PET_MC_Simulator_GEANT4.git
cd PET_MC_Simulator_GEANT4

# 3. Build
mkdir build && cd build
cmake ..
make -j4

# 4. Run
./bin/PET_Simulator
```

### Option B: Full Setup (90 minutes) - GEANT4 + ROOT

#### 2a. Install GEANT4
```bash
# Create installation directory
mkdir -p ~/Developer/physics-libraries
cd ~/Developer/physics-libraries

# Clone and build GEANT4
git clone https://github.com/Geant4/geant4.git
cd geant4
git checkout geant4-11.2
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DGEANT4_INSTALL_DATA=ON \
      -DGEANT4_BUILD_MULTITHREADED=ON \
      -DCMAKE_INSTALL_PREFIX=~/Developer/physics-libraries/geant4-install \
      ..

make -j$(sysctl -n hw.ncpu)  # Takes 20-40 minutes
make install
```

#### 2b. Set Environment Variables
Add to `~/.zshrc`:
```bash
export GEANT4_INSTALL=$HOME/Developer/physics-libraries/geant4-install
export PATH=$GEANT4_INSTALL/bin:$PATH
export LD_LIBRARY_PATH=$GEANT4_INSTALL/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$GEANT4_INSTALL/lib:$DYLD_LIBRARY_PATH
```

Then reload:
```bash
source ~/.zshrc
geant4-config --version  # Verify
```

#### 2c. Install ROOT
```bash
brew install root
root --version  # Verify
```

#### 2d. Build Simulator
```bash
cd /path/to/PET_MC_Simulator_GEANT4
rm -rf build && mkdir build && cd build
cmake ..
make -j4
```

## 🚀 Usage

### Run Simulator
```bash
./build/bin/PET_Simulator
```

**Expected Output:**


### Visualize Results (ROOT)
```bash
root PET_Simulation_Results.root

# Inside ROOT:
root> hEnergyLoss->Draw("L")      # Energy loss vs distance
root> hBraggPeak->Draw("L")       # Bragg peak curve
root> hDecayC11->Draw("L")        # C-11 decay
root> hActivityMap->Draw("colz")  # 2D activity map
root> gCrossSection->Draw("AP")   # Cross-section graph

# Save plots
root> hBraggPeak->SaveAs("bragg_peak.pdf")
root> .quit
```

## 🔬 Physics Implementation

### Bethe-Bloch Formula
Calculates energy loss of charged particles through matter:
-dE/dx = K * z² * (Z/A) * (1/β²) * [ln(2mₑβ²γ²Tₘₐₓ/I²) - 2β² - δ]

### Nuclear Reactions
Proton-induced reactions for PET isotope production:
- p + ¹²C → ¹¹C + n (σ ≈ 20 mb @ 100 MeV)
- p + ¹⁶O → ¹⁵O + n (σ ≈ 12 mb @ 100 MeV)

### PET Isotopes
| Isotope | Half-life | Decay Mode | Use |
|---------|-----------|-----------|-----|
| ¹¹C | 20.4 min | β⁺ | Brain imaging |
| ¹⁵O | 2.03 min | β⁺ | Flow measurement |
| ¹⁸F | 110 min | β⁺ | FDG-PET imaging |
| ¹³N | 10.0 min | β⁺ | Ammonia imaging |

## 📊 Validation

### Expected Physics Results (150 MeV proton in water)
- **Range:** ~140-150 mm ✓
- **Bragg Peak:** Sharp increase near end of range ✓
- **C-11 Decay:** Exponential with t₁/₂ = 20.4 min ✓
- **Cross-sections:** 5-20 millibarn in 50-200 MeV ✓

## 🛠️ Troubleshooting

### GEANT4 Not Found
```bash
# Verify installation
geant4-config --version

# If error, check environment variables
echo $GEANT4_INSTALL

# If empty, reload ~/.zshrc
source ~/.zshrc
```

### ROOT Not Found
```bash
# Verify installation
root --version

# If error:
brew install root
```

### Build Errors
```bash
# Clean and rebuild
rm -rf build
cmake -S . -B build
cmake --build build -j4
```

## 📚 References

- **GEANT4 Documentation:** https://geant4.web.cern.ch/
- **ROOT Documentation:** https://root.cern/
- **Paganetti, H. (2012):** "Proton Beam Therapy" (textbook standard)
- **Linz, U. (2019):** "Ion Beam Therapy" (comprehensive coverage)

## 🎓 What This Demonstrates

✅ Monte Carlo physics simulation methodology  
✅ C++17 modern software development  
✅ CMake cross-platform builds  
✅ Medical physics domain knowledge  
✅ Scientific visualization (ROOT)  
✅ GEANT4 integration for realistic physics  
✅ Professional code documentation  

## 📄 License

MIT License - See LICENSE file for details

## 👤 Author

**Dr. Binish Batool**
- Email: binish.maadii@gmail.com
- LinkedIn: linkedin.com/in/binish-batool
- GitHub: BinishMaadii
- Portfolio: sites.google.com/view/binishbatoolphd/bio

## 🤝 Contributing

Contributions welcome! Areas for enhancement:
- [ ] Full GEANT4 tracking engine
- [ ] DICOM CT image import
- [ ] Multi-particle batch simulation
- [ ] Machine learning for fast prediction
- [ ] Advanced ROOT analysis tools
- [ ] Clinical treatment plan integration

---

**Last Updated:** 2026-06-28  
**Status:** Active Development - Phase 2 Complete ✅



