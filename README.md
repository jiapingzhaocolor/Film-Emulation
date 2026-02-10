
# OpenDRT Film Pipeline (OFX)

This builds an OpenFX plugin that:
1) Interprets the input using **Input Gamut** + **Input Transfer**
2) Converts into **DaVinci Wide Gamut + DaVinci Intermediate**
3) Applies optional stages:
   - Negative (luma LUT): Cthulhu / Lilith / Tsathoggua / Yig + blend + enable
   - Color Separation: Hydra / Oorn / Zhar + blend + enable
   - Print: Kodak + blend + enable
4) Outputs Rec.709 (power 2.4) baseline, or Kodak print blended on top.

## Build (local)

```bash
git clone <your repo>
cd <repo>
git clone --depth 1 https://github.com/ofxa/openfx external/openfx
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The built module is `OpenDRTFilmPipeline.ofx` in `build/`.

## Packaging

OpenFX hosts expect a `.ofx.bundle` folder. The GitHub Actions workflow creates that bundle as an artifact.
