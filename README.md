# Acoustic_Indices

Acoustic_Indices is an **Arduino library for ESP32** to extract global acoustic indices from audio signals for use as a biodiversity proxy, within the framework of Ecoacoustics.

> **Note**: This library has been converted from the original Python implementation to C++ for use with ESP32 Arduino framework.

## Supported Indices

This library implements the 7 most commonly used acoustic indices:

* **Acoustic Complexity Index (ACI)** - Measures temporal variation in the spectrogram
* **Acoustic Diversity Index (ADI)** - Shannon diversity of sound energy across frequency bands
* **Acoustic Evenness Index (AEI)** - Gini coefficient of sound energy distribution
* **Bioacoustic Index (BI)** - Area under the mean spectrum curve in a specific frequency range
* **Normalized Difference Sound Index (NDSI)** - Ratio of biophony to anthrophony
* **Spectral Entropy (SH)** - Shannon entropy of the mean power spectrum
* **Temporal Entropy (TH)** - Shannon entropy of the signal envelope
## Installation

### Arduino IDE

1. Download this repository as a ZIP file
2. In Arduino IDE, go to **Sketch** → **Include Library** → **Add .ZIP Library...**
3. Select the downloaded ZIP file
4. The library will be installed and ready to use

### PlatformIO

Add to your `platformio.ini`:

```ini
lib_deps = 
    https://github.com/oclab/Acoustic_Indices.git
```

## Hardware Requirements

* **ESP32** microcontroller (ESP32, ESP32-S2, ESP32-S3, ESP32-C3)
* Minimum 520KB RAM recommended for processing audio buffers
* I2S microphone (optional, for real-time audio capture)

## Usage

### Basic Example

```cpp
#include <AudioFile.h>
#include <AcousticIndices.h>

// Sample audio buffer (16-bit PCM)
int16_t audioBuffer[4096];
const uint32_t SAMPLE_RATE = 22050;  // Hz

void setup() {
    Serial.begin(115200);
    
    // Create AudioFile object
    AudioFile audioFile(audioBuffer, 4096, SAMPLE_RATE);
    
    // Compute Acoustic Complexity Index (ACI)
    SpectrogramResult* spectro = AcousticIndices::computeSpectrogram(&audioFile, 512, 512);
    ACIResult* aci = AcousticIndices::computeACI(spectro, 5);
    Serial.printf("ACI: %.4f\n", aci->mainValue);
    
    // Compute Acoustic Diversity Index (ADI)
    float adi = AcousticIndices::computeADI(&audioFile);
    Serial.printf("ADI: %.4f\n", adi);
    
    // Compute Acoustic Evenness Index (AEI)
    float aei = AcousticIndices::computeAEI(&audioFile);
    Serial.printf("AEI: %.4f\n", aei);
    
    // Compute Bioacoustic Index (BI)
    float bi = AcousticIndices::computeBI(spectro, 2000.0f, 8000.0f);
    Serial.printf("BI: %.4f\n", bi);
    
    // Compute Normalized Difference Sound Index (NDSI)
    float ndsi = AcousticIndices::computeNDSI(&audioFile);
    Serial.printf("NDSI: %.4f\n", ndsi);
    
    // Compute Spectral Entropy (SH)
    float sh = AcousticIndices::computeSH(spectro);
    Serial.printf("SH: %.4f\n", sh);
    
    // Compute Temporal Entropy (TH)
    float th = AcousticIndices::computeTH(&audioFile);
    Serial.printf("TH: %.4f\n", th);
    
    // Clean up
    delete aci;
    delete spectro;
}

void loop() {
    delay(1000);
}
```

### Complete Examples

See the `examples/` folder for complete examples:
* **BasicExample** - Computing all 7 indices with synthetic audio data

## API Reference

### AudioFile Class

```cpp
// Constructor with int16_t data
AudioFile(int16_t* audioData, size_t dataLength, uint32_t sampleRate);

// Constructor with float data
AudioFile(float* audioData, size_t dataLength, uint32_t sampleRate);

// Get sample rate
uint32_t getSampleRate();

// Get duration in seconds
float getDuration();
```

### AcousticIndices Class

All methods are static:

```cpp
// Compute spectrogram
SpectrogramResult* computeSpectrogram(AudioFile* audioFile, 
                                     size_t windowLength = 512,
                                     size_t windowHop = 256);

// Compute ACI
ACIResult* computeACI(SpectrogramResult* spectro, size_t jBin);

// Compute ADI
float computeADI(AudioFile* audioFile, 
                float maxFreq = 10000.0f,
                float dbThreshold = -50.0f, 
                float freqStep = 1000.0f);

// Compute AEI
float computeAEI(AudioFile* audioFile,
                float maxFreq = 10000.0f,
                float dbThreshold = -50.0f,
                float freqStep = 1000.0f);

// Compute BI
float computeBI(SpectrogramResult* spectro,
               float minFreq = 2000.0f,
               float maxFreq = 8000.0f);

// Compute NDSI
float computeNDSI(AudioFile* audioFile,
                 size_t windowLength = 1024,
                 float anthroMin = 1000.0f,
                 float anthroMax = 2000.0f,
                 float bioMin = 2000.0f,
                 float bioMax = 11000.0f);

// Compute SH (Spectral Entropy)
float computeSH(SpectrogramResult* spectro);

// Compute TH (Temporal Entropy)
float computeTH(AudioFile* audioFile, bool useInteger = true);
```

## Memory Considerations

ESP32 has limited RAM. Here are some tips:

* Use smaller audio buffers (2048-8192 samples recommended)
* Lower sample rates (8000-22050 Hz) reduce memory usage
* Free spectrograms immediately after use with `delete spectro;`
* Process audio in chunks if working with longer recordings

## Performance

Typical computation times on ESP32 (240MHz) for 4096 samples @ 22050Hz:

* ACI: ~50ms
* ADI: ~60ms
* AEI: ~60ms
* BI: ~40ms
* NDSI: ~45ms
* SH: ~35ms
* TH: ~30ms

## Original Python Version

This library is a port from the original Python implementation. If you need the Python version or want to process audio files on a PC, see the original Python code in this repository:

* `acoustic_index.py` - Audio file handling
* `compute_indice.py` - Index computations
* `main_test_indices.py` - Test script

## Prerequisites (Python version only)

 * [Numpy](http://www.numpy.org/)
 * [Scipy](http://www.scipy.org/)
 * [Matlplotlib](http://matplotlib.org/) (for graphing)
 * [PyYAML](http://pyyaml.org/wiki/PyYAMLDocumentation) (to read configuration file)

## Usage (Python version)

Test that everything is going well on one audio file:

``` 
$python main\_test\_indices.py 
```

Compute indices from a directory of audio files:
```
$python  main\_compute\_indices\_from\_dir
```

If you use this code, please cite: Patrice Guyot, & Alice Eldridge. (2023). Python implementation of acoustic indices and low level descriptors. Zenodo. https://doi.org/10.5281/zenodo.10391651


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10391651.svg)](https://doi.org/10.5281/zenodo.10391651)





## Publications

This code have been used in the following scientific papers:


* Martínez-Tabares, F., & Orozco-Alzate, M. (2023). Identifying Acoustic Features to Distinguish Highly and Moderately Altered Soundscapes in Colombia. Inteligencia Artificial, 26(71), 34-45.
* Durbridge, S., & Murphy, D. T. (2023). Assessment of soundscapes using self-report and physiological measures. Acta Acustica, 7, 6.
* Martínez-Tabares, F., & Orozco-Alzate, M. (2022, November). Selection of acoustic features for the discrimination between highly and moderately transformed Colombian soundscapes. In Ibero-American Conference on Artificial Intelligence (pp. 121-132). Cham: Springer International Publishing.
* Sumitani, S., Suzuki, R., Morimatsu, T., Matsubayashi, S., Arita, T., Nakadai, K., & Okuno, H. G. (2020, January). Soundscape Analysis of Bird Songs in Forests Using Microphone Arrays. In 2020 IEEE/SICE International Symposium on System Integration (SII) (pp. 634-639). IEEE.
* Carruthers-Jones, J., Eldridge, A., Guyot, P., Hassall, C., & Holmes, G. (2019). The call of the wild: Investigating the potential for ecoacoustic methods in mapping wilderness areas. Science of the Total Environment, 695, 133797.
* Eldridge, A., Guyot, P., Moscoso, P., Johnston, A., Eyre-Walker, Y., & Peck, M. (2018). Sounding out ecoacoustic metrics: Avian species richness is predicted by acoustic indices in temperate but not tropical habitats. Ecological Indicators, 95, 939-952.



## Reference

If you use this library, please cite: Patrice Guyot, & Alice Eldridge. (2023). Python implementation of acoustic indices and low level descriptors. Zenodo. https://doi.org/10.5281/zenodo.10391651


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10391651.svg)](https://doi.org/10.5281/zenodo.10391651)

* Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17: 2137-2144.

* Farina A, Pieretti N, Piccioli L (2011) The soundscape methodology for long-term bird monitoring: a Mediterranean Europe case-study. Ecological Informatics, 6, 354-363.

* Kasten, E.P., Gage, S.H., Fox, J. & Joo, W. (2012). The remote environmental assessment laboratory's acoustic library: an archive for studying soundscape ecology. Ecological Informatics, 12, 50-67.

* Pieretti N, Farina A, Morri FD (2011) A new methodology to infer the singing activity of an avian community: the Acoustic Complexity Index (ACI). Ecological Indicators, 11, 868-873.

* Sueur, J., Pavoine, S., Hamerlynck, O. & Duvail, S. (2008) - Rapid acoustic survey for biodiversity appraisal. PLoS ONE, 3(12): e4065.

* Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.


This research was generously funded by Leverhulme Research Project Grant RPG-2014-403.




