/**
 * @file AcousticIndices.h
 * @brief Compute acoustic indices for biodiversity assessment
 * @author Patrice Guyot, Alice Eldridge, Mika Peck
 * @version 1.0.0
 * 
 * This library computes various acoustic indices used in soundscape ecology
 * and ecoacoustics research for biodiversity assessment.
 * 
 * References:
 * - Pieretti et al. (2011) - Acoustic Complexity Index
 * - Villanueva-Rivera et al. (2011) - Acoustic Diversity/Evenness Index  
 * - Boelman et al. (2007) - Bioacoustic Index
 * - Kasten et al. (2012) - Normalized Difference Sound Index
 * - Sueur et al. (2008) - Spectral/Temporal Entropy
 */

#ifndef ACOUSTIC_INDICES_H
#define ACOUSTIC_INDICES_H

#include <Arduino.h>
#include "AudioFile.h"

/**
 * @struct SpectrogramResult
 * @brief Structure to hold spectrogram data
 */
struct SpectrogramResult {
    float** data;           ///< 2D array [frequency][time]
    size_t numFreqBins;     ///< Number of frequency bins
    size_t numTimeFrames;   ///< Number of time frames
    float* frequencies;     ///< Array of frequencies (Hz)
    
    ~SpectrogramResult() {
        if (data != nullptr) {
            for (size_t i = 0; i < numFreqBins; i++) {
                delete[] data[i];
            }
            delete[] data;
        }
        if (frequencies != nullptr) {
            delete[] frequencies;
        }
    }
};

/**
 * @struct ACIResult
 * @brief Acoustic Complexity Index result
 */
struct ACIResult {
    float mainValue;        ///< Overall ACI value
    float* temporalValues;  ///< ACI values per time segment
    size_t numValues;       ///< Number of temporal values
    
    ~ACIResult() {
        if (temporalValues != nullptr) {
            delete[] temporalValues;
        }
    }
};

/**
 * @struct WaveSNRResult
 * @brief Wave Signal-to-Noise Ratio results
 */
struct WaveSNRResult {
    float SNR;                      ///< Signal-to-Noise Ratio (dB)
    float acousticActivity;         ///< Fraction of active frames
    int countAcousticEvents;        ///< Number of acoustic events
    float averageDuration;          ///< Average duration of events (seconds)
};

/**
 * @class AcousticIndices
 * @brief Main class for computing acoustic indices
 */
class AcousticIndices {
public:
    /**
     * @brief Compute spectrogram from audio signal
     * @param audioFile Audio file object
     * @param windowLength FFT window length (samples)
     * @param windowHop Hop size (samples)
     * @param scaleAudio Use float signal (-1 to 1) if true
     * @param square Square magnitude if true
     * @return SpectrogramResult structure
     */
    static SpectrogramResult* computeSpectrogram(
        AudioFile* audioFile,
        size_t windowLength = 512,
        size_t windowHop = 256,
        bool scaleAudio = true,
        bool square = true
    );
    
    /**
     * @brief Compute Acoustic Complexity Index (ACI)
     * @param spectro Spectrogram result
     * @param jBin Temporal window size (in spectrogram frames)
     * @return ACIResult structure
     */
    static ACIResult* computeACI(SpectrogramResult* spectro, size_t jBin);
    
    /**
     * @brief Compute Bioacoustic Index (BI)
     * @param spectro Spectrogram result
     * @param minFreq Minimum frequency (Hz)
     * @param maxFreq Maximum frequency (Hz)
     * @return BI value
     */
    static float computeBI(
        SpectrogramResult* spectro,
        float minFreq = 2000.0f,
        float maxFreq = 8000.0f
    );
    
    /**
     * @brief Compute Spectral Entropy (SH)
     * @param spectro Spectrogram result
     * @return Spectral entropy value
     */
    static float computeSH(SpectrogramResult* spectro);
    
    /**
     * @brief Compute Temporal Entropy (TH)
     * @param audioFile Audio file object
     * @param useInteger Use integer signal if true
     * @return Temporal entropy value
     */
    static float computeTH(AudioFile* audioFile, bool useInteger = true);
    
    /**
     * @brief Compute Normalized Difference Sound Index (NDSI)
     * @param audioFile Audio file object
     * @param windowLength Window length for Welch's method
     * @param anthroMin Minimum anthrophony frequency (Hz)
     * @param anthroMax Maximum anthrophony frequency (Hz)
     * @param bioMin Minimum biophony frequency (Hz)
     * @param bioMax Maximum biophony frequency (Hz)
     * @return NDSI value
     */
    static float computeNDSI(
        AudioFile* audioFile,
        size_t windowLength = 1024,
        float anthroMin = 1000.0f,
        float anthroMax = 2000.0f,
        float bioMin = 2000.0f,
        float bioMax = 11000.0f
    );
    
    /**
     * @brief Compute Acoustic Diversity Index (ADI)
     * @param audioFile Audio file object
     * @param maxFreq Maximum frequency (Hz)
     * @param dbThreshold Minimum dB threshold
     * @param freqStep Frequency step (Hz)
     * @return ADI value
     */
    static float computeADI(
        AudioFile* audioFile,
        float maxFreq = 10000.0f,
        float dbThreshold = -50.0f,
        float freqStep = 1000.0f
    );
    
    /**
     * @brief Compute Acoustic Evenness Index (AEI)
     * @param audioFile Audio file object
     * @param maxFreq Maximum frequency (Hz)
     * @param dbThreshold Minimum dB threshold
     * @param freqStep Frequency step (Hz)
     * @return AEI value (Gini coefficient)
     */
    static float computeAEI(
        AudioFile* audioFile,
        float maxFreq = 10000.0f,
        float dbThreshold = -50.0f,
        float freqStep = 1000.0f
    );

private:
    /**
     * @brief Compute FFT (real to complex)
     * @param input Input real signal
     * @param n Signal length (must be power of 2)
     * @param outputReal Output real part
     * @param outputImag Output imaginary part
     */
    static void computeFFT(float* input, size_t n, float* outputReal, float* outputImag);
    
    /**
     * @brief Apply Hanning window
     * @param data Input/output data
     * @param length Window length
     */
    static void applyHanningWindow(float* data, size_t length);
    
    /**
     * @brief Apply Hamming window
     * @param data Input/output data
     * @param length Window length
     */
    static void applyHammingWindow(float* data, size_t length);
    
    /**
     * @brief Compute Hilbert envelope
     * @param signal Input signal
     * @param length Signal length
     * @param output Output envelope
     */
    static void computeHilbertEnvelope(float* signal, size_t length, float* output);
    
    /**
     * @brief Compute Gini coefficient
     * @param values Array of values
     * @param n Number of values
     * @return Gini coefficient
     */
    static float computeGini(float* values, size_t n);
    
    /**
     * @brief Sort array in-place
     * @param arr Array to sort
     * @param n Array length
     */
    static void sortArray(float* arr, size_t n);
    
    /**
     * @brief Compute power spectral density using Welch's method
     * @param signal Input signal
     * @param length Signal length
     * @param sampleRate Sample rate (Hz)
     * @param windowLength Window length
     * @param outputPSD Output PSD array
     * @param outputFreqs Output frequency array
     * @param numFreqs Number of frequency bins
     */
    static void welchPSD(
        float* signal,
        size_t length,
        float sampleRate,
        size_t windowLength,
        float* outputPSD,
        float* outputFreqs,
        size_t* numFreqs
    );
};

#endif // ACOUSTIC_INDICES_H
