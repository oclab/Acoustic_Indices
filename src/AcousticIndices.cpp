/**
 * @file AcousticIndices.cpp
 * @brief Implementation of acoustic indices computations
 */

#include "AcousticIndices.h"
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

// Helper function: Check if number is power of 2
static bool isPowerOf2(size_t n) {
    return (n != 0) && ((n & (n - 1)) == 0);
}

// Helper function: Find next power of 2
static size_t nextPowerOf2(size_t n) {
    if (isPowerOf2(n)) return n;
    size_t power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

// Simple FFT implementation (Cooley-Tukey algorithm)
void AcousticIndices::computeFFT(float* input, size_t n, float* outputReal, float* outputImag) {
    if (!isPowerOf2(n)) {
        Serial.println("ERROR: FFT size must be power of 2");
        return;
    }
    
    // Bit-reversal permutation
    size_t j = 0;
    for (size_t i = 0; i < n; i++) {
        outputReal[i] = input[j];
        outputImag[i] = 0.0f;
        
        size_t m = n >> 1;
        while (m <= j && m > 0) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    // FFT computation
    for (size_t s = 1; s <= log2(n); s++) {
        size_t m = 1 << s;
        float wm_real = cos(-2.0f * M_PI / m);
        float wm_imag = sin(-2.0f * M_PI / m);
        
        for (size_t k = 0; k < n; k += m) {
            float w_real = 1.0f;
            float w_imag = 0.0f;
            
            for (size_t j = 0; j < m/2; j++) {
                size_t t_idx = k + j + m/2;
                size_t u_idx = k + j;
                
                float t_real = w_real * outputReal[t_idx] - w_imag * outputImag[t_idx];
                float t_imag = w_real * outputImag[t_idx] + w_imag * outputReal[t_idx];
                
                float u_real = outputReal[u_idx];
                float u_imag = outputImag[u_idx];
                
                outputReal[u_idx] = u_real + t_real;
                outputImag[u_idx] = u_imag + t_imag;
                outputReal[t_idx] = u_real - t_real;
                outputImag[t_idx] = u_imag - t_imag;
                
                float w_real_new = w_real * wm_real - w_imag * wm_imag;
                w_imag = w_real * wm_imag + w_imag * wm_real;
                w_real = w_real_new;
            }
        }
    }
}

void AcousticIndices::applyHanningWindow(float* data, size_t length) {
    for (size_t i = 0; i < length; i++) {
        float window = 0.5f * (1.0f - cos(2.0f * M_PI * i / (length - 1)));
        data[i] *= window;
    }
}

void AcousticIndices::applyHammingWindow(float* data, size_t length) {
    for (size_t i = 0; i < length; i++) {
        float window = 0.54f - 0.46f * cos(2.0f * M_PI * i / (length - 1));
        data[i] *= window;
    }
}

void AcousticIndices::sortArray(float* arr, size_t n) {
    // Simple bubble sort (for small arrays)
    for (size_t i = 0; i < n - 1; i++) {
        for (size_t j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                float temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

float AcousticIndices::computeGini(float* values, size_t n) {
    if (n == 0) return 0.0f;
    
    // Create a copy and sort
    float* sorted = new float[n];
    memcpy(sorted, values, n * sizeof(float));
    sortArray(sorted, n);
    
    // Compute Gini coefficient
    float sum = 0.0f;
    float weightedSum = 0.0f;
    
    for (size_t i = 0; i < n; i++) {
        sum += sorted[i];
        weightedSum += sorted[i] * (i + 1);
    }
    
    delete[] sorted;
    
    if (sum == 0.0f) return 0.0f;
    
    float G = 2.0f * weightedSum / sum - (n + 1);
    return G / n;
}

SpectrogramResult* AcousticIndices::computeSpectrogram(
    AudioFile* audioFile,
    size_t windowLength,
    size_t windowHop,
    bool scaleAudio,
    bool square
) {
    if (audioFile == nullptr) return nullptr;
    
    // Ensure windowLength is power of 2
    if (!isPowerOf2(windowLength)) {
        windowLength = nextPowerOf2(windowLength);
        Serial.printf("Window length adjusted to %d (power of 2)\n", windowLength);
    }
    
    float* signal = scaleAudio ? audioFile->getSignalFloat() : nullptr;
    int16_t* signalInt = !scaleAudio ? audioFile->getSignalInt() : nullptr;
    size_t signalLength = audioFile->getDataLength();
    
    // Calculate number of frames
    size_t numFrames = (signalLength - windowLength) / windowHop + 1;
    size_t halfWindowLength = windowLength / 2;
    
    // Allocate result
    SpectrogramResult* result = new SpectrogramResult();
    result->numFreqBins = halfWindowLength;
    result->numTimeFrames = numFrames;
    result->data = new float*[halfWindowLength];
    for (size_t i = 0; i < halfWindowLength; i++) {
        result->data[i] = new float[numFrames];
    }
    
    // Compute frequencies
    float nyquist = audioFile->getNyquist();
    result->frequencies = new float[halfWindowLength];
    for (size_t i = 0; i < halfWindowLength; i++) {
        result->frequencies[i] = (float)i * nyquist / halfWindowLength;
    }
    
    // Allocate buffers for FFT
    float* frame = new float[windowLength];
    float* fftReal = new float[windowLength];
    float* fftImag = new float[windowLength];
    
    // Process each frame
    for (size_t frameIdx = 0; frameIdx < numFrames; frameIdx++) {
        size_t startIdx = frameIdx * windowHop;
        
        // Extract frame
        for (size_t i = 0; i < windowLength; i++) {
            if (scaleAudio) {
                frame[i] = signal[startIdx + i];
            } else {
                frame[i] = (float)signalInt[startIdx + i];
            }
        }
        
        // Apply Hanning window
        applyHanningWindow(frame, windowLength);
        
        // Compute FFT
        computeFFT(frame, windowLength, fftReal, fftImag);
        
        // Compute magnitude and store
        for (size_t i = 0; i < halfWindowLength; i++) {
            float magnitude = sqrt(fftReal[i] * fftReal[i] + fftImag[i] * fftImag[i]);
            if (square) {
                result->data[i][frameIdx] = magnitude * magnitude;
            } else {
                result->data[i][frameIdx] = magnitude;
            }
        }
    }
    
    delete[] frame;
    delete[] fftReal;
    delete[] fftImag;
    
    return result;
}

ACIResult* AcousticIndices::computeACI(SpectrogramResult* spectro, size_t jBin) {
    if (spectro == nullptr) return nullptr;
    
    size_t numTimeFrames = spectro->numTimeFrames;
    size_t numFreqBins = spectro->numFreqBins;
    
    // Calculate number of segments
    size_t numSegments = (numTimeFrames - 10) / jBin;
    if (numSegments == 0) return nullptr;
    
    ACIResult* result = new ACIResult();
    result->numValues = numSegments;
    result->temporalValues = new float[numSegments];
    result->mainValue = 0.0f;
    
    for (size_t segIdx = 0; segIdx < numSegments; segIdx++) {
        size_t startTime = segIdx * jBin;
        float aciValue = 0.0f;
        
        for (size_t freqIdx = 0; freqIdx < numFreqBins; freqIdx++) {
            float sumDiff = 0.0f;
            float sumTotal = 0.0f;
            
            for (size_t t = startTime; t < startTime + jBin - 1; t++) {
                float diff = fabs(spectro->data[freqIdx][t + 1] - spectro->data[freqIdx][t]);
                sumDiff += diff;
                sumTotal += spectro->data[freqIdx][t];
            }
            
            if (sumTotal > 0.0f) {
                aciValue += sumDiff / sumTotal;
            }
        }
        
        result->temporalValues[segIdx] = aciValue;
        result->mainValue += aciValue;
    }
    
    return result;
}

float AcousticIndices::computeSH(SpectrogramResult* spectro) {
    if (spectro == nullptr) return 0.0f;
    
    size_t numFreqBins = spectro->numFreqBins;
    size_t numTimeFrames = spectro->numTimeFrames;
    
    // Compute mean spectrum
    float* meanSpectrum = new float[numFreqBins];
    float totalSum = 0.0f;
    
    for (size_t i = 0; i < numFreqBins; i++) {
        meanSpectrum[i] = 0.0f;
        for (size_t j = 0; j < numTimeFrames; j++) {
            meanSpectrum[i] += spectro->data[i][j];
        }
        totalSum += meanSpectrum[i];
    }
    
    // Normalize
    if (totalSum > 0.0f) {
        for (size_t i = 0; i < numFreqBins; i++) {
            meanSpectrum[i] /= totalSum;
        }
    }
    
    // Compute entropy
    float entropy = 0.0f;
    for (size_t i = 0; i < numFreqBins; i++) {
        if (meanSpectrum[i] > 0.0f) {
            entropy += meanSpectrum[i] * log2(meanSpectrum[i]);
        }
    }
    
    delete[] meanSpectrum;
    
    return -entropy / log2(numFreqBins);
}

float AcousticIndices::computeTH(AudioFile* audioFile, bool useInteger) {
    if (audioFile == nullptr) return 0.0f;
    
    size_t length = audioFile->getDataLength();
    float* signal = new float[length];
    
    // Copy signal
    if (useInteger) {
        int16_t* sigInt = audioFile->getSignalInt();
        for (size_t i = 0; i < length; i++) {
            signal[i] = (float)sigInt[i];
        }
    } else {
        float* sigFloat = audioFile->getSignalFloat();
        memcpy(signal, sigFloat, length * sizeof(float));
    }
    
    // Compute Hilbert envelope (simplified - using absolute value)
    float* envelope = new float[length];
    for (size_t i = 0; i < length; i++) {
        envelope[i] = fabs(signal[i]);
    }
    
    // Normalize envelope
    float sumEnv = 0.0f;
    for (size_t i = 0; i < length; i++) {
        sumEnv += envelope[i];
    }
    
    if (sumEnv > 0.0f) {
        for (size_t i = 0; i < length; i++) {
            envelope[i] /= sumEnv;
        }
    }
    
    // Compute entropy
    float entropy = 0.0f;
    for (size_t i = 0; i < length; i++) {
        if (envelope[i] > 0.0f) {
            entropy += envelope[i] * log2(envelope[i]);
        }
    }
    
    delete[] signal;
    delete[] envelope;
    
    return -entropy / log2(length);
}

float AcousticIndices::computeBI(
    SpectrogramResult* spectro,
    float minFreq,
    float maxFreq
) {
    if (spectro == nullptr) return 0.0f;
    
    // Find frequency bins
    size_t minBin = 0, maxBin = 0;
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        if (spectro->frequencies[i] >= minFreq && minBin == 0) {
            minBin = i;
        }
        if (spectro->frequencies[i] <= maxFreq) {
            maxBin = i;
        }
    }
    
    if (minBin > 0) minBin--; // Match R code
    
    // Find max value in spectrogram
    float maxVal = 0.0f;
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            if (spectro->data[i][j] > maxVal) {
                maxVal = spectro->data[i][j];
            }
        }
    }
    
    if (maxVal == 0.0f) return 0.0f;
    
    // Compute mean spectrum in dB
    float* meanSpectrumdB = new float[spectro->numFreqBins];
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        float sum = 0.0f;
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            float val = spectro->data[i][j] / maxVal;
            if (val > 0.0f) {
                sum += pow(10.0f, 20.0f * log10(val) / 10.0f);
            }
        }
        float meanVal = sum / spectro->numTimeFrames;
        meanSpectrumdB[i] = 10.0f * log10(meanVal);
    }
    
    // Extract segment and normalize
    float minVal = meanSpectrumdB[minBin];
    for (size_t i = minBin + 1; i < maxBin; i++) {
        if (meanSpectrumdB[i] < minVal) {
            minVal = meanSpectrumdB[i];
        }
    }
    
    // Compute area
    float area = 0.0f;
    float freqBandHz = spectro->frequencies[1] - spectro->frequencies[0];
    for (size_t i = minBin; i < maxBin; i++) {
        area += (meanSpectrumdB[i] - minVal) / freqBandHz;
    }
    
    delete[] meanSpectrumdB;
    
    return area;
}

float* AcousticIndices::computeZCR(
    AudioFile* audioFile,
    size_t windowLength,
    size_t windowHop,
    size_t* numValues
) {
    if (audioFile == nullptr) return nullptr;
    
    int16_t* signal = audioFile->getSignalInt();
    size_t length = audioFile->getDataLength();
    
    size_t numFrames = (length - windowLength) / windowHop + 1;
    float* zcr = new float[numFrames];
    
    for (size_t frameIdx = 0; frameIdx < numFrames; frameIdx++) {
        size_t startIdx = frameIdx * windowHop;
        int crossings = 0;
        
        for (size_t i = startIdx; i < startIdx + windowLength - 1; i++) {
            if ((signal[i] >= 0 && signal[i + 1] < 0) ||
                (signal[i] < 0 && signal[i + 1] >= 0)) {
                crossings++;
            }
        }
        
        zcr[frameIdx] = (float)crossings / windowLength;
    }
    
    if (numValues != nullptr) {
        *numValues = numFrames;
    }
    
    return zcr;
}

float* AcousticIndices::computeRMSEnergy(
    AudioFile* audioFile,
    size_t windowLength,
    size_t windowHop,
    bool useInteger,
    size_t* numValues
) {
    if (audioFile == nullptr) return nullptr;
    
    size_t length = audioFile->getDataLength();
    size_t numFrames = (length - windowLength) / windowHop + 1;
    float* rms = new float[numFrames];
    
    if (useInteger) {
        int16_t* signal = audioFile->getSignalInt();
        for (size_t frameIdx = 0; frameIdx < numFrames; frameIdx++) {
            size_t startIdx = frameIdx * windowHop;
            float sum = 0.0f;
            
            for (size_t i = startIdx; i < startIdx + windowLength; i++) {
                float val = (float)signal[i];
                sum += val * val;
            }
            
            rms[frameIdx] = sqrt(sum / windowLength);
        }
    } else {
        float* signal = audioFile->getSignalFloat();
        for (size_t frameIdx = 0; frameIdx < numFrames; frameIdx++) {
            size_t startIdx = frameIdx * windowHop;
            float sum = 0.0f;
            
            for (size_t i = startIdx; i < startIdx + windowLength; i++) {
                sum += signal[i] * signal[i];
            }
            
            rms[frameIdx] = sqrt(sum / windowLength);
        }
    }
    
    if (numValues != nullptr) {
        *numValues = numFrames;
    }
    
    return rms;
}

float* AcousticIndices::computeSpectralCentroid(
    SpectrogramResult* spectro,
    size_t* numValues
) {
    if (spectro == nullptr) return nullptr;
    
    size_t numFrames = spectro->numTimeFrames;
    float* centroid = new float[numFrames];
    
    for (size_t t = 0; t < numFrames; t++) {
        float sumWeighted = 0.0f;
        float sumMagnitude = 0.0f;
        
        for (size_t f = 0; f < spectro->numFreqBins; f++) {
            float magnitude = spectro->data[f][t];
            sumWeighted += magnitude * spectro->frequencies[f];
            sumMagnitude += magnitude;
        }
        
        centroid[t] = (sumMagnitude > 0.0f) ? (sumWeighted / sumMagnitude) : 0.0f;
    }
    
    if (numValues != nullptr) {
        *numValues = numFrames;
    }
    
    return centroid;
}

float AcousticIndices::computeADI(
    AudioFile* audioFile,
    float maxFreq,
    float dbThreshold,
    float freqStep
) {
    if (audioFile == nullptr) return 0.0f;
    if (freqStep == 0.0f || maxFreq == 0.0f) return 0.0f;  // Prevent division by zero
    
    // Compute spectrogram for ADI
    float freqBandHz = maxFreq / freqStep;
    size_t windowLength = (size_t)(audioFile->getSampleRate() / freqBandHz);
    
    SpectrogramResult* spectro = computeSpectrogram(audioFile, windowLength, windowLength, true, false);
    if (spectro == nullptr) return 0.0f;
    
    // Find max value
    float maxVal = 0.0f;
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            if (spectro->data[i][j] > maxVal) {
                maxVal = spectro->data[i][j];
            }
        }
    }
    
    if (maxVal == 0.0f) {
        delete spectro;
        return 0.0f;
    }
    
    // Calculate frequency bands
    size_t numBands = (size_t)(maxFreq / freqStep);
    float* values = new float[numBands];
    
    float actualFreqBandHz = spectro->frequencies[1] - spectro->frequencies[0];
    
    for (size_t band = 0; band < numBands; band++) {
        float minF = band * freqStep;
        float maxF = (band + 1) * freqStep;
        
        size_t minBin = (size_t)(minF / actualFreqBandHz);
        size_t maxBin = (size_t)(maxF / actualFreqBandHz);
        
        if (maxBin >= spectro->numFreqBins) {
            maxBin = spectro->numFreqBins - 1;
        }
        
        size_t countAboveThreshold = 0;
        size_t totalCount = 0;
        
        for (size_t i = minBin; i <= maxBin; i++) {
            for (size_t j = 0; j < spectro->numTimeFrames; j++) {
                float db = 20.0f * log10(spectro->data[i][j] / maxVal);
                if (db > dbThreshold) {
                    countAboveThreshold++;
                }
                totalCount++;
            }
        }
        
        values[band] = (totalCount > 0) ? ((float)countAboveThreshold / totalCount) : 0.0f;
    }
    
    // Remove zeros
    size_t nonZeroCount = 0;
    for (size_t i = 0; i < numBands; i++) {
        if (values[i] != 0.0f) {
            nonZeroCount++;
        }
    }
    
    if (nonZeroCount == 0) {
        delete[] values;
        delete spectro;
        return 0.0f;
    }
    
    float* nonZeroValues = new float[nonZeroCount];
    size_t idx = 0;
    for (size_t i = 0; i < numBands; i++) {
        if (values[i] != 0.0f) {
            nonZeroValues[idx++] = values[i];
        }
    }
    
    // Compute Shannon diversity
    float sum = 0.0f;
    for (size_t i = 0; i < nonZeroCount; i++) {
        sum += nonZeroValues[i];
    }
    
    float shannon = 0.0f;
    for (size_t i = 0; i < nonZeroCount; i++) {
        float p = nonZeroValues[i] / sum;
        shannon += -p * log(p);
    }
    
    delete[] values;
    delete[] nonZeroValues;
    delete spectro;
    
    return shannon;
}

float AcousticIndices::computeAEI(
    AudioFile* audioFile,
    float maxFreq,
    float dbThreshold,
    float freqStep
) {
    if (audioFile == nullptr) return 0.0f;
    if (freqStep == 0.0f || maxFreq == 0.0f) return 0.0f;  // Prevent division by zero
    
    // Compute spectrogram for AEI
    float freqBandHz = maxFreq / freqStep;
    size_t windowLength = (size_t)(audioFile->getSampleRate() / freqBandHz);
    
    SpectrogramResult* spectro = computeSpectrogram(audioFile, windowLength, windowLength, true, false);
    if (spectro == nullptr) return 0.0f;
    
    // Find max value
    float maxVal = 0.0f;
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            if (spectro->data[i][j] > maxVal) {
                maxVal = spectro->data[i][j];
            }
        }
    }
    
    if (maxVal == 0.0f) {
        delete spectro;
        return 0.0f;
    }
    
    // Calculate frequency bands
    size_t numBands = (size_t)(maxFreq / freqStep);
    float* values = new float[numBands];
    
    float actualFreqBandHz = spectro->frequencies[1] - spectro->frequencies[0];
    
    for (size_t band = 0; band < numBands; band++) {
        float minF = band * freqStep;
        float maxF = (band + 1) * freqStep;
        
        size_t minBin = (size_t)(minF / actualFreqBandHz);
        size_t maxBin = (size_t)(maxF / actualFreqBandHz);
        
        if (maxBin >= spectro->numFreqBins) {
            maxBin = spectro->numFreqBins - 1;
        }
        
        size_t countAboveThreshold = 0;
        size_t totalCount = 0;
        
        for (size_t i = minBin; i <= maxBin; i++) {
            for (size_t j = 0; j < spectro->numTimeFrames; j++) {
                float db = 20.0f * log10(spectro->data[i][j] / maxVal);
                if (db > dbThreshold) {
                    countAboveThreshold++;
                }
                totalCount++;
            }
        }
        
        values[band] = (totalCount > 0) ? ((float)countAboveThreshold / totalCount) : 0.0f;
    }
    
    // Compute Gini coefficient
    float gini = computeGini(values, numBands);
    
    delete[] values;
    delete spectro;
    
    return gini;
}

float AcousticIndices::computeNDSI(
    AudioFile* audioFile,
    size_t windowLength,
    float anthroMin,
    float anthroMax,
    float bioMin,
    float bioMax
) {
    if (audioFile == nullptr) return 0.0f;
    
    // This is a simplified version - full Welch's method would be more complex
    // For now, compute a simple PSD estimate
    
    SpectrogramResult* spectro = computeSpectrogram(audioFile, windowLength, windowLength/2, true, true);
    if (spectro == nullptr) return 0.0f;
    
    // Compute mean PSD
    float* meanPSD = new float[spectro->numFreqBins];
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        meanPSD[i] = 0.0f;
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            meanPSD[i] += spectro->data[i][j];
        }
        meanPSD[i] /= spectro->numTimeFrames;
    }
    
    // Find bins for anthrophony and biophony
    float anthroSum = 0.0f;
    float bioSum = 0.0f;
    
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        float freq = spectro->frequencies[i];
        if (freq >= anthroMin && freq <= anthroMax) {
            anthroSum += meanPSD[i];
        }
        if (freq >= bioMin && freq <= bioMax) {
            bioSum += meanPSD[i];
        }
    }
    
    delete[] meanPSD;
    delete spectro;
    
    if (anthroSum + bioSum == 0.0f) return 0.0f;
    
    return (bioSum - anthroSum) / (bioSum + anthroSum);
}

WaveSNRResult AcousticIndices::computeWaveSNR(
    AudioFile* audioFile,
    size_t frameLength,
    float minDB,
    size_t windowSmoothing,
    float activityThresholdDB,
    size_t histNumBins,
    float dbRange,
    float N
) {
    WaveSNRResult result;
    result.SNR = 0.0f;
    result.acousticActivity = 0.0f;
    result.countAcousticEvents = 0;
    result.averageDuration = 0.0f;
    
    if (audioFile == nullptr) return result;
    
    float* signal = audioFile->getSignalFloat();
    size_t length = audioFile->getDataLength();
    
    // Compute wave envelope
    size_t numFrames = (length - frameLength) / frameLength + 1;
    float* waveEnv = new float[numFrames];
    
    for (size_t i = 0; i < numFrames; i++) {
        size_t startIdx = i * frameLength;
        float maxVal = 0.0f;
        
        for (size_t j = 0; j < frameLength; j++) {
            float absVal = fabs(signal[startIdx + j]);
            if (absVal > maxVal) {
                maxVal = absVal;
            }
        }
        
        waveEnv[i] = 20.0f * log10(maxVal + 1e-10f);
    }
    
    // Find minimum
    float minVal = waveEnv[0];
    for (size_t i = 1; i < numFrames; i++) {
        if (waveEnv[i] < minVal) {
            minVal = waveEnv[i];
        }
    }
    
    if (minVal < minDB) minVal = minDB;
    
    // Compute histogram
    size_t* hist = new size_t[histNumBins];
    for (size_t i = 0; i < histNumBins; i++) {
        hist[i] = 0;
    }
    
    float binWidth = dbRange / histNumBins;
    for (size_t i = 0; i < numFrames; i++) {
        if (waveEnv[i] >= minVal && waveEnv[i] < minVal + dbRange) {
            size_t bin = (size_t)((waveEnv[i] - minVal) / binWidth);
            if (bin >= histNumBins) bin = histNumBins - 1;
            hist[bin]++;
        }
    }
    
    // Find modal intensity
    size_t modalBin = 0;
    size_t maxCount = hist[0];
    for (size_t i = 1; i < histNumBins; i++) {
        if (hist[i] > maxCount) {
            maxCount = hist[i];
            modalBin = i;
        }
    }
    
    float backgroundNoise = minVal + modalBin * binWidth;
    
    // Compute SNR
    float maxEnv = waveEnv[0];
    for (size_t i = 1; i < numFrames; i++) {
        if (waveEnv[i] > maxEnv) {
            maxEnv = waveEnv[i];
        }
    }
    
    result.SNR = maxEnv - backgroundNoise;
    
    // Compute acoustic activity
    size_t activeFrames = 0;
    for (size_t i = 0; i < numFrames; i++) {
        if (waveEnv[i] > backgroundNoise + activityThresholdDB) {
            activeFrames++;
        }
    }
    
    result.acousticActivity = (float)activeFrames / numFrames;
    
    delete[] waveEnv;
    delete[] hist;
    
    return result;
}

int AcousticIndices::computeNBPeaks(
    SpectrogramResult* spectro,
    float freqBand,
    bool normalization,
    float slopeLeft,
    float slopeRight
) {
    if (spectro == nullptr) return 0;
    
    // Compute mean spectrum
    float* meanSpec = new float[spectro->numFreqBins];
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        meanSpec[i] = 0.0f;
        for (size_t j = 0; j < spectro->numTimeFrames; j++) {
            meanSpec[i] += spectro->data[i][j];
        }
        meanSpec[i] /= spectro->numTimeFrames;
    }
    
    // Normalize if requested
    if (normalization) {
        float maxVal = meanSpec[0];
        for (size_t i = 1; i < spectro->numFreqBins; i++) {
            if (meanSpec[i] > maxVal) {
                maxVal = meanSpec[i];
            }
        }
        if (maxVal > 0.0f) {
            for (size_t i = 0; i < spectro->numFreqBins; i++) {
                meanSpec[i] /= maxVal;
            }
        }
    }
    
    // Find peaks with slope threshold
    bool* isPeak = new bool[spectro->numFreqBins];
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        isPeak[i] = false;
    }
    
    for (size_t i = 1; i < spectro->numFreqBins - 1; i++) {
        bool leftSlope = meanSpec[i] > meanSpec[i - 1] + slopeLeft;
        bool rightSlope = meanSpec[i] > meanSpec[i + 1] + slopeRight;
        if (leftSlope && rightSlope) {
            isPeak[i] = true;
        }
    }
    
    // Remove close peaks (within freqBand)
    float freqBandHz = spectro->frequencies[1] - spectro->frequencies[0];
    size_t nbBin = (size_t)(freqBand / freqBandHz);
    
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        if (!isPeak[i]) continue;
        
        // Check for nearby peaks
        for (size_t j = i + 1; j < i + nbBin && j < spectro->numFreqBins; j++) {
            if (isPeak[j]) {
                // Keep the higher peak
                if (meanSpec[j] > meanSpec[i]) {
                    isPeak[i] = false;
                } else {
                    isPeak[j] = false;
                }
            }
        }
    }
    
    // Count peaks
    int numPeaks = 0;
    for (size_t i = 0; i < spectro->numFreqBins; i++) {
        if (isPeak[i]) {
            numPeaks++;
        }
    }
    
    delete[] meanSpec;
    delete[] isPeak;
    
    return numPeaks;
}
