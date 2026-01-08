/**
 * @file AudioFile.h
 * @brief Audio file processing class for acoustic indices
 * @author Patrice Guyot, Alice Eldridge, Mika Peck
 * @version 1.0.0
 * 
 * Acoustic_Indices is an Arduino library to extract global acoustic indices from 
 * audio signals for use as biodiversity proxies within the framework of Ecoacoustics.
 */

#ifndef AUDIOFILE_H
#define AUDIOFILE_H

#include <Arduino.h>

/**
 * @class AudioFile
 * @brief Class for handling audio data and computing acoustic indices
 */
class AudioFile {
public:
    /**
     * @brief Constructor for AudioFile
     * @param audioData Pointer to audio data buffer (int16_t samples)
     * @param dataLength Length of audio data
     * @param sampleRate Sample rate in Hz
     */
    AudioFile(int16_t* audioData, size_t dataLength, uint32_t sampleRate);
    
    /**
     * @brief Constructor for AudioFile with float data
     * @param audioData Pointer to audio data buffer (float samples between -1.0 and 1.0)
     * @param dataLength Length of audio data
     * @param sampleRate Sample rate in Hz
     */
    AudioFile(float* audioData, size_t dataLength, uint32_t sampleRate);
    
    /**
     * @brief Destructor
     */
    ~AudioFile();
    
    /**
     * @brief Get sample rate
     * @return Sample rate in Hz
     */
    uint32_t getSampleRate() const { return sampleRate; }
    
    /**
     * @brief Get Nyquist frequency
     * @return Nyquist frequency in Hz
     */
    float getNyquist() const { return sampleRate / 2.0f; }
    
    /**
     * @brief Get duration of audio
     * @return Duration in seconds
     */
    float getDuration() const { return (float)dataLength / sampleRate; }
    
    /**
     * @brief Get data length
     * @return Number of samples
     */
    size_t getDataLength() const { return dataLength; }
    
    /**
     * @brief Get integer signal data
     * @return Pointer to int16_t array
     */
    int16_t* getSignalInt() const { return signalInt; }
    
    /**
     * @brief Get float signal data
     * @return Pointer to float array
     */
    float* getSignalFloat() const { return signalFloat; }
    
    /**
     * @brief Convert PCM int16 to float (-1.0 to 1.0)
     * @param input Input int16_t value
     * @return Float value between -1.0 and 1.0
     */
    static float pcm2float(int16_t input);
    
    /**
     * @brief Convert float to PCM int16
     * @param input Input float value (-1.0 to 1.0)
     * @return int16_t value
     */
    static int16_t float2pcm(float input);

private:
    int16_t* signalInt;      ///< Integer signal data
    float* signalFloat;      ///< Float signal data (normalized -1.0 to 1.0)
    size_t dataLength;       ///< Length of audio data
    uint32_t sampleRate;     ///< Sample rate in Hz
    bool ownsIntData;        ///< Whether this object owns the int data
    bool ownsFloatData;      ///< Whether this object owns the float data
    
    /**
     * @brief Initialize from integer data
     */
    void initFromInt(int16_t* audioData, size_t dataLength, uint32_t sampleRate);
    
    /**
     * @brief Initialize from float data
     */
    void initFromFloat(float* audioData, size_t dataLength, uint32_t sampleRate);
};

#endif // AUDIOFILE_H
