/**
 * @file AudioFile.cpp
 * @brief Audio file processing implementation
 */

#include "AudioFile.h"
#include <cmath>

AudioFile::AudioFile(int16_t* audioData, size_t dataLength, uint32_t sampleRate) {
    initFromInt(audioData, dataLength, sampleRate);
}

AudioFile::AudioFile(float* audioData, size_t dataLength, uint32_t sampleRate) {
    initFromFloat(audioData, dataLength, sampleRate);
}

AudioFile::~AudioFile() {
    if (ownsIntData && signalInt != nullptr) {
        delete[] signalInt;
    }
    if (ownsFloatData && signalFloat != nullptr) {
        delete[] signalFloat;
    }
}

void AudioFile::initFromInt(int16_t* audioData, size_t len, uint32_t sr) {
    this->dataLength = len;
    this->sampleRate = sr;
    this->signalInt = audioData;
    this->ownsIntData = false;
    
    // Convert to float
    this->signalFloat = new float[len];
    this->ownsFloatData = true;
    
    for (size_t i = 0; i < len; i++) {
        signalFloat[i] = pcm2float(audioData[i]);
    }
}

void AudioFile::initFromFloat(float* audioData, size_t len, uint32_t sr) {
    this->dataLength = len;
    this->sampleRate = sr;
    this->signalFloat = audioData;
    this->ownsFloatData = false;
    
    // Convert to int16
    this->signalInt = new int16_t[len];
    this->ownsIntData = true;
    
    for (size_t i = 0; i < len; i++) {
        signalInt[i] = float2pcm(audioData[i]);
    }
}

float AudioFile::pcm2float(int16_t input) {
    // Convert int16 (-32768 to 32767) to float (-1.0 to 1.0)
    return (float)input / 32768.0f;
}

int16_t AudioFile::float2pcm(float input) {
    // Clamp to [-1.0, 1.0]
    if (input > 1.0f) input = 1.0f;
    if (input < -1.0f) input = -1.0f;
    
    // Convert float to int16
    return (int16_t)(input * 32767.0f);
}
