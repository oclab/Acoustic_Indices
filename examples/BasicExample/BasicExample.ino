/**
 * @file BasicExample.ino
 * @brief Basic example for computing 7 common acoustic indices on ESP32
 * 
 * This example demonstrates how to compute the most commonly used acoustic indices:
 * - Acoustic Complexity Index (ACI)
 * - Acoustic Diversity Index (ADI)
 * - Acoustic Evenness Index (AEI)
 * - Bioacoustic Index (BI)
 * - Normalized Difference Sound Index (NDSI)
 * - Spectral Entropy (SH)
 * - Temporal Entropy (TH)
 * 
 * Hardware: ESP32
 */

#include <AudioFile.h>
#include <AcousticIndices.h>

// Sample audio configuration
const uint32_t SAMPLE_RATE = 22050;  // 22.05 kHz sample rate
const size_t BUFFER_SIZE = 4096;     // Audio buffer size

// Audio buffer (simulated data)
int16_t audioBuffer[BUFFER_SIZE];

void setup() {
    Serial.begin(115200);
    while (!Serial) {
        delay(10);
    }
    
    Serial.println("\n=== Acoustic Indices Library - Basic Example ===");
    Serial.println("ESP32 Arduino Framework\n");
    
    // Generate synthetic audio data (sine wave + noise for demo)
    generateTestAudio(audioBuffer, BUFFER_SIZE, SAMPLE_RATE);
    
    // Create AudioFile object
    AudioFile audioFile(audioBuffer, BUFFER_SIZE, SAMPLE_RATE);
    
    Serial.println("Audio file created:");
    Serial.printf("  Sample Rate: %d Hz\n", audioFile.getSampleRate());
    Serial.printf("  Duration: %.2f seconds\n", audioFile.getDuration());
    Serial.printf("  Samples: %d\n\n", audioFile.getDataLength());
    
    // Compute all 7 acoustic indices
    Serial.println("Computing Acoustic Indices...\n");
    
    // 1. Acoustic Complexity Index (ACI)
    Serial.println("1. Computing Acoustic Complexity Index (ACI)...");
    SpectrogramResult* spectroACI = AcousticIndices::computeSpectrogram(&audioFile, 512, 512, false, false);
    if (spectroACI != nullptr) {
        size_t jBin = 5 * SAMPLE_RATE / 512; // 5 seconds in frames
        ACIResult* aciResult = AcousticIndices::computeACI(spectroACI, jBin);
        if (aciResult != nullptr) {
            Serial.printf("   ACI main_value: %.4f\n\n", aciResult->mainValue);
            delete aciResult;
        }
        delete spectroACI;
    }
    
    // 2. Acoustic Diversity Index (ADI)
    Serial.println("2. Computing Acoustic Diversity Index (ADI)...");
    float adiValue = AcousticIndices::computeADI(&audioFile, 10000.0f, -50.0f, 1000.0f);
    Serial.printf("   ADI main_value: %.4f\n\n", adiValue);
    
    // 3. Acoustic Evenness Index (AEI)
    Serial.println("3. Computing Acoustic Evenness Index (AEI)...");
    float aeiValue = AcousticIndices::computeAEI(&audioFile, 10000.0f, -50.0f, 1000.0f);
    Serial.printf("   AEI main_value: %.4f\n\n", aeiValue);
    
    // 4. Bioacoustic Index (BI)
    Serial.println("4. Computing Bioacoustic Index (BI)...");
    SpectrogramResult* spectroBI = AcousticIndices::computeSpectrogram(&audioFile, 512, 256, true, false);
    if (spectroBI != nullptr) {
        float biValue = AcousticIndices::computeBI(spectroBI, 2000.0f, 8000.0f);
        Serial.printf("   BI main_value: %.4f\n\n", biValue);
        delete spectroBI;
    }
    
    // 5. Normalized Difference Sound Index (NDSI)
    Serial.println("5. Computing Normalized Difference Sound Index (NDSI)...");
    float ndsiValue = AcousticIndices::computeNDSI(&audioFile, 1024, 1000.0f, 2000.0f, 2000.0f, 11000.0f);
    Serial.printf("   NDSI main_value: %.4f\n\n", ndsiValue);
    
    // 6. Spectral Entropy (SH)
    Serial.println("6. Computing Spectral Entropy (SH)...");
    SpectrogramResult* spectroSH = AcousticIndices::computeSpectrogram(&audioFile, 512, 256, true, false);
    if (spectroSH != nullptr) {
        float shValue = AcousticIndices::computeSH(spectroSH);
        Serial.printf("   SH main_value: %.4f\n\n", shValue);
        delete spectroSH;
    }
    
    // 7. Temporal Entropy (TH)
    Serial.println("7. Computing Temporal Entropy (TH)...");
    float thValue = AcousticIndices::computeTH(&audioFile, true);
    Serial.printf("   TH main_value: %.4f\n\n", thValue);
    
    Serial.println("=== All indices computed successfully! ===");
    
    // Print memory info
    Serial.printf("\nFree heap: %d bytes\n", ESP.getFreeHeap());
}

void loop() {
    // Nothing to do in loop
    delay(1000);
}

/**
 * @brief Generate test audio data (sine wave with harmonics)
 */
void generateTestAudio(int16_t* buffer, size_t length, uint32_t sampleRate) {
    const float frequency1 = 440.0f;  // A4 note
    const float frequency2 = 880.0f;  // A5 note
    const float frequency3 = 1760.0f; // A6 note
    
    for (size_t i = 0; i < length; i++) {
        float t = (float)i / sampleRate;
        
        // Create complex signal with multiple frequencies
        float signal = 0.3f * sin(2.0f * PI * frequency1 * t) +
                      0.2f * sin(2.0f * PI * frequency2 * t) +
                      0.1f * sin(2.0f * PI * frequency3 * t);
        
        // Add some noise
        signal += 0.05f * ((float)random(-1000, 1000) / 1000.0f);
        
        // Convert to int16
        buffer[i] = (int16_t)(signal * 32767.0f * 0.7f);
    }
}
