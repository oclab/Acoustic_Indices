/**
 * @file I2S_Microphone.ino
 * @brief Example using I2S microphone to capture real-time audio and compute acoustic indices
 * 
 * This example demonstrates real-time audio capture from an I2S microphone (e.g., INMP441)
 * and computes the 7 common acoustic indices periodically.
 * 
 * Hardware connections (INMP441 example):
 * - VDD to 3.3V
 * - GND to GND
 * - SD to GPIO32 (I2S Data)
 * - WS to GPIO25 (I2S Word Select / LRCK)
 * - SCK to GPIO33 (I2S Bit Clock / SCLK)
 * - L/R to GND (left channel) or 3.3V (right channel)
 * 
 * Hardware: ESP32
 */

#include <driver/i2s.h>
#include <AudioFile.h>
#include <AcousticIndices.h>

// I2S Configuration
#define I2S_NUM         I2S_NUM_0
#define I2S_SAMPLE_RATE 22050
#define I2S_BCLK_PIN    33
#define I2S_LRC_PIN     25
#define I2S_DOUT_PIN    -1  // Not used for input
#define I2S_DIN_PIN     32

// Audio buffer configuration
#define BUFFER_SIZE     8192  // Audio samples to capture
#define CAPTURE_TIME_MS 5000  // Capture duration in milliseconds

int16_t audioBuffer[BUFFER_SIZE];

void setup() {
    Serial.begin(115200);
    while (!Serial) {
        delay(10);
    }
    
    Serial.println("\n=== Acoustic Indices - I2S Microphone Example ===");
    Serial.println("ESP32 with I2S Microphone\n");
    
    // Initialize I2S
    if (!initI2S()) {
        Serial.println("ERROR: Failed to initialize I2S!");
        while (1) {
            delay(1000);
        }
    }
    
    Serial.println("I2S initialized successfully!");
    Serial.printf("Sample Rate: %d Hz\n", I2S_SAMPLE_RATE);
    Serial.printf("Buffer Size: %d samples\n", BUFFER_SIZE);
    Serial.printf("Capture Duration: %d ms\n\n", CAPTURE_TIME_MS);
    
    delay(1000);
}

void loop() {
    Serial.println("=====================================");
    Serial.println("Capturing audio from microphone...");
    
    // Capture audio from I2S microphone
    if (!captureAudio(audioBuffer, BUFFER_SIZE)) {
        Serial.println("ERROR: Failed to capture audio!");
        delay(5000);
        return;
    }
    
    Serial.println("Audio captured successfully!");
    
    // Create AudioFile object
    AudioFile audioFile(audioBuffer, BUFFER_SIZE, I2S_SAMPLE_RATE);
    
    Serial.printf("Duration: %.2f seconds\n\n", audioFile.getDuration());
    
    // Compute all 7 acoustic indices
    Serial.println("Computing Acoustic Indices...\n");
    
    unsigned long startTime = millis();
    
    // 1. Acoustic Complexity Index (ACI)
    Serial.print("1. ACI... ");
    SpectrogramResult* spectroACI = AcousticIndices::computeSpectrogram(&audioFile, 512, 512, false, false);
    if (spectroACI != nullptr) {
        size_t jBin = 5 * I2S_SAMPLE_RATE / 512; // 5 seconds
        ACIResult* aciResult = AcousticIndices::computeACI(spectroACI, jBin);
        if (aciResult != nullptr) {
            Serial.printf("%.4f\n", aciResult->mainValue);
            delete aciResult;
        }
        delete spectroACI;
    }
    
    // 2. Acoustic Diversity Index (ADI)
    Serial.print("2. ADI... ");
    float adiValue = AcousticIndices::computeADI(&audioFile, 10000.0f, -50.0f, 1000.0f);
    Serial.printf("%.4f\n", adiValue);
    
    // 3. Acoustic Evenness Index (AEI)
    Serial.print("3. AEI... ");
    float aeiValue = AcousticIndices::computeAEI(&audioFile, 10000.0f, -50.0f, 1000.0f);
    Serial.printf("%.4f\n", aeiValue);
    
    // 4. Bioacoustic Index (BI)
    Serial.print("4. BI... ");
    SpectrogramResult* spectroBI = AcousticIndices::computeSpectrogram(&audioFile, 512, 256, true, false);
    if (spectroBI != nullptr) {
        float biValue = AcousticIndices::computeBI(spectroBI, 2000.0f, 8000.0f);
        Serial.printf("%.4f\n", biValue);
        delete spectroBI;
    }
    
    // 5. Normalized Difference Sound Index (NDSI)
    Serial.print("5. NDSI... ");
    float ndsiValue = AcousticIndices::computeNDSI(&audioFile, 1024, 1000.0f, 2000.0f, 2000.0f, 11000.0f);
    Serial.printf("%.4f\n", ndsiValue);
    
    // 6. Spectral Entropy (SH)
    Serial.print("6. SH... ");
    SpectrogramResult* spectroSH = AcousticIndices::computeSpectrogram(&audioFile, 512, 256, true, false);
    if (spectroSH != nullptr) {
        float shValue = AcousticIndices::computeSH(spectroSH);
        Serial.printf("%.4f\n", shValue);
        delete spectroSH;
    }
    
    // 7. Temporal Entropy (TH)
    Serial.print("7. TH... ");
    float thValue = AcousticIndices::computeTH(&audioFile, true);
    Serial.printf("%.4f\n", thValue);
    
    unsigned long elapsed = millis() - startTime;
    Serial.printf("\nTotal computation time: %lu ms\n", elapsed);
    Serial.printf("Free heap: %d bytes\n\n", ESP.getFreeHeap());
    
    // Wait before next capture
    delay(CAPTURE_TIME_MS);
}

/**
 * @brief Initialize I2S interface
 */
bool initI2S() {
    i2s_config_t i2s_config = {
        .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_RX),
        .sample_rate = I2S_SAMPLE_RATE,
        .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format = I2S_CHANNEL_FMT_ONLY_LEFT,
        .communication_format = I2S_COMM_FORMAT_STAND_I2S,
        .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
        .dma_buf_count = 8,
        .dma_buf_len = 1024,
        .use_apll = false,
        .tx_desc_auto_clear = false,
        .fixed_mclk = 0
    };
    
    i2s_pin_config_t pin_config = {
        .bck_io_num = I2S_BCLK_PIN,
        .ws_io_num = I2S_LRC_PIN,
        .data_out_num = I2S_DOUT_PIN,
        .data_in_num = I2S_DIN_PIN
    };
    
    esp_err_t err;
    
    // Install and start I2S driver
    err = i2s_driver_install(I2S_NUM, &i2s_config, 0, NULL);
    if (err != ESP_OK) {
        Serial.printf("Failed to install I2S driver: %d\n", err);
        return false;
    }
    
    err = i2s_set_pin(I2S_NUM, &pin_config);
    if (err != ESP_OK) {
        Serial.printf("Failed to set I2S pins: %d\n", err);
        return false;
    }
    
    // Clear I2S buffer
    err = i2s_zero_dma_buffer(I2S_NUM);
    if (err != ESP_OK) {
        Serial.printf("Failed to clear I2S buffer: %d\n", err);
        return false;
    }
    
    return true;
}

/**
 * @brief Capture audio from I2S microphone
 */
bool captureAudio(int16_t* buffer, size_t bufferSize) {
    size_t bytesRead = 0;
    size_t totalBytesRead = 0;
    size_t bytesToRead = bufferSize * sizeof(int16_t);
    
    // Read audio data from I2S
    esp_err_t err = i2s_read(I2S_NUM, buffer, bytesToRead, &bytesRead, portMAX_DELAY);
    
    if (err != ESP_OK) {
        Serial.printf("I2S read error: %d\n", err);
        return false;
    }
    
    if (bytesRead != bytesToRead) {
        Serial.printf("Warning: Expected %d bytes, read %d bytes\n", bytesToRead, bytesRead);
    }
    
    // Some I2S microphones need sample shifting
    // Uncomment if your microphone needs it:
    // for (size_t i = 0; i < bufferSize; i++) {
    //     buffer[i] = buffer[i] >> 8;  // Shift right 8 bits
    // }
    
    return true;
}
