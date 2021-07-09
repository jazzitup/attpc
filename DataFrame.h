#ifndef __DataFrame_h
#define __DataFrame_h

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>

struct Header {
    int8_t metaType;             // 1
    int8_t frameSize[3];         // 4
    int8_t frameItemSource;      // 5
    int8_t frameType[2];         // 7
    int8_t revision;             // 8
    int8_t headerSize[2];        // 10
    int8_t itemSize[2];          // 12
    int8_t nItems[4];            // 16
    unsigned char eventTime[6];  // 22
    unsigned char eventIdx[4];   // 26
    int8_t coboIdx;              // 27
    int8_t asadIdx;              // 28
    int8_t readOffset[2];        // 30
    int8_t status;               // 31
    int8_t hitPat_0[9];          // 40
    int8_t hitPat_1[9];          // 49
    int8_t hitPat_2[9];          // 58
    int8_t hitPat_3[9];          // 67
    int8_t multip_0[2];          // 69
    int8_t multip_1[2];          // 71
    int8_t multip_2[2];          // 73
    int8_t multip_3[2];          // 75
    int8_t windowOut[4];         // 79
    int8_t lastCell_0[2];        // 81
    int8_t lastCell_1[2];        // 83
    int8_t lastCell_2[2];        // 85
    int8_t lastCell_3[2];        // 87
    int8_t trash[41];            // 128
};

class DataFrame {
   private:
    std::ifstream grawFile_;
    Header header_;
    uint32_t item_;
    unsigned int eventIdx_ = -1;
    int ADC_[4][68][512] = {0};

   public:
    DataFrame();
    ~DataFrame();
    // graw
    void OpenGrawFile(const std::string&);
    void CloseGrawFile();
    //
    bool Decode();
    void ReadHeader();
    void ReadItem();
    unsigned int GetEventIdx();
    unsigned long long GetEventTime();
    int GetItems();
    // Frame Item Format : aacc cccc | cbbb bbbb | bb00 ssss | ssss ssss
    // Read in reverse order : ssss ssss | bb00 ssss | cbbb bbbb | aacc cccc
    // a : agetIdx, c : chanIdx, b : buckIdx, s : sample
    int GetAgetIdx();
    int GetChanIdx();
    int GetBuckIdx();
    int GetSample();
    bool IsFPNChannel(int chanIdx);
    int GetADC(int agetIdx, int chanIdx, int buckIdx);
    void Clear();
    void PrintHeader();

    bool IsEvenNoise(int agetIdx, int chanIdx);
};

#endif