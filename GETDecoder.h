#ifndef __GETDecoder_h
#define __GETDecoder_h

#include <fstream>
#include <iostream>
#include <string>
#include <queue>

struct Header {
    uint8_t metaType;         // 1
    uint8_t frameSize[3];     // 4
    uint8_t frameItemSource;  // 5
    uint8_t frameType[2];     // 7
    uint8_t revision;         // 8
    uint8_t headerSize[2];    // 10
    uint8_t itemSize[2];      // 12
    uint8_t nItems[4];        // 16
    uint8_t eventTime[6];     // 22
    uint8_t eventId[4];       // 26
    uint8_t coboId;           // 27
    uint8_t asadId;           // 28
    uint8_t readOffset[2];    // 30
    uint8_t status;           // 31
    uint8_t hitPat_0[9];      // 40
    uint8_t hitPat_1[9];      // 49
    uint8_t hitPat_2[9];      // 58
    uint8_t hitPat_3[9];      // 67
    uint8_t multip_0[2];      // 69
    uint8_t multip_1[2];      // 71
    uint8_t multip_2[2];      // 73
    uint8_t multip_3[2];      // 75
    uint8_t windowOut[4];     // 79
    uint8_t lastCell_0[2];    // 81
    uint8_t lastCell_1[2];    // 83
    uint8_t lastCell_2[2];    // 85
    uint8_t lastCell_3[2];    // 87
    uint8_t trash[41];        // 128
};
class GETDecoder {
   private:
    Header header_;
    uint32_t item_;
    std::ifstream grawFile_;
    int year_;
    int month_;
    int day_;
    int hour_;
    int minute_;
    float second_;
    std::queue<std::string> fileQueue_;
    int ADC_[4][68][512];
    int maxADC_[4][68];
    int minADC_[4][68];
    int timeToMax_[4][68];
    int timeToMin_[4][68];
    uint64_t prevTime_, diffTime_;

    void ReadHeader();
    void ReadItem();
    int GetItems();
    int GetAgetId();
    int GetChanId();
    int GetBuckId();
    int GetSample();

   public:
    GETDecoder();
    void Open(const std::string& fileName);
    void OpenFromList(const std::string& fileListName);
    //
    bool Run();
    int GetEventId();
    int GetADC(int agetId, int chanId, int buckId);
    int GetMaxADC(int agetId, int chanId);
    int GetMinADC(int agetId, int chanId);
    int GetTimeToMax(int agetId, int chanId);
    int GetTimeToMin(int agetId, int chanId);
    uint64_t GetEventTime();
    uint64_t GetDiffTime();
};
#endif