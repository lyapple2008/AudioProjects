#include "WavReader.h"
#include "VAD.h"

int main(int argc, char *argv[])
{
    const OsInt8 *filename = "C:\\CloudMusic\\De-Esse.wav";
    WavReader reader;
    reader.Open(filename);
    OsInt32 sampleRate = reader.GetSampleRate();
    OsInt32 channels = reader.GetChannels();
    OsInt32 bitsPerSample = reader.GetBitsPerSample();
    OsInt32 dataLen = reader.GetDataSize();

    OsInt32 winSize = 20 * sampleRate / 1000;
    OsInt32 winNum = 6;
    OsInt32 overlap = 50;

    OsInt16 *inbuf = new OsInt16[winSize];
    assert(inbuf);

    LTSD vadDetector(winSize, winNum, overlap);



    if (bitsPerSample != 16) {
        printf("Only support 16 bits per sample!!!\n");
        return -1;
    }

    int leftBytes = 0;
    int frameCnt = 0;
    while (dataLen > 0) {
        int offset = winSize - leftBytes;
        for (int i = 0; i < leftBytes; i++) {
            inbuf[i] = inbuf[offset + i];
        }
        reader.Read(inbuf + leftBytes, offset);

        int bytes = vadDetector.readFrame(inbuf, winSize);
        dataLen -= bytes;
        leftBytes = winSize - bytes;

        double result = vadDetector.computeVADRate();
        printf("%d -- %f\n", frameCnt++, result);
    }

    delete inbuf;

    return 0;
}