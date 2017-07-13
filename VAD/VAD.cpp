#include "VAD.h"
#include <assert.h>

#define M_PI       3.14159265358979323846

LTSD::LTSD(int winSize, int winNum, int overlap)
{
    m_winSize = winSize;
    m_winNum = winNum;
    m_overlap = overlap;
    m_updateIndex = 0;
    m_overlapSize = m_winSize * m_overlap / 100;

    m_frameBuf = new double[m_winSize];
    assert(m_frameBuf);
    memset(m_frameBuf, 0, sizeof(double) * m_winSize);
    
    m_winFrameBuf = new double[m_winSize];
    assert(m_winFrameBuf);
    memset(m_winFrameBuf, 0, sizeof(double) * m_winSize);

    m_hammWin = new double[m_winSize];
    assert(m_hammWin);
    memset(m_hammWin, 0, sizeof(double) * m_winSize);
    initHammingWin();

    m_tempFFT = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_winSize);
    assert(m_tempFFT);
    m_plan = fftw_plan_dft_r2c_1d(m_winSize, m_winFrameBuf, m_tempFFT, FFTW_ESTIMATE);

    m_mags.setZero(m_winSize, 2 * m_winNum + 1);
    m_ltse.setZero(m_winSize);
}

LTSD::~LTSD()
{
    delete m_frameBuf;
    delete m_hammWin;
    fftw_free(m_tempFFT);
}

// return actually bytes read
int LTSD::readFrame(short *in, int inSize)
{
    // update frame buffer
    int shift = m_winSize - m_overlapSize;
    if (shift > inSize) {
        return 0;
    }
    int i = 0;
    for (; i < m_overlapSize; i++) {
        m_frameBuf[i] = m_frameBuf[shift + i];
    }
    for (; i < m_winSize; i++) {
        m_frameBuf[i] = *in;
        in++;
    }

    return shift;
}

double LTSD::computeVADRate()
{
    updateMagMatrix();
    computeLTSE();
    computeAverNoiseMag();
    double result = computeLTSD();

    return result;
}

void LTSD::initHammingWin()
{
    double alpha = 0.54;
    double beta = 1 - alpha;

    double coef = 1 / (m_winSize - 1);
    for (int i = 0; i < m_winSize; i++) {
        m_hammWin[i] = alpha - beta * cos(2 * M_PI * i) * coef;
    }
}

void LTSD::updateMagMatrix()
{
    for (int i = 0; i < m_winSize; i++) {
        m_winFrameBuf[i] = m_frameBuf[i] * m_hammWin[i];
    }

    fftw_execute(m_plan);
    for (int i = 0; i < m_winSize; i++) {
        double energy = pow(m_tempFFT[i][0], 2) + pow(m_tempFFT[i][1], 2);
        m_mags(i, m_updateIndex) = sqrt(energy);
    }

    m_updateIndex++;
    if (m_updateIndex == 2 * m_winNum + 1) {
        m_updateIndex = 0;
    }
}

void LTSD::computeLTSE()
{
    m_ltse.setZero(m_winSize);
    for (int i = 0; i < 2 * m_winNum + 1; i++) {
        for (int j = 0; j < m_winSize; j++) {
            if (m_mags(j, i) > m_ltse(j)) {
                m_ltse(j) = m_mags(j, i);
            }
        }
    }
}

void LTSD::computeAverNoiseMag()
{
    m_avMag.setZero(m_winSize);

    for (int i = 0; i < 2 * m_winNum + 1; i++) {
        for (int j = 0; j < m_winSize; j++) {
            m_avMag(j) += m_mags(j, i);
        }
    }

    for (int i = 0; i < m_winSize; i++) {
        m_avMag(i) = m_avMag(i) / (2*m_winNum + 1);
    }
}

double LTSD::computeLTSD()
{
    double ltsd = 0.0;
    
    for (int i = 0; i < m_winSize; i++) {
        ltsd += pow(m_ltse(i), 2) + pow(m_avMag(i), 2);
    }
    ltsd = 10 * log(ltsd / m_winSize);

    return ltsd;
}



