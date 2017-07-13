#include <fftw3/fftw3.h>
#include <Eigen/Dense>

class LTSD
{
public:
    LTSD(int winSize, int winNum, int overlap);
    ~LTSD();

    int readFrame(short *in, int inSize);
    double computeVADRate();
private:
    void initHammingWin();
    void updateMagMatrix();
    void computeLTSE();
    void computeAverNoiseMag();
    double computeLTSD();
private:
    int m_winSize;
    int m_winNum;
    int m_overlap;
    int m_updateIndex;
    int m_overlapSize;

    double *m_frameBuf;
    double *m_winFrameBuf;
    double *m_hammWin;

    // fft
    fftw_plan m_plan;
    fftw_complex *m_tempFFT;

    Eigen::MatrixXd m_mags;
    Eigen::VectorXd m_ltse;
    Eigen::VectorXd m_avMag;
};