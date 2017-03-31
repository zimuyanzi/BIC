/**
* 
*/

#ifndef __SPEECH_RECOGNITON__
#define __SPEECH_RECOGNITON__

#include "main.h"
#include "bussiness/mfcc.h"
#include "opencv2/opencv.hpp"

using namespace cv;

#define ABS(a, b) a > b ? a - b : b - a

typedef enum _ERROR_CODE
{
	RECOGNITION_CODE_SUCCESS       = 0,              /** 成功 */
	RECOGNITION_CODE_FAILED        = 1,              /** 失败 */

	RECOGNITION_CODE_CONTINUE      = 3,              /** 任务继续 */
	RECOGNITION_CODE_FINISH        = 4,              /** 任务结束 */
}ERROR_CODE;

typedef struct _SChangePoint
{
	time_t    tChangeTime;          /** 说话人改变点时间 */
	double    dValue;               /** 说话人改变点值 */
}SChangePoint;

typedef struct _SAudioData
{
	uchar *   pAudio;     /** 音频原始数据 */
	int32_t   nSize;      /** 音频数据大小 */
	time_t    tPts;       /** 音频数据对应时间 */
	time_t    tStart;     /** 音频数据开始时间 */
	time_t    tEnd;       /** 音频数据结束时间 */
	_SAudioData()
	{
		pAudio = NULL;
	}
	~_SAudioData()
	{
		if (pAudio)
		{
			delete pAudio;
			pAudio = NULL;
		}
	}
}SAudioData;

typedef struct _SSilenceInfo
{
	int32_t nStartIdx;     /** 静音区间开始下标 */
	int32_t nEndIdx;       /** 静音区间结束下标 */
	int32_t nStartPos;     /** 静音区间开始位置 */
	int32_t nEndPos;       /** 静音区间结束位置 */
	int32_t nTotalSize;    /** 静音区间总大小 */
	time_t  tStartTime;
	time_t  tEndTime;
	_SSilenceInfo()
	{
		nStartIdx = -1;
		nEndIdx = -1;
		nStartPos = -1;
		nEndPos = -1;
		nTotalSize = 0;
	}
}SSilenceInfo;

typedef struct _SFeatureTime
{
	time_t tTime;                      /** 当前特征对应的时间 */
	vector<double> vectMfccFeature;    /** mfcc特征 */
}SFeatureTime;

void testSpeech();

class CSpeechRecognition
{
public:
	CSpeechRecognition(int32_t nSampleRate, int32_t nByte, int32_t nChannels);
	CSpeechRecognition();
	~CSpeechRecognition();
private:
	void set();
	bool initAudioParam(int32_t nSampleRate, int32_t nByte, int32_t nChannels);
	bool initFFT();
	bool initMFCC();
	void release();

public:
	/**
	/*  获取说话人改变点
	/*	
	/*  @param pAudio：          [ IN ]  推入的原始音频数据
	/*  @param nLength:          [ IN ]  音频数据的长度
	/*  @param tPts:             [ IN ]  当前音频数据对应的时间
	/*  @param vectChgPoint:     [ OUT ] 话人改变点时间，有检出则不为空
	/*  
	/*  @return                  result: 0 - success, !0 - failed
	*/
	int32_t getSpeakerChgPoint(IN uchar *pAudio, IN int32_t nDataSize, IN time_t tPts, OUT vector<SChangePoint> &vectChgPoint);

	/**
	/* 初始化音频参数
	*/
	bool init(IN int32_t nSampleRate, IN int32_t nByte, IN int32_t nChannels);

private:
	/** 推音频数据 */
	int32_t pushAudioData(IN uchar *pAudio, IN int32_t nDataSize);
	/** 推入原始音频数据 */
	int32_t pushOriAudioData(IN uchar *pAudio, IN int32_t nDataSize, IN time_t tPts);
	/** 移除静音部分 */
	int32_t delSlienceData();
	/** 设置静音位置 */
	int32_t setSilencePos(int32_t nCurIdx, int32_t nCurPos, int32_t nCurSize, time_t tCurTime, SSilenceInfo &silenceInfo);
	/** 重置静音位置 */
	int32_t resetSilencePos(SSilenceInfo &silenceInfo);
	/** 删除检出的静音片段 */
	int32_t delSlienceData(SSilenceInfo &silenceInfo);
	/** 将去除静音的音频推入列表 */
	int32_t pushAudioToLst(bool bSilenceEnd, SSilenceInfo &silenceInfo);
	/** 判断是否是静音点 */
	bool checkSilencePoint(IN uchar *pAudio, IN int32_t nSize);
	/** 判断是否是静音结束 */
	bool isSilenceEndPos(uchar *pAudio, int32_t nCurIdx, int32_t nCurPos, int32_t nCurSize, time_t tCurTime, SSilenceInfo &silenceInfo);
	/** 计算说话人改变点 */
	bool calcChangePoint(OUT vector<SChangePoint> &vectChgPoint);
	/** 查找峰值 */
	void findPeaks(IN int nLength, IN double * pBIC, IN int nMinPeakDistance, IN int nThreshold, OUT vector<int> &vIndex);
	/** 提取MFCC特征 */
	bool extractMFCCFeature();
	/** 提取mfcc特征 */
	bool extractMFCCFeature(time_t tCurTime);

public:
	/** 读文件 */
	bool readFile(string strFilePath, vector<uchar> &vectBuff);
	/** 遍历当前路径下所有pcm文件 */
	bool getAllFilePathByRootPath(string strRootPath, map<time_t, string> &vectFilePath);
	/** 写文件 */
	void writeFile(int32_t nIdx, uchar *pAudio, int32_t nSize);
	/** 读取真实跳变点 */
	void readRealPoint();
	/** 统计结果 */
	void staticReault(vector<SChangePoint> &vectChgPoint);

public:
	/** pts转毫秒 */
	time_t ptsToMs(time_t tPts);
	/** 毫秒转pts */
	time_t msToPts(time_t tMs);
	
private:
	
	bool                        m_bInitFlag;                                         /** 初始化标志 */

	int32_t                     m_nByte;                                             /** 每个音频点占字节 */
	int32_t                     m_nSampleRate;                                       /** 音频采样率 */
	int32_t                     m_nChannels;                                         /** 音频通道数 */
	uchar *                     m_pAudio;                                            /** 音频数据 */
	const static int32_t        m_nAudioBufferSize = 128*1024;                       /** 音频buffer大小 */
	int32_t                     m_nReadPos;                                          /** 读取音频的位置 */
	int32_t                     m_nWritePos;                                         /** 写音频的位置 */

	list<SAudioData *>          m_lstOriAudioData;                                   /** 原始音频数据 */
	list<SAudioData *>          m_lstAudio;                                          /** 无静音区间的音频 */

	

	fft_plan                    m_pFFTPlan;                                          /** fft plan */
	fft_complex *               m_pFFTIn;                                            /** fft 输入数据指针 */
	fft_complex *               m_pFFTOut;                                           /** fft 输出数据指针 */
	CMFCC *                     m_pmfcc;                                             /** mfcc计算器 */

	list<SFeatureTime>          m_lstFeature;                                        /** 未处理的特征信息 */
	int32_t                     m_nComFeatureNum;                                    /** 对照特征数，用以获取2s内最大峰值参照 */

public:
	float                       m_fBICWindowDuration;                                /** BIC窗口对应特征时长,单位秒 */
	float                       m_fBICStep;                                          /** BIC窗口步长，单位秒 */
	int32_t                     m_nBICPenaltyFactor;                                 /** BIC惩罚因子 */
	int32_t                     m_nThreshold;                                        /** BIC阈值 */
	const static int32_t        m_nFFTPointNum = 1024;                               /** fft变换点数 */
	const static int32_t        m_nStep = 200;                                       /** 语音帧间隔 */

public:

	list<time_t>                m_lstRealPoint;                                      /** 真实的改变点 */
	list<SChangePoint>          m_lstMissPoint;
	list<SChangePoint>          m_lstHitPoint;
	list<SChangePoint>          m_lstErrorPoint;
};

#endif