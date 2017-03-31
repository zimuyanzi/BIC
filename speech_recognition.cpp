#include "bussiness/speech_recognition.h"
#include "log4cplus/XtLogger.h"
#include "vdl/vdl_dir.h"
#include "vdl/vdl_file.h"

using namespace log4cplus;
using namespace log4cplus_common;


void testSpeech()
{
	CSpeechRecognition speechRec(48000, 2, 1);
	int32_t nMaxBuff = 1000;
	map<time_t, string> vectFilePath;
	map<time_t, string>::iterator itFilePath;
	time_t tUtcTime = 0;
	double dUtcTime = 0;
	speechRec.readRealPoint();
	speechRec.getAllFilePathByRootPath("E:\\homed\\homed_project\\homed_project_splictVoice\\audio\\cctv1", vectFilePath);
	for (itFilePath = vectFilePath.begin(); itFilePath != vectFilePath.end(); itFilePath ++)
	{
		/** 读文件 */
		vector<uchar> vectBuff;
		if (speechRec.readFile(itFilePath->second, vectBuff) == false)
		{
			return;
		}

		int32_t nPosition = 0;
		int32_t nDataSize = 0;
		while (nPosition < vectBuff.size())
		{
			vector<SChangePoint> vectChgPoint;
			nDataSize = std::min(nMaxBuff, (int32_t)vectBuff.size() - nPosition);
			speechRec.getSpeakerChgPoint(&vectBuff[nPosition], nDataSize, (time_t)dUtcTime, vectChgPoint);

			/** 统计结果 */
			speechRec.staticReault(vectChgPoint);

			nPosition += nDataSize;
			dUtcTime = dUtcTime + nMaxBuff * 1000.0 / (48000*2);
		}
	}

	/** 打印参数 */
	FILE *pf = fopen("result.txt", "w+");
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_nBICWindowDuration :" << speechRec.m_fBICWindowDuration);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_fBICStep :" << speechRec.m_fBICStep);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_nBICPenaltyFactor :" << speechRec.m_nBICPenaltyFactor);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_nThreshold :" << speechRec.m_nThreshold);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_nFFTPointNum :" << speechRec.m_nFFTPointNum);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] m_nStep :" << speechRec.m_nStep);
	fprintf(pf, "BIC窗口间隔 : %f\n", speechRec.m_fBICWindowDuration);
	fprintf(pf, "BIC步长 : %f\n", speechRec.m_fBICStep);
	fprintf(pf, "惩罚因子 : %lld\n", speechRec.m_nBICPenaltyFactor);
	fprintf(pf, "BIC阈值 : %d\n", speechRec.m_nThreshold);
	fprintf(pf, "m_nFFTPointNum : %d\n", speechRec.m_nFFTPointNum);
	fprintf(pf, "m_nStep : %d\n", speechRec.m_nStep);

	/** 打印统计信息 */
	int32_t nTotalNum = 0;
	nTotalNum += speechRec.m_lstErrorPoint.size();
	nTotalNum += speechRec.m_lstHitPoint.size();
	nTotalNum += speechRec.m_lstMissPoint.size();

	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] nTotalNum :" << nTotalNum);
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] success :" << speechRec.m_lstHitPoint.size());
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] error :" << speechRec.m_lstErrorPoint.size());
	LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] miss :" << speechRec.m_lstMissPoint.size());
	fprintf(pf, "nTotalNum : %d\n", nTotalNum);
	fprintf(pf, "success : %d\n", speechRec.m_lstHitPoint.size());
	fprintf(pf, "error : %d\n", speechRec.m_lstErrorPoint.size());
	fprintf(pf, "miss : %d\n", speechRec.m_lstMissPoint.size());

	list<SChangePoint>::iterator itPoint;
	itPoint = speechRec.m_lstHitPoint.begin();
	while (itPoint != speechRec.m_lstHitPoint.end())
	{
		LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] success "
			<< ", time :" << itPoint->tChangeTime
			<< ", value :" << itPoint->dValue);
		fprintf(pf, "success : %lld ms, value: %.2f\n", itPoint->tChangeTime, itPoint->dValue);
		itPoint ++;
	}

	itPoint = speechRec.m_lstErrorPoint.begin();
	while (itPoint != speechRec.m_lstErrorPoint.end())
	{
		LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] error  "
			<< ", time :" << itPoint->tChangeTime
			<< ", value :" << itPoint->dValue);
		fprintf(pf, "error : %lld ms, value: %.2f\n", itPoint->tChangeTime, itPoint->dValue);
		itPoint ++;
	}

	itPoint = speechRec.m_lstMissPoint.begin();
	while (itPoint != speechRec.m_lstMissPoint.end())
	{
		LOG4CPLUS_INFO(g_oServerLogger, "[testSpeech] miss  "
			<< ", time :" << itPoint->tChangeTime
			<< ", value :" << itPoint->dValue);
		fprintf(pf, "miss : %lld ms, value: %.2f\n", itPoint->tChangeTime, itPoint->dValue);
		itPoint ++;
	}

	fclose(pf);
}

CSpeechRecognition::CSpeechRecognition()
{
	set();

	if (initFFT() == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::CSpeechRecognition] initFFT failed");
		goto EXIT_END;
	}


	m_bInitFlag = false;

	return;

EXIT_END:
	
	release();
}

CSpeechRecognition::CSpeechRecognition(int32_t nSampleRate, int32_t nByte, int32_t nChannels)
{
	set();

	if (initAudioParam(nSampleRate, nByte, nChannels) == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::CSpeechRecognition] initAudioParam failed");
		goto EXIT_END;
	}

	if (initFFT() == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::CSpeechRecognition] initFFT failed");
		goto EXIT_END;
	}

	if (initMFCC() == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::CSpeechRecognition] initMFCC failed");
		goto EXIT_END;
	}

	m_bInitFlag = true;

	return;
EXIT_END:
	
	release();
}

CSpeechRecognition::~CSpeechRecognition()
{
	release();
}

void CSpeechRecognition::set()
{
	m_bInitFlag = false;
	m_pFFTIn = NULL;
	m_pFFTOut = NULL;
	m_pFFTPlan = NULL;
	m_pmfcc = NULL;
	m_nComFeatureNum = 0;
	m_fBICWindowDuration = 3;
	m_fBICStep = 0.05f;
	m_nBICPenaltyFactor = 3;
	m_nThreshold = -100;
}

bool CSpeechRecognition::initAudioParam(int32_t nSampleRate, int32_t nByte, int32_t nChannels)
{
	/** 初始化音频参数 */
	if (nSampleRate <= 0 || nByte <= 0 || nChannels <= 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::initAudioParam] param error"
			<< ", nSampleRate :" << nSampleRate
			<< ", nByte :" << nByte
			<< ", nChannels :" << nChannels);
		return false;
	}
	m_nSampleRate = nSampleRate;
	m_nByte = nByte;
	m_nChannels = nChannels;

	m_pAudio = new uchar[m_nAudioBufferSize];
	m_nReadPos = 0;
	m_nWritePos = 0;

	return true;
}

bool CSpeechRecognition::initFFT()
{
	/** 初始化fft多线程安全函数 */
	fft_make_planner_thread_safe();
	/** 初始化 fftw 相关 */
	m_pFFTIn = (fft_complex*) fft_malloc(sizeof(fft_complex) * (m_nFFTPointNum + 1));
	m_pFFTOut = (fft_complex*) fft_malloc(sizeof(fft_complex) * (m_nFFTPointNum + 1));
	m_pFFTPlan = fft_plan_dft_1d(m_nFFTPointNum, m_pFFTIn, m_pFFTOut, FFTW_FORWARD, FFTW_MEASURE);

	return true;
}

bool CSpeechRecognition::initMFCC()
{
	/** 初始化mfcc */
	m_pmfcc = new CMFCC(m_nFFTPointNum, m_nByte, m_nSampleRate);
	if (m_pmfcc == NULL)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::initMFCC] new CMFCC failed");
		return false;
	}
	int nRect = m_pmfcc->init(m_pFFTPlan, m_pFFTIn, m_pFFTOut);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::initMFCC] m_pmfcc->init failed");
		return false;
	}

	return true;
}

void CSpeechRecognition::release()
{
	/** 释放fftw相关内存 */
	if (m_pFFTPlan)
	{
		fft_destroy_plan(m_pFFTPlan);
		m_pFFTPlan = NULL;
	}

	if (m_pFFTIn)
	{
		fft_free(m_pFFTIn);
		m_pFFTIn = NULL;
	}

	if (m_pFFTOut)
	{
		fft_free(m_pFFTOut);
		m_pFFTOut = NULL;
	}

	if (m_pmfcc)
	{
		delete m_pmfcc;
		m_pmfcc = NULL;
	}
	m_bInitFlag = false;
}

bool CSpeechRecognition::init(IN int32_t nSampleRate, IN int32_t nByte, IN int32_t nChannels)
{
	m_bInitFlag = false;

	if (initAudioParam(nSampleRate, nByte, nChannels) == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::init] initAudioParam failed");
		goto EXIT_END;
	}

	if (initMFCC() == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::init] initMFCC failed");
		goto EXIT_END;
	}

	m_bInitFlag = true;

EXIT_END:
	release();

	return true;
}

int32_t CSpeechRecognition::getSpeakerChgPoint(IN uchar *pAudio, IN int32_t nDataSize, IN time_t tPts, OUT vector<SChangePoint> &vectChgPoint)
{
	int32_t nRect = RECOGNITION_CODE_SUCCESS;
	vector< SFeatureTime > vectMFCCFeature;

	/** 检验参数 */
	if (pAudio == NULL || nDataSize <= 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getSpeakerChgPoint] param error"
			<< ", nDataSize :" << nDataSize);
		return RECOGNITION_CODE_FAILED;
	}

	/** 判断音频参数是否初始化成功 */
	if (m_bInitFlag == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getSpeakerChgPoint] audio data not init");
		return RECOGNITION_CODE_FAILED;
	}
	
	/** 将原始音频数据推入 */
	nRect = pushOriAudioData(pAudio, nDataSize, tPts);
	if (nRect != RECOGNITION_CODE_SUCCESS)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getSpeakerChgPoint] pushOriAudioData failed");
		return RECOGNITION_CODE_FAILED;
	}

	/** 隔一定时间则去掉音频数据中的静音部分 */
	//if (tPts < 1740000)
	if (m_lstOriAudioData.size() < 500)
	{
		return RECOGNITION_CODE_SUCCESS;
	}

#if 0
	/** 保存音频数据 */
	char strName[128];
	sprintf(strName, "src_audio.pcm");
	void *pf = vdl_file_open(strName, VDL_FILE_WRITE_ONLY);
	list<SAudioData *>::iterator itAudioData;
	itAudioData = m_lstOriAudioData.begin();
	while (itAudioData != m_lstOriAudioData.end())
	{
		vdl_file_write(pf, (*itAudioData)->pAudio, (*itAudioData)->nSize);
		itAudioData ++;
	}
	vdl_file_close(pf);

#endif

	//while (delSlienceData() == RECOGNITION_CODE_CONTINUE);

	list<SAudioData *>::iterator itAudioData;
	itAudioData = m_lstOriAudioData.begin();
	while (itAudioData != m_lstOriAudioData.end())
	{
		m_lstAudio.push_back(*itAudioData);
		itAudioData  = m_lstOriAudioData.erase(itAudioData);
	}

#if 0
	/** 保存音频数据 */
	sprintf(strName, "remove_slience_audio.pcm");
	pf = vdl_file_open(strName, VDL_FILE_WRITE_ONLY);
	itAudioData = m_lstAudio.begin();
	while (itAudioData != m_lstAudio.end())
	{
		vdl_file_write(pf, (*itAudioData)->pAudio, (*itAudioData)->nSize);
		itAudioData ++;
	}
	vdl_file_close(pf);

#endif

	/** 提取mfcc特征 */
	if (extractMFCCFeature() == false)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getSpeakerChgPoint] extractMFCCFeature failed");
		return RECOGNITION_CODE_FAILED;
	}

	/** 计算说话人改变点 */
	bool bRect = calcChangePoint(vectChgPoint);
	if (!bRect)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getSpeakerChgPoint] calcChangePoint failed");
		return RECOGNITION_CODE_FAILED;
	}

	return RECOGNITION_CODE_SUCCESS;
}

bool CSpeechRecognition::calcChangePoint(OUT vector<SChangePoint> &vectChgPoint)
{
	int32_t nWidth = 0, nHeight = 0;
	vector<time_t> vectTime;
	list<SFeatureTime>::iterator itFeature;

	if (m_lstFeature.size() == 0)
	{
		/** 无特征则不处理 */
		return false;
	}

	itFeature = m_lstFeature.begin();
	nWidth = (int32_t)itFeature->vectMfccFeature.size();
	nHeight = (int32_t)m_lstFeature.size();

	/** 存mfcc特征的矩阵 */
	Mat matFeature = Mat::zeros(nHeight, nWidth, CV_32FC1);

	/** 将特征写入Mat */
	int32_t nRow = 0;
	itFeature = m_lstFeature.begin();
	while (itFeature != m_lstFeature.end())
	{
		float *pData = (float*)matFeature.data + nRow * nWidth;
		for (int i = 0; i < nWidth; i ++)
		{
			pData[i] = (float)itFeature->vectMfccFeature[i];
		}

		vectTime.push_back(itFeature->tTime);
		nRow ++;
		itFeature ++;
	}

	//35*128对应着0.1秒的采样数据，BIC的窗口固定为2秒，对应690*128的采样数据
	int32_t nWinFeatureNum = m_nSampleRate * m_fBICWindowDuration/ m_nStep;
	int32_t nStepFeatureNum = (int32_t)(m_nSampleRate * m_fBICStep) / m_nStep;

	if (matFeature.rows < nWinFeatureNum*3)
	{
		return RECOGNITION_CODE_SUCCESS;
	}

	int nLength = (int)floor(float(matFeature.rows - nWinFeatureNum*2)/nStepFeatureNum);
	//GLR对应4秒窗口长度的惩罚因子，通过实验获得
	int nNameda = m_nBICPenaltyFactor;
	//特征维数
	int nFeatureDim = nWidth;
	double *pBIC = new double[nLength];
	for (int i = 0;i < nLength;i++)
	{
		Rect rc(0, i*nStepFeatureNum, nFeatureDim, nWinFeatureNum*2);
		Rect rc1(0, i*nStepFeatureNum, nFeatureDim, nWinFeatureNum);
		Rect rc2(0, i*nStepFeatureNum + nWinFeatureNum, nFeatureDim, nWinFeatureNum);
		Mat matDataWin = matFeature(rc).clone();
		Mat matDataWin1 = matFeature(rc1).clone();
		Mat matDataWin2 = matFeature(rc2).clone();

		int ctype = std::max(CV_32F, matDataWin.depth());
		Mat mCovar,mMean;
		calcCovarMatrix( matDataWin, mCovar, mMean, CV_COVAR_NORMAL|CV_COVAR_ROWS, ctype );
		mCovar = mCovar*1.0 / (nWinFeatureNum*2);
		Mat mCovar1,mMean1;
		calcCovarMatrix( matDataWin1, mCovar1, mMean1, CV_COVAR_NORMAL|CV_COVAR_ROWS, ctype );
		Mat mCovar2,mMean2;
		calcCovarMatrix( matDataWin2, mCovar2, mMean2, CV_COVAR_NORMAL|CV_COVAR_ROWS, ctype );
		mCovar1 = mCovar1*1.0 / nWinFeatureNum;
		mCovar2 = mCovar2*1.0 / nWinFeatureNum;

		double tmp = determinant(mCovar);
		tmp = determinant(mCovar1);
		tmp = determinant(mCovar2);

		pBIC[i] = nWinFeatureNum*log10(determinant(mCovar)) 
			- nWinFeatureNum/2*log10(determinant(mCovar1)) 
			- nWinFeatureNum/2*log10(determinant(mCovar2)) 
			- 0.5*nNameda*(nFeatureDim + 0.5*nFeatureDim*(nFeatureDim+1))*log10(float(nWinFeatureNum*2));

	}

	vector<int> vIndex;
	findPeaks(nLength, pBIC, nWinFeatureNum / nStepFeatureNum, m_nThreshold, vIndex);

	/** 计算重叠区域开始、结束下标 */
	int32_t nStratIdx = m_nComFeatureNum, nEndIdx = 0;
	nEndIdx = nLength * nStepFeatureNum - nWinFeatureNum;

	/** 保存峰值 */
	for (int i = 0; i < vIndex.size(); i ++)
	{
		SChangePoint tmpChgPoint;

		/** 判断是否需要加入结果 */
		int32_t nTimeIdx = vIndex[i] * nStepFeatureNum + nWinFeatureNum;
		if (nTimeIdx < m_nComFeatureNum)
		{
			continue;
		}
		if (nTimeIdx > nEndIdx)
		{
			break;
		}

		tmpChgPoint.tChangeTime = vectTime[nTimeIdx];
		tmpChgPoint.dValue = pBIC[vIndex[i]];
		vectChgPoint.push_back(tmpChgPoint);

		LOG4CPLUS_INFO(g_oServerLogger, "[CSpeechRecognition::calcChangePoint] "
			<< " tChangeTime :" << tmpChgPoint.tChangeTime
			<< " dValue :" << tmpChgPoint.dValue);
	}

	/** 删除已经处理过的mfcc结果 */
	m_nComFeatureNum = nLength * nStepFeatureNum - 2 * nWinFeatureNum;
	nRow = 0;
	itFeature = m_lstFeature.begin();
	while (itFeature != m_lstFeature.end())
	{
		if (nRow < m_nComFeatureNum)
		{
			itFeature = m_lstFeature.erase(itFeature);
		}
		else
		{
			break;
		}
		nRow ++;
	}

	m_nComFeatureNum = nWinFeatureNum;

	if (pBIC)
	{
		delete pBIC;
		pBIC = NULL;
	}

	return true;
}

int32_t CSpeechRecognition::pushOriAudioData(IN uchar *pAudio, IN int32_t nDataSize, IN time_t tPts)
{
	if (pAudio == NULL || nDataSize <= 0 || tPts < 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::pushOriAudioData] param error"
			<< ", pAudio :" << (void *)pAudio
			<< ", nDataSize :" << nDataSize
			<< ", tPts :" << tPts);
		return RECOGNITION_CODE_FAILED;
	}

	SAudioData *pAudioData = new SAudioData;

	pAudioData->pAudio = new uchar[nDataSize];
	memcpy(pAudioData->pAudio, pAudio, nDataSize);
	pAudioData->nSize = nDataSize;
	pAudioData->tPts = tPts;

	m_lstOriAudioData.push_back(pAudioData);

	return RECOGNITION_CODE_SUCCESS;
}

int32_t CSpeechRecognition::setSilencePos(int32_t nCurIdx, int32_t nCurPos, int32_t nCurSize, time_t tCurTime, SSilenceInfo &silenceInfo)
{
	if (silenceInfo.nStartIdx == -1 && silenceInfo.nStartPos == -1)
	{
		silenceInfo.nStartIdx = nCurIdx;
		silenceInfo.nStartPos = nCurPos;
		silenceInfo.tStartTime = tCurTime;
		silenceInfo.nTotalSize = 0;
	}
	else
	{
		silenceInfo.nTotalSize += nCurSize;
	}

	silenceInfo.tEndTime = tCurTime;
	silenceInfo.nEndIdx = nCurIdx;
	silenceInfo.nEndPos = nCurPos + nCurSize;

	return RECOGNITION_CODE_SUCCESS;
}

int32_t CSpeechRecognition::resetSilencePos(SSilenceInfo &silenceInfo)
{
	silenceInfo.nStartIdx = -1;
	silenceInfo.nStartPos = -1;
	silenceInfo.nEndIdx = -1;
	silenceInfo.nEndPos = -1;
	silenceInfo.nTotalSize = 0;
	return RECOGNITION_CODE_SUCCESS;
}

bool CSpeechRecognition::isSilenceEndPos(uchar *pAudio, int32_t nCurIdx, int32_t nCurPos, int32_t nCurSize, time_t tCurTime, SSilenceInfo &silenceInfo)
{
	const int32_t nMinSilenceDuration = 30;//ms
	const int32_t nMinSilenceSize = m_nSampleRate * m_nByte * nMinSilenceDuration / 1000;

	if (checkSilencePoint(pAudio, nCurSize))
	{
		setSilencePos(nCurIdx, nCurPos, nCurSize, tCurTime, silenceInfo);
	}
	else
	{
		if (silenceInfo.nTotalSize > nMinSilenceSize)
		{
			return true;
		}
		else
		{
			resetSilencePos(silenceInfo);
		}
	}

	return false;
}

int32_t CSpeechRecognition::delSlienceData()
{
	list<SAudioData *>::iterator itAudioData;
	bool bHaveMute = false;
	SSilenceInfo silenInfo;
	int32_t nCurIdx = 0, nCurPos = 0;
	const int32_t nTotalPoint = 128;
	const int32_t nTotalSize = nTotalPoint * m_nByte;
	uchar *pAudio = NULL;
	int32_t nLastLeftSize = 0;
	bool bSilenceEnd = false;

	itAudioData = m_lstOriAudioData.begin();
	while (itAudioData != m_lstOriAudioData.end())
	{
		for (nCurPos = 0; nCurPos < (*itAudioData)->nSize; nCurPos += nTotalSize)
		{
			if (nCurPos + nTotalSize <= (*itAudioData)->nSize)
			{
				if (isSilenceEndPos((*itAudioData)->pAudio + nCurPos, nCurIdx, nCurPos, nTotalSize, (*itAudioData)->tPts, silenInfo))
				{
					bSilenceEnd = true;
					break;
				}
			}
			else
			{
				break;
			}
		}

		if (bSilenceEnd)
		{
			break;
		}

		nCurIdx ++;
		itAudioData ++;
	}

	/** 将结果推入音频 */
	return pushAudioToLst(bSilenceEnd, silenInfo);
}

int32_t CSpeechRecognition::pushAudioToLst(bool bSilenceEnd, SSilenceInfo &silenceInfo)
{
	list<SAudioData *>::iterator itAudioData;
	int32_t nCurIdx = 0;

	/** 还未检测到静音结束点 */
	if (bSilenceEnd == false)
	{
		/** 将静音前面的部分推入结果 */
		itAudioData = m_lstOriAudioData.begin();
		while (itAudioData != m_lstOriAudioData.end())
		{
			if (nCurIdx < silenceInfo.nStartIdx || silenceInfo.nStartIdx == -1)
			{
				m_lstAudio.push_back((*itAudioData));
				itAudioData = m_lstOriAudioData.erase(itAudioData);
				nCurIdx ++;
				continue;
			}
			else
			{
				break;
			}
		}

		return RECOGNITION_CODE_FINISH;
	}

	/** 检测到静音点 */
	delSlienceData(silenceInfo);

	return RECOGNITION_CODE_CONTINUE;
}

int32_t CSpeechRecognition::delSlienceData(SSilenceInfo &silenceInfo)
{
	list<SAudioData *>::iterator itAudioData;
	int32_t nCurIdx = 0;

	itAudioData = m_lstOriAudioData.begin();
	while (itAudioData != m_lstOriAudioData.end())
	{
		if (nCurIdx < silenceInfo.nStartIdx)
		{
			/** 直接推入静音区间 */
			m_lstAudio.push_back((*itAudioData));
			itAudioData = m_lstOriAudioData.erase(itAudioData);
		}
		else if (nCurIdx == silenceInfo.nStartIdx)
		{
			//LOG4CPLUS_INFO(g_oServerLogger, "[CSpeechRecognition::delSlienceData] silence point start, time :" << (*itAudioData)->tPts);

			/** 分割静音开始部分 */
			uchar *pAudio = new uchar[silenceInfo.nStartPos + 1];
			memcpy(pAudio, (*itAudioData)->pAudio, silenceInfo.nStartPos);
			if ((*itAudioData)->pAudio)
			{
				delete (*itAudioData)->pAudio;
				(*itAudioData)->pAudio = NULL;
			}
			(*itAudioData)->pAudio = pAudio;
			(*itAudioData)->nSize = silenceInfo.nStartPos;
			m_lstAudio.push_back((*itAudioData));
			itAudioData = m_lstOriAudioData.erase(itAudioData);

		}
		else if (nCurIdx > silenceInfo.nStartIdx && nCurIdx < silenceInfo.nEndIdx)
		{
			/** 静音区间直接删掉 */
			if ((*itAudioData))
			{
				delete (*itAudioData);
				(*itAudioData) = NULL;
			}
			itAudioData = m_lstOriAudioData.erase(itAudioData);
		}
		else if (nCurIdx == silenceInfo.nEndIdx)
		{
			//LOG4CPLUS_INFO(g_oServerLogger, "[CSpeechRecognition::delSlienceData] silence point end, time :" << (*itAudioData)->tPts);

			/** 分割静音结束部分 */
			int32_t nSize = (*itAudioData)->nSize - silenceInfo.nEndPos;
			uchar *pAudio = new uchar[nSize + 1];
			memcpy(pAudio, (*itAudioData)->pAudio + silenceInfo.nEndPos, nSize);
			if ((*itAudioData)->pAudio)
			{
				delete (*itAudioData)->pAudio;
				(*itAudioData)->pAudio = NULL;
			}
			(*itAudioData)->pAudio = pAudio;
			(*itAudioData)->nSize = nSize;
			m_lstAudio.push_back((*itAudioData));
			itAudioData = m_lstOriAudioData.erase(itAudioData);
		}
		else
		{
			break;
		}

		nCurIdx ++;
	}

	return RECOGNITION_CODE_SUCCESS;
}

bool CSpeechRecognition::checkSilencePoint(IN uchar *pAudio, IN int32_t nSize)
{
	const double dMuteAvg = 100;
	double dSum = 0;
	int32_t nPointNum = nSize / m_nByte;
	short *sfAudio = NULL;
	double *pdAudio = NULL;

	if (nSize % m_nByte != 0 || pAudio == NULL || nSize <= 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::checkSilencePoint] param error");
		return false;
	}

	switch(m_nByte)
	{
	case 2:
		sfAudio = (short *)pAudio;
		for (int i = 0; i < nPointNum; i ++)
		{
			dSum += abs(sfAudio[i]);
		}
		break;
	case 4:
	case 8:
		pdAudio = (double *)pAudio;
		for (int i = 0; i < nPointNum; i ++)
		{
			dSum += pdAudio[i];
		}
		break;
	default:
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::checkSilencePoint] m_nByte error :" << m_nByte);
		return false;
	}

	dSum = dSum / nPointNum;

	if (dSum < dMuteAvg)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool CSpeechRecognition::extractMFCCFeature(time_t tCurTime)
{
	int nSamplePoint = m_nFFTPointNum * m_nByte;
	int nStep = m_nStep * m_nByte;

	while (m_nWritePos - m_nReadPos >= nSamplePoint)
	{
		SFeatureTime tmpFeature;
		vector<double> vectTmp;
		tmpFeature.tTime = tCurTime;
		m_pmfcc->reset();
		m_pmfcc->extractMFCCFeature(m_pAudio + m_nReadPos, nSamplePoint, vectTmp);

		for (int i = vectTmp.size() - 12; i < vectTmp.size(); i ++ )
		{
			tmpFeature.vectMfccFeature.push_back(vectTmp[i]);
		}

		m_lstFeature.push_back(tmpFeature);
		m_nReadPos += nStep;
	}

	return true;
}

int32_t CSpeechRecognition::pushAudioData(IN uchar *pAudio, IN int32_t nDataSize)
{
	if (m_nWritePos + nDataSize >= m_nAudioBufferSize)
	{
		for (int i = 0; i < m_nWritePos - m_nReadPos; i ++)
		{
			m_pAudio[i] = m_pAudio[i + m_nReadPos];
		}

		m_nWritePos = m_nWritePos - m_nReadPos;
		m_nReadPos = 0;
	}

	memcpy(m_pAudio + m_nWritePos, pAudio, nDataSize);
	m_nWritePos += nDataSize;

	return 0;
}

bool CSpeechRecognition::extractMFCCFeature()
{
	int nSamplePoint = m_nFFTPointNum * m_nByte;
	int nStep = m_nStep * m_nByte;

	double dMax = 0;
	list<SAudioData *>::iterator itAudioData;

	if (m_lstAudio.size() == 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::extractMFCCFeature] m_lstAudio.size() == 0");
		return false;
	}

	itAudioData = m_lstAudio.begin();
	while (itAudioData != m_lstAudio.end())
	{
		/** 推数据 */
		pushAudioData((*itAudioData)->pAudio, (*itAudioData)->nSize);
		/** 先读数据 */
		if (m_nWritePos - m_nReadPos > nSamplePoint)
		{
			extractMFCCFeature((*itAudioData)->tPts);
		}

		if ((*itAudioData))
		{
			delete (*itAudioData);
			(*itAudioData) = NULL;
		}

		itAudioData = m_lstAudio.erase(itAudioData);
	}

	return true;
}

void CSpeechRecognition::findPeaks(IN int nLength, IN double * pBIC, IN int nMinPeakDistance, IN int nThreshold, OUT vector<int> &vIndex) 
{
	list<int> lstIdx;

	for (int i = 1;i < nLength-1;i++ )  //先找局部极大值点
	{
		if (pBIC[i] < nThreshold)
		{
			continue;
		}
		if (pBIC[i] > pBIC[i-1] && pBIC[i] > pBIC[i+1])
		{
			lstIdx.push_back(i);
		}
	}

	//对相邻两个极值点之间的距离比nMinPeakDistance小的予以删除（只删除极值较小的那个极值点）
	while(true) 
	{
		bool bDelete = false;
		list<int>::iterator iter1, iter2;

		for (iter1 = lstIdx.begin(); iter1 != lstIdx.end();)
		{
			double dValue = pBIC[*iter1];

			iter2 = iter1;
			iter2 ++;

			if (iter2 == lstIdx.end())
			{
				break;
			}

			if (abs(*iter2 - *iter1) < nMinPeakDistance)
			{
				bDelete = true;
				if (dValue < pBIC[*iter2])
				{
					iter1 = lstIdx.erase(iter1);
					continue;
				}
				else
				{
					iter2 = lstIdx.erase(iter2);
					continue;
				}
			}
			else
			{
				iter1 ++;
			}
		}
		if (!bDelete)
		{
			break;
		}
	}

	list<int>::iterator iter;
	for (iter = lstIdx.begin(); iter != lstIdx.end(); iter ++)
	{
		vIndex.push_back(*iter);
	}
}

bool CSpeechRecognition::readFile(string strFilePath, vector<uchar> &vectBuff)
{
	int32_t nRect = 0;
	int64_t nFielSize = 0, tModTime = 0;
	void *pFile = NULL;

	/** 文件是否存在 */
	nRect = vdl_file_exist(strFilePath.c_str());
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::readFile] vdl_file_exist failed"
			<< ", strFilePath :" << strFilePath);
		return false;
	}

	/** 文件大小 */
	nRect = vdl_file_get_info(strFilePath.c_str(), &nFielSize, &tModTime);
	if (nRect != 0 || nFielSize == 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::readFile] vdl_file_get_info failed"
			<< ", strFilePath :" << strFilePath);
		return false;
	}

	if (nFielSize >= 0xFFFFFFF)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::readFile] nFielSize too big"
			<< ", nFielSize :" << nFielSize);
		return false;
	}

	/** 打开文件 */
	pFile = vdl_file_open(strFilePath.c_str(), VDL_FILE_READ_ONLY);
	if (pFile == NULL)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::readFile] vdl_file_open failed"
			<< ", strFilePath :" << strFilePath);
		return false;
	}

	vectBuff.clear();
	vectBuff.resize(nFielSize, 0);

	/** 读取文件 */
	nRect = vdl_file_read(pFile, &vectBuff[0], (int32_t)nFielSize);
	if (nRect == -1)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::readFile] vdl_file_read failed"
			<< ", strFilePath :" << strFilePath);
		vectBuff.clear();
		vdl_file_close(pFile);
		return false;
	}

	vdl_file_close(pFile);

	return true;
}

bool CSpeechRecognition::getAllFilePathByRootPath(string strRootPath, map<time_t, string> &vectFilePath)
{
	char szFile[256];
	int32_t nSize = 256, nType = 0;
	bool bRect = false;
	int32_t nRect = 0;
	void * hDir = NULL;

	if (vdl_dir_exist(NULL, strRootPath.c_str()) != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getAllFilePathByRootPath] vdl_dir_exist not exist"
			<< ", strRootPath :" << strRootPath);
		return false;
	}

	hDir = vdl_dir_open( NULL, strRootPath.c_str(), NULL );
	if(!hDir)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CSpeechRecognition::getAllFilePathByRootPath] vdl_dir_open failed, str :" << strRootPath);
		return false;
	}

	nRect = vdl_dir_get_first( hDir, &nType, szFile, nSize);
	while(nRect == 0)
	{
		if(nType == DFTYPE_FILE && strstr(szFile, ".pcm"))
		{
			int64_t nSize = 0;
			time_t tModTime = 0;
			vdl_file_get_info(szFile, &nSize, &tModTime);
			vectFilePath[tModTime] = szFile;
			//vectFilePath.push_back(szFile);
		}	
		nRect = vdl_dir_get_next( hDir, &nType, szFile, nSize);
	}
	vdl_dir_close(hDir);

	return true;
}

time_t CSpeechRecognition::ptsToMs(time_t tPts)
{
	return tPts * 1000 / 90;
}

time_t CSpeechRecognition::msToPts(time_t tMs)
{
	return tMs * 90 / 1000;
}

void CSpeechRecognition::writeFile(int32_t nIdx, uchar *pAudio, int32_t nSize)
{
	char strName[128];
	sprintf(strName, "bbbb_%d.pcm", nIdx);
	void *pf = vdl_file_open(strName, VDL_FILE_WRITE_ONLY);
	vdl_file_write(pf, pAudio, nSize);

	vdl_file_close(pf);
}

void CSpeechRecognition::readRealPoint()
{
	string strFileName = "point.txt";

	FILE *pf = fopen(strFileName.c_str(), "r+");
	while (!feof(pf))
	{
		time_t tTime;
		fscanf(pf, "%lld", &tTime);
		tTime = tTime * 1000;
		m_lstRealPoint.push_back(tTime);
	}

}

void CSpeechRecognition::staticReault(vector<SChangePoint> &vectChgPoint)
{
	list<time_t>::iterator itPoint;
	const time_t tDiff = 3000;

	for (int i = 0; i < vectChgPoint.size(); i ++)
	{
		itPoint = m_lstRealPoint.begin();
		while (itPoint != m_lstRealPoint.end())
		{
			time_t tAbs = ABS((*itPoint), vectChgPoint[i].tChangeTime);
			if (tAbs >= tDiff)
			{
				if ((*itPoint) > vectChgPoint[i].tChangeTime)
				{
					m_lstErrorPoint.push_back(vectChgPoint[i]);
					break;
				}
				else
				{
					SChangePoint tmp;
					tmp.tChangeTime = *itPoint;
					tmp.dValue = 0;
					m_lstMissPoint.push_back(tmp);
					itPoint = m_lstRealPoint.erase(itPoint);
					continue;
				}
			}
			else
			{
				m_lstHitPoint.push_back(vectChgPoint[i]);
				itPoint = m_lstRealPoint.erase(itPoint);
				break;
			}
			itPoint ++;
		}
	}
}