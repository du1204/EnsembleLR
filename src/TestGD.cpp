/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TestGD.h"
#include "EvaluatorUtils.h"
#include "Ciphertext.h"
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>
#include <NTL/BasicThreadPool.h>
#include "CipherGD.h"
#include "GD.h"

long TestGD::suggestLogN(long lambda, long logQ) {
	long NBnd = ceil(logQ * (lambda + 110) / 3.6);
	double logNBnd = log2((double)NBnd);
	return (long)ceil(logNBnd);
}

void TestGD::testEncNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
		bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", train samples = " << sampleDimTrain << ", test samples = " << sampleDimTest << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	GD::normalizezData2(zDataTrain, zDataTest, factorDim, sampleDimTrain, sampleDimTest);

	timeutils.start("Encrypting zData...");
	cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
	timeutils.stop("zData encryption");

	timeutils.start("Encrypting wData and vData...");
	if(isInitZero) {
		cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
	} else {
		cipherGD.encWVDataAverage(encWData, encVData, encZData, cnum, sBits, bBits);
	}
	timeutils.stop("wData and vData encryption");

	if(isInitZero) {
		GD::initialWDataVDataZero(pwData, pvData, factorDim);
		GD::initialWDataVDataZero(twData, tvData, factorDim);
	} else {
		GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
		GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
	}

	double alpha0, alpha1, eta, gamma;
	double enccor, encauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long iter = 0; iter < numIter; ++iter) {
		cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		cout << "encWData.logq before: " << encWData[0].logq << endl;
		timeutils.start("Enc NLGD");
		cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
		timeutils.stop("Enc NLGD");
		cout << "encWData.logq after: " << encWData[0].logq << endl;

		cout << "----ENCRYPTED-----" << endl;
		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
		GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
		cout << "------------------" << endl;

		GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);
//		cout << "-------TRUE-------" << endl;
//		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
//		cout << "------------------" << endl;

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
	}
	cout << "----ENCRYPTED-----" << endl;
	cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
	GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
	cout << "------------------" << endl;

	cout << "-------TRUE-------" << endl;
	GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
	cout << "------------------" << endl;

//	GD::calculateMSE(twData, cwData, factorDim);
//	GD::calculateNMSE(twData, cwData, factorDim);

}

void TestGD::testPlainNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
		bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", train samples = " << sampleDimTrain << ", test samples = " << sampleDimTest << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	TimeUtils timeutils;

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	GD::normalizezData2(zDataTrain, zDataTest, factorDim, sampleDimTrain, sampleDimTest);

	if(isInitZero) {
		GD::initialWDataVDataZero(pwData, pvData, factorDim);
		GD::initialWDataVDataZero(twData, tvData, factorDim);
	} else {
		GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
		GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
	}

	double alpha0, alpha1, eta, gamma;
	double plaincor, plainauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long iter = 0; iter < numIter; ++iter) {
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
		GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
	}


//	GD::calculateMSE(twData, pwData, factorDim);
//	GD::calculateNMSE(twData, pwData, factorDim);

}

void TestGD::testEncNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);

	double enccor, encauc, truecor, trueauc, plainauc, plaincor;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0;

	for (long fnum = 0; fnum < fold; ++fnum) {
		cout << " !!! START " << fnum + 1 << " FOLD !!! " << endl;

		for (long i = 0; i < sampleDimTest; ++i) {
			zDataTest[i] = zData[fnum * sampleDimTest + i];
		}
		for (long j = 0; j < fnum; ++j) {
			for (long i = 0; i < sampleDimTest; ++i) {
				zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
			}
		}
		for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
			zDataTrain[i - sampleDimTest] = zData[i];
		}

		timeutils.start("Encrypting zData...");
		cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
		timeutils.stop("zData encryption");

		timeutils.start("Encrypting wData and vData...");
		if(isInitZero) {
			cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
		} else {
			cipherGD.encWVDataAverage(encWData, encVData, encZData, cnum, sBits, bBits);
		}
		timeutils.stop("wData and vData encryption");

		if(isInitZero) {
			GD::initialWDataVDataZero(pwData, pvData, factorDim);
			GD::initialWDataVDataZero(twData, tvData, factorDim);
		} else {
			GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
			GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

			cout << "encWData.logq before: " << encWData[0].logq << endl;
			timeutils.start("Enc NLGD");
			cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
			timeutils.stop("Enc NLGD");
			cout << "encWData.logq after: " << encWData[0].logq << endl;

//			cout << "----ENCRYPTED-----" << endl;
//			cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
//			GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest);
//			cout << "------------------" << endl;

			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

//			cout << "-------TRUE-------" << endl;
//			GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, correctness, auc);
//			cout << "------------------" << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}
		cout << "----ENCRYPTED-----" << endl;
		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
		GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
		cout << "------------------" << endl;
		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

		averenccor += enccor;
		averencauc += encauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

//		GD::calculateMSE(twData, cwData, factorDim);
//		GD::calculateNMSE(twData, cwData, factorDim);

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Average Encrypted correctness: " << averenccor /fold << "%" << endl;
	cout << "Average Encrypted AUC: " << averencauc /fold << endl;
	cout << "Average True correctness: " << avertruecor /fold << "%" << endl;
	cout << "Average True AUC: " << avertrueauc /fold << endl;
}

void TestGD::testPlainNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	TimeUtils timeutils;

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);

	double plaincor, plainauc, truecor, trueauc;
	double averplaincor = 0, averplainauc = 0, avertruecor = 0, avertrueauc = 0;

	for (long fnum = 0; fnum < fold; ++fnum) {
		cout << " !!! START " << fnum + 1 << " FOLD !!! " << endl;

		for (long i = 0; i < sampleDimTest; ++i) {
			zDataTest[i] = zData[fnum * sampleDimTest + i];
		}
		for (long j = 0; j < fnum; ++j) {
			for (long i = 0; i < sampleDimTest; ++i) {
				zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
			}
		}
		for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
			zDataTrain[i - sampleDimTest] = zData[i];
		}

		if(isInitZero) {
			GD::initialWDataVDataZero(pwData, pvData, factorDim);
			GD::initialWDataVDataZero(twData, tvData, factorDim);
		} else {
			GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
			GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		}
		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

		averplaincor += plaincor;
		averplainauc += plainauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Average Plain correctness: " << averplaincor / fold << "%" << endl;
	cout << "Average Plain AUC: " << averplainauc / fold << endl;
	cout << "Average True correctness: " << avertruecor / fold << "%" << endl;
	cout << "Average True AUC: " << avertrueauc / fold << endl;
}


void TestGD::testEncEnsembleNLGD(long fold, double** zData, long factorDim, long sampleDim,
	bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	//SHOULD USE POWER-OF-2 part!!!!!!
	double bound = 0.5;
	long part = 4;
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim <<", fold = " << fold << ", part = " << part << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", initial bound = " << bound << endl;

	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch =  1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long sampleDimEnsemble = (1 << sdimBits) / part;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim * part];
	double* pvData = new double[factorDim * part];

	double* twData = new double[factorDim * part];
	double* tvData = new double[factorDim * part];

	double* cwData = new double[factorDim];
	double* pwData_ensemble = new double[factorDim]();
	double* twData_ensemble = new double[factorDim]();

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);


	double enccor, encauc, truecor, trueauc, plainauc, plaincor;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0, averplaincor = 0, averplainauc = 0;

	
	for (long fnum = 0; fnum < fold; ++fnum) {
		cout << " !!! START Ensemble NLGD for " << fnum + 1 << " FOLD !!! " << endl;

		for (long i = 0; i < sampleDimTest; ++i) {
			zDataTest[i] = zData[fnum * sampleDimTest + i];
		}
		for (long j = 0; j < fnum; ++j) {
			for (long i = 0; i < sampleDimTest; ++i) {
				zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
			}
		}
		for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
			zDataTrain[i - sampleDimTest] = zData[i];
		}


		timeutils.start("Encrypting zData...");
		cipherGD.encZData_ensemble(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
		timeutils.stop("zData encryption");


		pwData = EvaluatorUtils::randomRealArray_withsign(factorDim * part, bound);
		//	pwData = EvaluatorUtils::randomRealArray(factorDim * part, bound);

		for (long i = 0; i < factorDim * part; i++) {
			pvData[i] = pwData[i];
			tvData[i] = pwData[i];
			twData[i] = pwData[i];
		}
		cout << "initial vector = ";
		for(int i = 0; i < part; i++){
			for(int j = 0; j < factorDim; j++){
				cout << pwData[i * factorDim + j] << ", ";
			}
		}
		cout << endl;
		//encrypt wdata and vdata
		for(long i = 0; i < cnum; i++){
			complex<double>* pvec = new complex<double>[batch * part]();
			for(long j = 0; j < part; j++){
				for(long k = 0; k < batch; k++){
					if(i * batch + k < factorDim) pvec[j * batch + k].real(pwData[j * factorDim + i * batch + k]);
				}
			}
			timeutils.start("Encrypting wData and vData...");
			encWData[i] = scheme.encrypt(pvec, batch * part, wBits, logQ);
			encVData[i] = encWData[i];
			timeutils.stop("wData and vData encryption");
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		double** zData_tmp = new double*[(1 << sdimBits)];
		double** zData_ensemble = new double*[(1 << sdimBits)];
		for(int i = 0; i < (1 << sdimBits); i++){
			zData_tmp[i] = new double[factorDim];
			for(int j = 0; j < factorDim; j++){
				if(i < sampleDimTrain)
					zData_tmp[i][j] = zDataTrain[i][j];
				else
					zData_tmp[i][j] = zDataTrain[i - sampleDimTrain][j];
			}
		}
		for(int i = 0; i < part; i++){
			for(int k = 0; k < sampleDimEnsemble; k++){
				zData_ensemble[i * sampleDimEnsemble + k] = new double[factorDim];
				for(int j = 0; j < factorDim; j++){
					zData_ensemble[i * sampleDimEnsemble + k][j] = zData_tmp[part * k + i][j];  
				}
			}
		}

		for (long iter = 0; iter < numIter; ++iter) {
			cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimEnsemble : gammaUp / (iter - gammaDown) / sampleDimEnsemble;

			cout << "encWData.logq before: " << encWData[0].logq << endl;
			timeutils.start("Enc NLGD");
			cipherGD.encNLGDiteration_ensemble(part, kdeg, encZData, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
			timeutils.stop("Enc NLGD");
			cout << "encWData.logq after: " << encWData[0].logq << endl;

			for(int i = 0; i < part; i++){
				GD::plainNLGDiteration(kdeg, zData_ensemble + i * sampleDimEnsemble, pwData + i * factorDim, pvData + i * factorDim, factorDim, sampleDimEnsemble, gamma, eta);
				GD::trueNLGDiteration(zData_ensemble + i * sampleDimEnsemble , twData + i * factorDim, tvData + i * factorDim, factorDim, sampleDimEnsemble, gamma, eta);
			}

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}
		cout << "----ENCRYPTED-----" << endl;

		long part_bits = (long)ceil(log2(part));
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			for (long l = bBits; l < bBits + part_bits; l++) {
				Ciphertext rot = scheme.leftRotateByPo2(encWData[i], l);
				scheme.addAndEqual(encWData[i], rot);
			}
		}
		NTL_EXEC_RANGE_END;

		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);

		for(int i = 0; i < factorDim; i++){
			for(int j = 0; j < part; j++){
				pwData_ensemble[i] += pwData[factorDim * j + i];
				twData_ensemble[i] += twData[factorDim * j + i];
			}
		}

		cout << " !!! STOP  Ensemble NLGD !!! " << endl;
		cout << "------------------" << endl;


		for(long i = 0; i < factorDim; ++i) {
			cwData[i] /= part;
			pwData_ensemble[i] /= part;
			twData_ensemble[i] /= part;
		}
		GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
		cout << "------------------" << endl;

		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData_ensemble, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData_ensemble, factorDim, sampleDimTest, avertruecor, trueauc);
		cout << "------------------" << endl;
		averenccor += enccor;
		averencauc += encauc;
		averplaincor += plaincor;
		averplainauc += plainauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Encrypted correctness: " << averenccor / fold << "%" << endl;
	cout << "Encrypted AUC: " << averencauc / fold << endl;
	cout << "Plain correctness: " << averplaincor / fold << "%" << endl;
	cout << "Plain AUC: " << averplainauc / fold << endl;
	cout << "True correctness: " << avertruecor / fold << "%" << endl;
	cout << "True AUC: " << avertrueauc / fold << endl;
}

void TestGD::testPlainEnsembleNLGD(long fold, double** zData, long factorDim, long sampleDim,
	bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	//SHOULD USE POWER-OF-2 part!!!!!!
	double bound = 2;
	long part = 4;
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim <<",fold = " << fold <<", part = " << part << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", initial bound = " << bound << endl;


	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));



	// long wBits = 30;
	// long pBits = 20;
	// long lBits = 5;
	// long aBits = 3;
	// long kBits = (long)ceil(log2(kdeg));

	// long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
	// 		(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	// long logN = TestGD::suggestLogN(80, logQ);
	// long bBits = min(logN - 1 - sdimBits, fdimBits);
	// long batch =  1 << bBits;
	// long sBits = sdimBits + bBits;
	// long slots =  1 << sBits;
	long sampleDimEnsemble = (1 << sdimBits) / part;

	// cout << "batch = " << batch << ", slots = " << slots << endl;

	// cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	// cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;

	double* pwData = new double[factorDim * part];
	double* pvData = new double[factorDim * part];

	double* twData = new double[factorDim * part];
	double* tvData = new double[factorDim * part];

	double* pwData_ensemble = new double[factorDim]();
	double* twData_ensemble = new double[factorDim]();

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);


	double truecor, trueauc, plainauc, plaincor;
	double avertruecor = 0, avertrueauc = 0, averplaincor = 0, averplainauc = 0;


	for (long fnum = 0; fnum < fold; ++fnum) {
		cout << " !!! START Ensemble NLGD for " << fnum + 1 << " FOLD !!! " << endl;

		for (long i = 0; i < sampleDimTest; ++i) {
			zDataTest[i] = zData[fnum * sampleDimTest + i];
		}
		for (long j = 0; j < fnum; ++j) {
			for (long i = 0; i < sampleDimTest; ++i) {
				zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
			}
		}
		for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
			zDataTrain[i - sampleDimTest] = zData[i];
		}


		pwData = EvaluatorUtils::randomRealArray_withsign(factorDim * part, bound);
	//	pwData = EvaluatorUtils::randomRealArray(factorDim * part, bound);

		for (long i = 0; i < factorDim * part; i++) {
			pvData[i] = pwData[i];
			tvData[i] = pwData[i];
			twData[i] = pwData[i];
		}
		// cout << "initial vector = ";
		// for(int i = 0; i < part; i++){
		// 	for(int j = 0; j < factorDim; j++){
		// 		cout << pwData[i * factorDim + j] << ", ";
		// 	}
		// }
		// cout << endl;


		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		double** zData_tmp = new double*[(1 << sdimBits)];
		double** zData_ensemble = new double*[(1 << sdimBits)];
		for(int i = 0; i < (1 << sdimBits); i++){
			zData_tmp[i] = new double[factorDim];
			for(int j = 0; j < factorDim; j++){
				if(i < sampleDimTrain)
					zData_tmp[i][j] = zDataTrain[i][j];
				else
					zData_tmp[i][j] = zDataTrain[i - sampleDimTrain][j];
			}
		}
		for(int i = 0; i < part; i++){
			for(int k = 0; k < sampleDimEnsemble; k++){
				zData_ensemble[i * sampleDimEnsemble + k] = new double[factorDim];
				for(int j = 0; j < factorDim; j++){
					zData_ensemble[i * sampleDimEnsemble + k][j] = zData_tmp[part * k + i][j];  
				}
			}
		}

		for (long iter = 0; iter < numIter; ++iter) {
	//		cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimEnsemble : gammaUp / (iter - gammaDown) / sampleDimEnsemble;

			for(int i = 0; i < part; i++){
				GD::plainNLGDiteration(kdeg, zData_ensemble + i * sampleDimEnsemble, pwData + i * factorDim, pvData + i * factorDim, factorDim, sampleDimEnsemble, gamma, eta);
				GD::trueNLGDiteration(zData_ensemble + i * sampleDimEnsemble , twData + i * factorDim, tvData + i * factorDim, factorDim, sampleDimEnsemble, gamma, eta);
			}

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
	//		cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}


		for(int i = 0; i < factorDim; i++){
			for(int j = 0; j < part; j++){
				pwData_ensemble[i] += pwData[factorDim * j + i];
				twData_ensemble[i] += twData[factorDim * j + i];
			}
		}

		cout << " !!! STOP  Ensemble NLGD !!! " << endl;
		cout << "------------------" << endl;


		for(long i = 0; i < factorDim; ++i) {
			pwData_ensemble[i] /= part;
			twData_ensemble[i] /= part;
		}

		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData_ensemble, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData_ensemble, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

	
		averplaincor += plaincor;
		averplainauc += plainauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Plain correctness: " << averplaincor / fold << "%" << endl;
	cout << "Plain AUC: " << averplainauc / fold << endl;
	cout << "True correctness: " << avertruecor / fold << "%" << endl;
	cout << "True AUC: " << avertrueauc / fold << endl;

}








void TestGD::testEncEnsembleNLGD2(long part, double** zData, double** zData_test, long factorDim, long sampleDim, long sampleDimTest,
	bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	//SHOULD USE POWER-OF-2 part!!!!!!
	double bound = 2;
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", part = " << part << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", initial bound = " << bound << endl;


	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch =  1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long sampleDimTrain = (1 << sdimBits) / part;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim * part];
	double* pvData = new double[factorDim * part];

	double* twData = new double[factorDim * part];
	double* tvData = new double[factorDim * part];

	double* cwData = new double[factorDim];
	double* pwData_ensemble = new double[factorDim]();
	double* twData_ensemble = new double[factorDim]();

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);


	double enccor, encauc, truecor, trueauc, plainauc, plaincor;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0, averplaincor = 0, averplainauc = 0;

	
	cout << " !!! START Ensemble NLGD !!! " << endl;


	timeutils.start("Encrypting zData...");
	cipherGD.encZData_ensemble(encZData, zData, slots, factorDim, sampleDim, batch, cnum, wBits, logQ);
	timeutils.stop("zData encryption");


	pwData = EvaluatorUtils::randomRealArray_withsign(factorDim * part, bound);
//	pwData = EvaluatorUtils::randomRealArray(factorDim * part, bound);

	for (long i = 0; i < factorDim * part; i++) {
		pvData[i] = pwData[i];
		tvData[i] = pwData[i];
		twData[i] = pwData[i];
	}
	// cout << "initial vector = ";
	// for(int i = 0; i < part; i++){
	// 	for(int j = 0; j < factorDim; j++){
	// 		cout << pwData[i * factorDim + j] << ", ";
	// 	}
	// }
	// cout << endl;
	//encrypt wdata and vdata
	for(long i = 0; i < cnum; i++){
		complex<double>* pvec = new complex<double>[batch * part]();
		for(long j = 0; j < part; j++){
			for(long k = 0; k < batch; k++){
				if(i * batch + k < factorDim) pvec[j * batch + k].real(pwData[j * factorDim + i * batch + k]);
			}
		}
		timeutils.start("Encrypting wData and vData...");
		encWData[i] = scheme.encrypt(pvec, batch * part, wBits, logQ);
		encVData[i] = encWData[i];
		timeutils.stop("wData and vData encryption");
	}

	//-----------------------------------------

	double alpha0, alpha1, eta, gamma;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	double** zData_tmp = new double*[(1 << sdimBits)];
	double** zData_ensemble = new double*[(1 << sdimBits)];
	for(int i = 0; i < (1 << sdimBits); i++){
		zData_tmp[i] = new double[factorDim];
		for(int j = 0; j < factorDim; j++){
			if(i < sampleDim)
				zData_tmp[i][j] = zData[i][j];
			else
				zData_tmp[i][j] = zData[i - sampleDim][j];
		}
	}
	for(int i = 0; i < part; i++){
		for(int k = 0; k < sampleDimTrain; k++){
			zData_ensemble[i * sampleDimTrain + k] = new double[factorDim];
			for(int j = 0; j < factorDim; j++){
				zData_ensemble[i * sampleDimTrain + k][j] = zData_tmp[part * k + i][j];  
			}
		}
	}

	for (long iter = 0; iter < numIter; ++iter) {
		cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		cout << "encWData.logq before: " << encWData[0].logq << endl;
		timeutils.start("Enc NLGD");
		cipherGD.encNLGDiteration_ensemble(part, kdeg, encZData, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
		timeutils.stop("Enc NLGD");
		cout << "encWData.logq after: " << encWData[0].logq << endl;

		for(int i = 0; i < part; i++){
			GD::plainNLGDiteration(kdeg, zData_ensemble + i * sampleDimTrain, pwData + i * factorDim, pvData + i * factorDim, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zData_ensemble + i * sampleDimTrain , twData + i * factorDim, tvData + i * factorDim, factorDim, sampleDimTrain, gamma, eta);
		}

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
	}
	cout << "----ENCRYPTED-----" << endl;

	long part_bits = (long)ceil(log2(part));
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < bBits + part_bits; l++) {
			Ciphertext rot = scheme.leftRotateByPo2(encWData[i], l);
			scheme.addAndEqual(encWData[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);

	for(int i = 0; i < factorDim; i++){
		for(int j = 0; j < part; j++){
			pwData_ensemble[i] += pwData[factorDim * j + i];
			twData_ensemble[i] += twData[factorDim * j + i];
		}
	}

	cout << " !!! STOP  Ensemble NLGD !!! " << endl;
	cout << "------------------" << endl;


	for(long i = 0; i < factorDim; ++i) {
		cwData[i] /= part;
		pwData_ensemble[i] /= part;
		twData_ensemble[i] /= part;
	}
	GD::calculateAUC(zData_test, cwData, factorDim, sampleDimTest, averenccor, averencauc);
	cout << "------------------" << endl;

	cout << "------PLAIN-------" << endl;
	GD::calculateAUC(zData_test, pwData_ensemble, factorDim, sampleDimTest, averplaincor, averplainauc);
	cout << "------------------" << endl;

	cout << "-------TRUE-------" << endl;
	GD::calculateAUC(zData_test, twData_ensemble, factorDim, sampleDimTest, avertruecor, avertrueauc);
	cout << "------------------" << endl;

	cout << "Encrypted correctness: " << averenccor  << "%" << endl;
	cout << "Encrypted AUC: " << averencauc << endl;
	cout << "True correctness: " << avertruecor << "%" << endl;
	cout << "True AUC: " << avertrueauc << endl;
}

void TestGD::testPlainEnsembleNLGD2(long part, double** zData, double** zData_test,long factorDim, long sampleDim, long sampleDimTest,
	bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	//SHOULD USE POWER-OF-2 part!!!!!!
	double bound = 1;
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", part = " << part << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", initial bound = " << bound << endl;


	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));

	// long wBits = 30;
	// long pBits = 20;
	// long lBits = 5;
	// long aBits = 3;
	// long kBits = (long)ceil(log2(kdeg));

	// long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
	// 		(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	// long logN = TestGD::suggestLogN(80, logQ);
	// long bBits = min(logN - 1 - sdimBits, fdimBits);
	// long batch =  1 << bBits;
	// long sBits = sdimBits + bBits;
	// long slots =  1 << sBits;
	long sampleDimTrain = (1 << sdimBits) / part;

	// cout << "batch = " << batch << ", slots = " << slots << endl;

	// cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	// cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;

	double* pwData = new double[factorDim * part];
	double* pvData = new double[factorDim * part];

	double* twData = new double[factorDim * part];
	double* tvData = new double[factorDim * part];

	double* pwData_ensemble = new double[factorDim]();
	double* twData_ensemble = new double[factorDim]();

	// GD::normalizeZData(zData, factorDim, sampleDim);
	// GD::shuffleZData(zData, factorDim, sampleDim);


	double enccor, encauc, truecor, trueauc, plainauc, plaincor;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0, averplaincor = 0, averplainauc = 0;

	
	cout << " !!! START Ensemble NLGD !!! " << endl;


	pwData = EvaluatorUtils::randomRealArray_withsign(factorDim * part, bound);
//	pwData = EvaluatorUtils::randomRealArray(factorDim * part, bound);

	for (long i = 0; i < factorDim * part; i++) {
		pvData[i] = pwData[i];
		tvData[i] = pwData[i];
		twData[i] = pwData[i];
	}
	// cout << "initial vector = ";
	// for(int i = 0; i < part; i++){
	// 	for(int j = 0; j < factorDim; j++){
	// 		cout << pwData[i * factorDim + j] << ", ";
	// 	}
	// }
	// cout << endl;


	//-----------------------------------------

	double alpha0, alpha1, eta, gamma;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	double** zData_tmp = new double*[(1 << sdimBits)];
	double** zData_ensemble = new double*[(1 << sdimBits)];
	for(int i = 0; i < (1 << sdimBits); i++){
		zData_tmp[i] = new double[factorDim];
		for(int j = 0; j < factorDim; j++){
			if(i < sampleDim)
				zData_tmp[i][j] = zData[i][j];
			else
				zData_tmp[i][j] = zData[i - sampleDim][j];
		}
	}
	for(int i = 0; i < part; i++){
		for(int k = 0; k < sampleDimTrain; k++){
			zData_ensemble[i * sampleDimTrain + k] = new double[factorDim];
			for(int j = 0; j < factorDim; j++){
				zData_ensemble[i * sampleDimTrain + k][j] = zData_tmp[part * k + i][j];  
			}
		}
	}

	for (long iter = 0; iter < numIter; ++iter) {
//		cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		for(int i = 0; i < part; i++){
			GD::plainNLGDiteration(kdeg, zData_ensemble + i * sampleDimTrain, pwData + i * factorDim, pvData + i * factorDim, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zData_ensemble + i * sampleDimTrain , twData + i * factorDim, tvData + i * factorDim, factorDim, sampleDimTrain, gamma, eta);
		}

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
//		cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;

		for(int i = 0; i < factorDim; i++){
			for(int j = 0; j < part; j++){
				pwData_ensemble[i] += pwData[factorDim * j + i];
				twData_ensemble[i] += twData[factorDim * j + i];
			}
		}
		for(long i = 0; i < factorDim; ++i) {
			pwData_ensemble[i] /= part;
			twData_ensemble[i] /= part;
		}
		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zData_test, pwData_ensemble, factorDim, sampleDimTest, averplaincor, averplainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zData_test, twData_ensemble, factorDim, sampleDimTest, avertruecor, avertrueauc);
		cout << "------------------" << endl;
		for(long i = 0; i < factorDim; ++i) {
			pwData_ensemble[i] = 0;
			twData_ensemble[i] = 0;
		}
	}



	cout << " !!! STOP  Ensemble NLGD !!! " << endl;
	cout << "------------------" << endl;


}

void TestGD::testNewNLGD(long part, double** zData, double** zData_test, long factorDim, long sampleDim, long sampleDimTest,
			bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", part = " << part << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long sampleDimTrain = sampleDim/part;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQ);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];
	double* temp = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	double* averwData = new double[factorDim]();
	double* averpwData = new double[factorDim]();
	double* avertwData = new double[factorDim]();

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	double **zDataTrain;

	zDataTrain = new double*[sampleDimTrain];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);

	double bound = 0.01;
	double enccor, encauc, truecor, trueauc, plainauc, plaincor;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0, averplaincor = 0, averplainauc = 0;

//	NTL_EXEC_RANGE(part, first, last);
	for (long fnum = 0; fnum < part; ++fnum) {
		cout << " !!! START " << fnum + 1 << " part !!! " << endl;


		for (long i = 0; i < sampleDimTrain; ++i) {
			zDataTrain[i] = zData[i + fnum * sampleDimTrain];
		}

		timeutils.start("Encrypting zData...");
		cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
		timeutils.stop("zData encryption");


		pwData = EvaluatorUtils::randomRealArray_withsign(factorDim, bound);
		for (long i = 0; i < factorDim; i++) {
			pvData[i] = pwData[i];
			tvData[i] = pwData[i];
			twData[i] = pwData[i];
			temp[i] = pwData[i];
		}
		//encrypt wdata and vdata
		for(long i = 0; i < cnum; i++){
			complex<double>* pvec = new complex<double>[batch]();
			for(long j = 0; j < batch; j++) {
				if(i * batch + j < factorDim) pvec[j].real(pwData[i * batch + j]);
			}
			encWData[i] = scheme.encrypt(pvec, batch, wBits, logQ);
		//	encWData[i].slots = slots;
			encVData[i] = encWData[i];
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

			cout << "encWData.logq before: " << encWData[0].logq << endl;
			timeutils.start("Enc NLGD");
			cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
			timeutils.stop("Enc NLGD");
			cout << "encWData.logq after: " << encWData[0].logq << endl;

			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}
		cout << "----ENCRYPTED-----" << endl;
		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
		for(long i = 0; i < factorDim; ++i) {
			averwData[i] += cwData[i];
			averpwData[i] += pwData[i];
			avertwData[i] += twData[i];
		}

		cout << " !!! STOP " << fnum + 1 << " PART !!! " << endl;
		cout << "------------------" << endl;
	}
//	NTL_EXEC_RANGE_END;

	for(long i = 0; i < factorDim; ++i) {
		averwData[i] /= part;
		averpwData[i] /= part;
		avertwData[i] /= part;
	}
	GD::calculateAUC(zData_test, averwData, factorDim, sampleDimTest, averenccor, averencauc);
	cout << "------------------" << endl;

	cout << "------PLAIN-------" << endl;
	GD::calculateAUC(zData_test, averpwData, factorDim, sampleDimTest, averplaincor, averplainauc);
	cout << "------------------" << endl;

	cout << "-------TRUE-------" << endl;
	GD::calculateAUC(zData_test, avertwData, factorDim, sampleDimTest, avertruecor, avertrueauc);
	cout << "------------------" << endl;

	cout << "Average Encrypted correctness: " << averenccor  << "%" << endl;
	cout << "Average Encrypted AUC: " << averencauc << endl;
	cout << "Average True correctness: " << avertruecor << "%" << endl;
	cout << "Average True AUC: " << avertrueauc << endl;
}



