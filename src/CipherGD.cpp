/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "CipherGD.h"
#include "GD.h"

#include <Ciphertext.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

void CipherGD::encZData(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ) {
	complex<double>* pzData = new complex<double>[slots];
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < sampleDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(zData[j][batch * i + l]);
			}
		}
		encZData[i] = scheme.encrypt(pzData, slots, wBits, logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < sampleDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pzData[batch * j + l].real(zData[j][batch * (cnum - 1) + l]);
		}
		for (long l = rest; l < batch; ++l) {
			pzData[batch * j + l] = 0;
		}
	}
	encZData[cnum - 1] = scheme.encrypt(pzData, slots, wBits, logQ);

	delete[] pzData;
}
void CipherGD::encZData_ensemble(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ) {
	complex<double>* pzData = new complex<double>[slots]();
	long sdimBits = (long)ceil(log2(sampleDim));
	long sdim = 1 << sdimBits;
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < sdim; ++j) {
			for (long l = 0; l < batch; ++l) {
				if(j < sampleDim) 
					pzData[batch * j + l].real(zData[j][batch * i + l]);
				else
					pzData[batch * j + l].real(zData[j - sampleDim][batch * i + l]);
			}
		}
		encZData[i] = scheme.encrypt(pzData, slots, wBits, logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < sdim; ++j) {
		for (long l = 0; l < rest; ++l) {
			if(j < sampleDim)
				pzData[batch * j + l].real(zData[j][batch * (cnum - 1) + l]);
			else
				pzData[batch * j + l].real(zData[j - sampleDim][batch * (cnum - 1) + l]);
		}
		for (long l = rest; l < batch; ++l) {
			pzData[batch * j + l] = 0;
		}
	}
	encZData[cnum - 1] = scheme.encrypt(pzData, slots, wBits, logQ);

	delete[] pzData;
}

void CipherGD::encWDataAverage(Ciphertext* encWData, Ciphertext* encZData, long cnum, long sBits, long bBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = encZData[i];
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encWData[i], l);
			scheme.addAndEqual(encWData[i], rot);
		}
		scheme.divByPo2AndEqual(encWData[i], sBits - bBits);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWVDataAverage(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encZData, long cnum, long sBits, long bBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = encZData[i];
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encWData[i], l);
			scheme.addAndEqual(encWData[i], rot);
		}
		scheme.divByPo2AndEqual(encWData[i], sBits - bBits);
		encVData[i] = encWData[i];
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWDataZero(Ciphertext* encWData, long cnum, long slots, long wBits, long logQ) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = scheme.encryptZeros(slots, wBits, logQ);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWVDataZero(Ciphertext* encWData, Ciphertext* encVData, long cnum, long slots, long wBits, long logQ) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = scheme.encryptZeros(slots, wBits, logQ);
		encVData[i] = encWData[i];
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWVDataRandom(Ciphertext* encWData, Ciphertext* encVData, long cnum, long slots, long wBits, long logQ, double bound) {
NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		double* randw = EvaluatorUtils::randomRealArray(slots, bound);
		encWData[i] = scheme.encrypt(randw, slots, wBits, logQ);
		encVData[i] = encWData[i];
		delete[] randw;
	}
	NTL_EXEC_RANGE_END;
}

ZZX CipherGD::generateAuxPoly(long slots, long batch, long pBits) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}

Ciphertext CipherGD::encInnerProduct(Ciphertext* encZData, Ciphertext* encWData, ZZX& poly, long cnum, long bBits, long wBits, long pBits) {
	Ciphertext* encIPvec = new Ciphertext[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
	//	cout << "hi6" << endl;
		encIPvec[i] = scheme.modDownTo(encZData[i], encWData[i].logq);
	//	cout << "hi7" << endl;
		scheme.multAndEqual(encIPvec[i], encWData[i]); // xy * w
	//	cout << "hi8" << endl;
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(encIPvec[i], l);
			scheme.addAndEqual(encIPvec[i], rot);
		}
	}
	NTL_EXEC_RANGE_END
	Ciphertext encIP = encIPvec[0];
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(encIP, encIPvec[i]);
	}

	scheme.reScaleByAndEqual(encIP, wBits);
	scheme.multByPolyAndEqual(encIP, poly, pBits);
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(encIP, l);
		scheme.addAndEqual(encIP, tmp);
	}
	delete[] encIPvec;
	return encIP;
}

void CipherGD::encSigmoid(long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits) {
	Ciphertext encIP2 = scheme.square(encIP);
	scheme.reScaleByAndEqual(encIP2, wBits);

	if(kdeg == 3) {
		scheme.addConstAndEqual(encIP2, degree3[1] / degree3[2], wBits - 2 * aBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			encGrad[i] = scheme.multByConst(encZData[i], gamma  * degree3[2], wBits + 3 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP2);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree3[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else if (kdeg == 5) {
		Ciphertext encIP4 = scheme.square(encIP2);
		scheme.multByConstAndEqual(encIP2, degree5[2] / degree5[3], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(encIP2, wBits);
		scheme.reScaleByAndEqual(encIP4, wBits);
		scheme.addAndEqual(encIP4, encIP2);
		scheme.addConstAndEqual(encIP4, degree5[1] / degree5[3], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree5[3], wBits + 5 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP4.logq);
			scheme.multAndEqual(encGrad[i], encIP4);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree5[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else {
		Ciphertext encIP4 = scheme.square(encIP2);
		scheme.reScaleByAndEqual(encIP4, wBits);
		Ciphertext encIP2c = scheme.multByConst(encIP2, degree7[3] / degree7[4], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(encIP2c, wBits);
		scheme.addAndEqual(encIP4, encIP2c);
		scheme.addConstAndEqual(encIP4, degree7[2] / degree7[4], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree7[1], wBits + aBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encIP.logq);
			scheme.multAndEqual(tmp, encIP);
			scheme.reScaleByAndEqual(tmp, wBits);

			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree7[0], wBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], tmp.logq);
			scheme.addAndEqual(tmp, encGrad[i]);

			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree7[4], wBits + 7 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP2);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP4);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;

	}

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext tmp = scheme.leftRotateByPo2(encGrad[i], l);
			scheme.addAndEqual(encGrad[i], tmp);
		}
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::enclinHinge(long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits) {
	// Inefficient algorithm : g(x) - g(x-1) 
	Ciphertext encIP2 = scheme.square(encIP);
	scheme.reScaleByAndEqual(encIP2, wBits);
	Ciphertext encIP2_addconst = scheme.addConst(encIP2, degree4[2] / degree4[3], wBits - 2 * aBits);
	scheme.addConstAndEqual(encIP, degree4[0] / degree4[1], wBits - aBits);

	//preparation for x-1
	Ciphertext encIP_minus = scheme.addConst(encIP, -1.0, wBits - aBits);
	Ciphertext encIP_minus2 = scheme.square(encIP_minus);
	scheme.reScaleByAndEqual(encIP_minus2, wBits);
	Ciphertext encIP_minus2_addconst = scheme.addConst(encIP_minus2, degree4[2] / degree4[3], wBits - 2 * aBits);
	scheme.addConstAndEqual(encIP_minus, degree4[0] / degree4[1], wBits - aBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		//a_1 + a_2 x^2 + a_3 x^4
		encGrad[i] = scheme.multByConst(encZData[i], gamma  * degree4[3], wBits + 4 * aBits);
		scheme.reScaleByAndEqual(encGrad[i], wBits);
		scheme.modDownToAndEqual(encGrad[i], encIP2.logq);
		scheme.multAndEqual(encGrad[i], encIP2);
		scheme.reScaleByAndEqual(encGrad[i], wBits);
		scheme.multAndEqual(encGrad[i], encIP2_addconst);
		scheme.reScaleByAndEqual(encGrad[i], wBits);


		Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree4[1], wBits);
		scheme.reScaleByAndEqual(tmp, wBits);
		scheme.modDownToAndEqual(tmp, encGrad[i].logq);
		scheme.addAndEqual(encGrad[i], tmp);

		//a_2 (x-1)^2 + a_3 (x-1)^4
		tmp = scheme.multByConst(encZData[i], gamma  * degree4[3], wBits + 4 * aBits);
		scheme.reScaleByAndEqual(tmp, wBits);
		scheme.modDownToAndEqual(tmp, encIP_minus2.logq);
		scheme.multAndEqual(tmp, encIP_minus2);
		scheme.reScaleByAndEqual(tmp, wBits);
		scheme.multAndEqual(tmp, encIP_minus2_addconst);
		scheme.reScaleByAndEqual(tmp, wBits);

		//g(x) - g(x-1)
		scheme.subAndEqual(encGrad[i], tmp);
	}
	NTL_EXEC_RANGE_END;

}

void CipherGD::encSigmoid_ensemble(long part, long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits) {
	Ciphertext encIP2 = scheme.square(encIP);
	scheme.reScaleByAndEqual(encIP2, wBits);

	if(kdeg == 3) {
		scheme.addConstAndEqual(encIP2, degree3[1] / degree3[2], wBits - 2 * aBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			encGrad[i] = scheme.multByConst(encZData[i], gamma  * degree3[2], wBits + 3 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP2);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree3[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else if (kdeg == 5) {
		Ciphertext encIP4 = scheme.square(encIP2);
		scheme.multByConstAndEqual(encIP2, degree5[2] / degree5[3], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(encIP2, wBits);
		scheme.reScaleByAndEqual(encIP4, wBits);
		scheme.addAndEqual(encIP4, encIP2);
		scheme.addConstAndEqual(encIP4, degree5[1] / degree5[3], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree5[3], wBits + 5 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP4.logq);
			scheme.multAndEqual(encGrad[i], encIP4);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree5[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else {
		Ciphertext encIP4 = scheme.square(encIP2);
		scheme.reScaleByAndEqual(encIP4, wBits);
		Ciphertext encIP2c = scheme.multByConst(encIP2, degree7[3] / degree7[4], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(encIP2c, wBits);
		scheme.addAndEqual(encIP4, encIP2c);
		scheme.addConstAndEqual(encIP4, degree7[2] / degree7[4], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Ciphertext tmp = scheme.multByConst(encZData[i], gamma * degree7[1], wBits + aBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, encIP.logq);
			scheme.multAndEqual(tmp, encIP);
			scheme.reScaleByAndEqual(tmp, wBits);

			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree7[0], wBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], tmp.logq);
			scheme.addAndEqual(tmp, encGrad[i]);

			encGrad[i] = scheme.multByConst(encZData[i], gamma * degree7[4], wBits + 7 * aBits);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.modDownToAndEqual(encGrad[i], encIP.logq);
			scheme.multAndEqual(encGrad[i], encIP);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP2);
			scheme.reScaleByAndEqual(encGrad[i], wBits);
			scheme.multAndEqual(encGrad[i], encIP4);
			scheme.reScaleByAndEqual(encGrad[i], wBits);

			scheme.modDownToAndEqual(tmp, encGrad[i].logq);
			scheme.addAndEqual(encGrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;

	}

	long partBits = (long)ceil(log2(part));
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits + partBits; l < sBits; ++l) {
			Ciphertext tmp = scheme.leftRotateByPo2(encGrad[i], l);
			scheme.addAndEqual(encGrad[i], tmp);
		}
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encLGDstep(Ciphertext* encWData, Ciphertext* encGrad, long cnum) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(encWData[i], encGrad[i].logq);
		scheme.subAndEqual(encWData[i], encGrad[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encMLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long pBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multByConstAndEqual(encVData[i], eta, pBits);
		scheme.reScaleByAndEqual(encVData[i], pBits);
		scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);
		scheme.addAndEqual(encVData[i], encGrad[i]);
		scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
		scheme.subAndEqual(encWData[i], encVData[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encNLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long pBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);
		Ciphertext ctmpw = scheme.sub(encVData[i], encGrad[i]);
		encVData[i] = scheme.multByConst(ctmpw, 1. - eta, pBits);
		scheme.reScaleByAndEqual(encVData[i], pBits);
		scheme.multByConstAndEqual(encWData[i], eta, pBits);
		scheme.reScaleByAndEqual(encWData[i], pBits);
		scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
		scheme.addAndEqual(encVData[i], encWData[i]);
		encWData[i] = ctmpw;
	}
	NTL_EXEC_RANGE_END;
}


void CipherGD::encLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, ZZX& poly, long cnum, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* encGrad = new Ciphertext[cnum];
	Ciphertext encIP = encInnerProduct(encZData, encWData, poly, cnum, bBits, wBits, pBits);
	encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits);
	encLGDstep(encWData, encGrad, cnum);
	delete[] encGrad;
}

void CipherGD::encMLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* encGrad = new Ciphertext[cnum];
	Ciphertext encIP = encInnerProduct(encZData, encWData, poly, cnum, bBits, wBits, pBits);
	encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits);
	encMLGDstep(encWData, encVData, encGrad, eta, cnum, wBits);
	delete[] encGrad;
}

void CipherGD::encNLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* encGrad = new Ciphertext[cnum];
	Ciphertext encIP = encInnerProduct(encZData, encVData, poly, cnum, bBits, wBits, pBits);
	scheme.reScaleByAndEqual(encIP, pBits + aBits);
	encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits);
//	enclinHinge(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits);
	// complex<double>* dcw = scheme.decrypt(secretKey, encGrad[0]);
	// for(int i = 0; i < 100; i++){
	// 	cout << dcw[i].real() << ",";
	// }
	// cout << endl;
	encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits);
	delete[] encGrad;
}

void CipherGD::encNLGDiteration_ensemble(long part, long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* encGrad = new Ciphertext[cnum];
// 	cout << "hi1" << endl;
	Ciphertext encIP = encInnerProduct(encZData, encVData, poly, cnum, bBits, wBits, pBits);
//	cout << "hi2" << endl;
	scheme.reScaleByAndEqual(encIP, pBits + aBits);
//	cout << "hi3" << endl;
	encSigmoid_ensemble(part, kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits);
//	cout << "hi4" << endl;
	encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits);
//	cout << "hi5" << endl;
	delete[] encGrad;
}

void CipherGD::decWData(double* wData, Ciphertext* encWData, long factorDim, long batch, long cnum, long wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		complex<double>* dcw = scheme.decrypt(secretKey, encWData[i]);
		for (long j = 0; j < batch; ++j) {
			wData[batch * i + j] = dcw[j].real();
		}
		delete[] dcw;
	}
	complex<double>* dcw = scheme.decrypt(secretKey, encWData[cnum-1]);
	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < rest; ++j) {
		wData[batch * (cnum - 1) + j] = dcw[j].real();
	}
	delete[] dcw;
}

// void CipherGD::decWData_ensemble(long part, double* wData, Ciphertext* encWData, long factorDim, long batch, long cnum, long wBits) {
// 	wData = new double[factorDim]();
// 	for (long i = 0; i < (cnum - 1); ++i) {
// 		complex<double>* dcw = scheme.decrypt(secretKey, encWData[i]);
// 		for(long j = 0; j < part; j++){
// 			for (long k = 0; k < batch; ++k) {
// 				wData[batch * i + k] += dcw[j * batch + k].real();
// 			}
// 		}
// 		delete[] dcw;
// 	}
// 	complex<double>* dcw = scheme.decrypt(secretKey, encWData[cnum-1]);
// 	long rest = factorDim - batch * (cnum - 1);
// 	for(long j = 0; j < part; j++){
// 		for (long k = 0; k < rest; ++k) {
// 			wData[batch * (cnum - 1) + k] += dcw[j * batch + k].real();
// 		}
// 	}
// 	delete[] dcw;
// }
