package myHMMV4;

import java.util.*;
import java.lang.Math;

public class DTW {

	/** number of states */
	public int numStates;

	/** size of output vocabulary */
	public int sigmaSize;

	/** size of output ordinary vocabulary */
	public int ordinaryVocabSize;

	/** weight of non ordinary (words in glossary) vocabulary */
	public double keywordWeight;
	
	/** extract feature from first N sentences; 0 means using all */
	public int segWordsMax;

	/** emission probabilities */
	public double b[][];

	/** language model smoothing */
	public double smoothing;

	/** model parameters : numStates * sigmaSize */
	public double distribution[][];

	/** align to last states (slides, bullet, etc.) */
	public boolean toTail;
	
	/** compute language model */
	List<HashMap<Integer, Double>> distrMaps;
	double sumLogP[];
	double sumCount[];

	public static double LZERO = -1.0E10; /* ~log(0) */
	public static double LSMALL = -0.99E10;
	public static double MINEARG = -708.3; /* lowest exp() arg = log(MINLARG) */
	public static double MINLARG = 2.45E-308; /* lowest log() arg = exp(MINEARG) */
	public static double ZERO = 1e-20;
	
	public DTW(int numStates, int sigmaSize, int ordinaryVocabSize,
			double keywordWeight, int segWordsMax, double smoothing,
			boolean... params) {
		assert params.length <= 1;
		this.numStates = numStates;
		this.sigmaSize = sigmaSize;
		this.ordinaryVocabSize = ordinaryVocabSize;
		this.keywordWeight = keywordWeight;
		this.segWordsMax = segWordsMax;
		this.smoothing = smoothing;
		this.toTail = params.length > 0 ? params[0] : true;

		setDistribution(new double[numStates][sigmaSize]);
	}
	
	public void setDistribution(double d[][]) {
		this.distribution = d;
		computelogLM();
	}
	
	public void setSegWordsMax(int segWordsMax) {
		this.segWordsMax = segWordsMax;
	}

	public void setDistributionSparse(int[][] oId, double[][] oCount) {
		double[][] gamma = new double[numStates][numStates];
		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++)
				gamma[i][j] = 0;
		}
		for (int i = 0; i < numStates; i++)
			gamma[i][i] = 1;
		this.distribution = sparseOToDistribution(gamma, oId, oCount);
		computelogLM();
	}
	
	public void setDistributionSparseWindowing(int[][] oId, double[][] oCount, double expDecay) {
		double[][] gamma = new double[numStates][numStates];		
		for (int i = 0; i < numStates; i++) {
			double normalize = 0.0;
			for (int j = 0; j < numStates; j++) {
				gamma[i][j] = Math.exp(-expDecay * Math.abs(i - j));
				normalize += gamma[i][j];
			}
			for (int j = 0; j < numStates; j++)
				gamma[i][j] /= normalize;
		}		
		
		this.distribution = sparseOToDistribution(gamma, oId, oCount);
		computelogLM();
	}
	
	/** Log version addition */
	public double LAdd(double x, double y) {
		double diff;
		if (x < y) {
			double temp = x;
			x = y;
			y = temp;
		} // make sure x > y
		diff = y - x; // diff < 0
		if (diff < MINEARG)
			return (x <= LSMALL) ? LZERO : x;
		else
			return x + Math.log(1.0 + Math.exp(diff));
	}

	double LOG(double a) {
		if (a < MINLARG)
			return LZERO;
		else
			return Math.log(a);
	}

	/**
	 * implementation of HMM training, viterbi version. Adaptation, log version,
	 * hard decision alignment.
	 * 
	 * @param o
	 *            the training set, featureSize * n_sentence
	 * @param steps
	 *            the number of steps
	 * @param adaptation
	 *            ratio between origin and adapted corpus, only emission
	 * @param adaptStructure
	 *            whether to adapt transition (a) and initial (pi) probability
	 * @param w
	 *            window size
	 */
	public void trainAdaptationLogHard(int[][] oId, double[][] oCount,
			int steps, double adaptation, int w) {
		int T = oId.length;
		int state_seq[] = new int[T];
		double[][] gamma = new double[numStates][T];
		double[][] distributionAdapt;

		double[][] distributionOrigin = new double[distribution.length][];
		for (int i = 0; i < distribution.length; i++)
			distributionOrigin[i] = distribution[i].clone();

		for (int s = 0; s < steps; s++) {
			b = calBjOtLog(oId, oCount);
			
			double[][] DTW = new double[numStates + 1][T + 1];
			int[][] path = new int[numStates + 1][T + 1];
			for (int n = 0; n <= numStates; n++) {
				for (int m = 0; m <= T; m++)
					DTW[n][m] = MINEARG * 20;				
			}
			DTW[0][0] = 0;
			
			for (int n = 1; n <= numStates; n++) {
				double center = (T * (2*n - 1)) / (double) (numStates * 2); 
				for (int m = Math.max(1, ((int) center) - w); m <= Math.min(T,
						((int) center) + w); m++) {
					/*Math.min(DTW[n-1][m],    // skip slides: NOT OK - 0
                            DTW[n][m-1],    // multiple sentences to one slides page: OK - 1
                            DTW[n-1][m-1]);    // next pages slides: OK - 2*/
					if (DTW[n-1][m-1] < DTW[n][m-1]) {
						DTW[n][m] = b[n-1][m-1] + DTW[n][m-1];
						path[n][m] = n;
					}
					else {
						DTW[n][m] = b[n-1][m-1] + DTW[n-1][m-1];
						path[n][m] = n-1;
					}
				}
			}			
			
			/* backtrack */
			state_seq[T - 1] = numStates - 1;
			for (int t = T - 2; t >= 0; t--) {
				state_seq[t] = path[state_seq[t + 1] + 1][t + 2] - 1;
			}

			for (int t = 0; t < T; t++) {
				for (int i = 0; i < numStates; i++)
					gamma[i][t] = 0;
				gamma[state_seq[t]][t] = 1;
			}
			distributionAdapt = sparseOToDistribution(gamma, oId, oCount);

			for (int i = 0; i < numStates; i++) {
				for (int d = 0; d < sigmaSize; d++) {
					distribution[i][d] = adaptation * distributionAdapt[i][d]
							+ (1 - adaptation) * distributionOrigin[i][d];
				}
			}
			computelogLM();
		}
	}

	/** Viterbi log version */
	public Map<Double, int[]> viterbiLog(int[][] oId, double[][] oCount, int w) {
		int T = oId.length;
		int state_seq[] = new int[T];
		int[][] path = new int[numStates + 1][T + 1];

//		double[][] gamma = new double[numStates][T];
		//double[][] distributionAdapt;
/*
		double[][] distributionOrigin = new double[distribution.length][];
		for (int i = 0; i < distribution.length; i++)
			distributionOrigin[i] = distribution[i].clone();
*/

		b = calBjOtLog(oId, oCount); /* log version TBD */
		double[][] DTW = new double[numStates + 1][T + 1];
			
		for (int n = 0; n <= numStates; n++) {
			for (int m = 0; m <= T; m++)
				DTW[n][m] = MINEARG * 20;
		}
		DTW[0][0] = 0;
			
		for (int n = 1; n <= numStates; n++) {
			double center = (T * (2 * n - 1)) / (double) (numStates * 2);
			for (int m = Math.max(1, ((int) center) - w); m <= Math.min(T,
					((int) center) + w); m++) {
				/*
				 * Math.min(DTW[n-1][m], // skip slides: NOT OK - 0 DTW[n][m-1], //
				 * multiple sentences to one slides page: OK - 1 DTW[n-1][m-1]);
				 * // next pages slides: OK - 2
				 */
				if (DTW[n - 1][m - 1] < DTW[n][m - 1]) {
					DTW[n][m] = b[n - 1][m - 1] + DTW[n][m - 1];
					path[n][m] = n;
				} else {
					DTW[n][m] = b[n - 1][m - 1] + DTW[n - 1][m - 1];
					path[n][m] = n - 1;
				}
			}
		}
			
		/* backtrack */
		state_seq[T - 1] = numStates - 1;
		for (int t = T - 2; t >= 0; t--) {
			state_seq[t] = path[state_seq[t + 1] + 1][t + 2] - 1;
		}

		Map<Double, int[]> resultMap = new HashMap<Double, int[]>();
		resultMap.put(DTW[numStates][T], state_seq);
		return resultMap;
	}

	public double[][] sparseOToDistribution(double[][] gamma, int[][] oId,
			double[][] oCount) {
		double distribution1[][] = new double[numStates][sigmaSize];
		int T = oId.length;
		
		for (int i = 0; i < numStates; i++) {
			for (int k = 0; k < sigmaSize; k++) {
				distribution1[i][k] = 0;
			}
		}
		
		for (int t = 0; t < T; t++) {
			int nWords = oId[t].length;
			if (segWordsMax != 0 && segWordsMax < oId[t].length) nWords = segWordsMax;			
			for (int q = 0; q < nWords; q++) {
				for (int i = 0; i < numStates; i++)
					distribution1[i][oId[t][q]] += oCount[t][q] * gamma[i][t];
			}
		}
		return distribution1;
	}

	/** calculate bjot, actually use KL instead of probability */
	public double[][] calBjOt(int[][] oId, double[][] oCount) {
		int T = oId.length;
		double bjot[][] = new double[numStates][T];
		double num;
		double zeroCntP;

		for (int t = 0; t < T; t++) {
			int nWords = oId[t].length;
			if (segWordsMax != 0 && segWordsMax < oId[t].length) nWords = segWordsMax;
			double qLength = 0;
			double logP;

			HashMap<Integer, Double> qmap = new HashMap<Integer, Double>();
			for (int i = 0; i < nWords; i++) {
				if (oId[t][i] < ordinaryVocabSize) {
					if (qmap.containsKey(oId[t][i])) {
						qmap.put(oId[t][i], qmap.get(oId[t][i]) + oCount[t][i]);
					}
					else { qmap.put(oId[t][i], oCount[t][i]); }
					qLength += oCount[t][i];
				}
				else {
					if (qmap.containsKey(oId[t][i])) {
						qmap.put(oId[t][i], qmap.get(oId[t][i]) + oCount[t][i] * keywordWeight);
					}
					else { qmap.put(oId[t][i], oCount[t][i] * keywordWeight); }
					qLength += oCount[t][i] * keywordWeight;
				}
			}

			for (int j = 0; j < numStates; j++) {
				num = 0;
				zeroCntP = sumLogP[j];
				HashMap<Integer, Double> map = distrMaps.get(j);

				for (Object key : qmap.keySet()) {
					if (map.containsKey(key)) {
						logP = map.get(key);
					} else {
						logP = Math.log((smoothing / sigmaSize)
								/ (sumCount[j] + smoothing));
					}
					num += ((qmap.get(key) + smoothing / sigmaSize) / (qLength + smoothing))
							* logP;
					zeroCntP -= logP;
				}

				num += ((smoothing / sigmaSize) / (qLength + smoothing))
						* zeroCntP;
				bjot[j][t] = Math.exp(num);
			}
		}

		return bjot;
	}

	/** calculate bjot, actually use KL instead of probability, log version */
	public double[][] calBjOtLog(int[][] oId, double[][] oCount) {
		int T = oId.length;
		double bjot[][] = new double[numStates][T];
		double num;
		double zeroCntP;

		for (int t = 0; t < T; t++) {
			int nWords = oId[t].length;
			if (segWordsMax != 0 && segWordsMax < oId[t].length) nWords = segWordsMax;
			double qLength = 0;
			double logP;

			HashMap<Integer, Double> qmap = new HashMap<Integer, Double>();
			for (int i = 0; i < nWords; i++) {
				if (oId[t][i] < ordinaryVocabSize) {
					if (qmap.containsKey(oId[t][i])) {
						qmap.put(oId[t][i], qmap.get(oId[t][i]) + oCount[t][i]);
					}
					else { qmap.put(oId[t][i], oCount[t][i]); }
					qLength += oCount[t][i];
				}
				else {
					if (qmap.containsKey(oId[t][i])) {
						qmap.put(oId[t][i], qmap.get(oId[t][i]) + oCount[t][i] * keywordWeight);
					}
					else { qmap.put(oId[t][i], oCount[t][i] * keywordWeight); }
					qLength += oCount[t][i] * keywordWeight;
				}
			}
						
			for (int j = 0; j < numStates; j++) {
				num = 0;
				zeroCntP = sumLogP[j];
				HashMap<Integer, Double> map = distrMaps.get(j);

				for (Object key : qmap.keySet()) {
					if (map.containsKey(key)) {
						logP = map.get(key);
					} else {
						logP = Math.log((smoothing / sigmaSize)
								/ (sumCount[j] + smoothing));
					}
					num += ((qmap.get(key) + smoothing / sigmaSize) / (qLength + smoothing))
							* logP;
					zeroCntP -= logP;
				}

				num += ((smoothing / sigmaSize) / (qLength + smoothing))
						* zeroCntP;
				if (Double.isNaN(num) || num < MINEARG)
					num = MINEARG;
				bjot[j][t] = num;
			}
		}

		return bjot;
	}
	
	public void computelogLM() {
		double num = 0;
		/* calculate document log LM */
		distrMaps = new ArrayList<HashMap<Integer, Double>>();
		sumLogP = new double[numStates];
		sumCount = new double[numStates];
	
		for (int i = 0; i < numStates; i++) {
			HashMap<Integer, Double> map = new HashMap<Integer, Double>();
			for (int d = 0; d < ordinaryVocabSize; d++) {
				sumCount[i] += distribution[i][d];
			}
			for (int d = ordinaryVocabSize; d < sigmaSize; d++) {
				sumCount[i] += keywordWeight * distribution[i][d];
			}
	
			for (int d = 0; d < ordinaryVocabSize; d++) {
				if (distribution[i][d] != 0) {
					num = Math.log((distribution[i][d] + smoothing / sigmaSize)
							/ (sumCount[i] + smoothing));
					map.put(d, num);
					sumLogP[i] += num;
				} else {
					sumLogP[i] += Math.log((smoothing / sigmaSize)
							/ (sumCount[i] + smoothing));
				}
			}
			for (int d = ordinaryVocabSize; d < sigmaSize; d++) {
				if (distribution[i][d] != 0) {
					num = Math
							.log((keywordWeight * distribution[i][d] + smoothing
									/ sigmaSize)
									/ (sumCount[i] + smoothing));
					map.put(d, num);
					sumLogP[i] += num;
				} else {
					sumLogP[i] += Math.log((smoothing / sigmaSize)
							/ (sumCount[i] + smoothing));
				}
			}
			distrMaps.add(map);
		}
	}
	
	public double[][] computeStatePosterior(int [] state_seq) {
		int T = state_seq.length;
		double[][] gamma = new double[numStates][T]; 
		for (int t = 0; t < T; t++) {
			for (int i = 0; i < numStates; i++)
				gamma[i][t] = 0;
			gamma[state_seq[t]][t] = 1;
		}
		return gamma;
	}

	/** divides two doubles. 0 / 0 = 0! */
	public double divide(double n, double d) {
		if (n == 0)
			return 0;
		else
			return n / d;
	}
}
