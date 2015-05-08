package myHMMV4;

import java.text.*;
import java.util.*;
import java.lang.Math;

/**
 * This class implements a Hidden Markov Model, as well as the Baum-Welch
 * Algorithm for training HMMs.
 * 
 */

public class HMM {
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

	/** initial state probabilities */
	public double pi[];

	/** transition probabilities */
	public double a[][];

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

	/**
	 * initializes an HMM.
	 * 
	 * @param numStates
	 *            number of states
	 * @param sigmaSize
	 *            size of output vocabulary
	 */
	public HMM(int numStates, int sigmaSize, int ordinaryVocabSize,
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

		setPi(new double[numStates]);
		setA(new double[numStates][numStates]);
		setDistribution(new double[numStates][sigmaSize]);
	}

	public void setA(double a[][]) {
		this.a = a;
	}

	public void setPi(double pi[]) {
		this.pi = pi;
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
	 * implementation of the Baum-Welch Algorithm for HMMs.
	 * 
	 * @param o
	 *            the training set, featureSize * n_sentence
	 * @param steps
	 *            the number of steps
	 */
	public void train(int[][] oId, double[][] oCount, int steps) {
		int T = oId.length;
		double[][] fwd;
		double[][] bwd;
		double[][] gamma;
		double[][][] xi;

		double pi1[] = new double[numStates];
		double a1[][] = new double[numStates][numStates];

		for (int s = 0; s < steps; s++) {
			/*
			 * calculation of Forward-Backward Variables from the current model
			 */
			b = calBjOt(oId, oCount);
			fwd = forwardProc();
			bwd = backwardProc();
			gamma = calGamma(fwd, bwd);
			xi = calXi(fwd, bwd);

			/* re-estimation of initial state probabilities */
			for (int i = 0; i < numStates; i++)
				pi1[i] = gamma[i][0];

			/* re-estimation of transition probabilities */
			for (int i = 0; i < numStates; i++) {
				for (int j = 0; j < numStates; j++) {
					double num = 0;
					double denom = 0;
					for (int t = 0; t <= T - 1; t++) {
						num += xi[i][j][t];
						denom += gamma[i][t];
					}
					a1[i][j] = divide(num, denom);
				}
			}

			/* re-estimation of emission probabilities */
			distribution = sparseOToDistribution(gamma, oId, oCount);
			pi = pi1;
			a = a1;
			computelogLM();
		}
	}

	/**
	 * implementation of the Baum-Welch Algorithm for HMMs. Adaptation, log
	 * version.
	 * 
	 * @param o
	 *            the training set, featureSize * n_sentence
	 * @param steps
	 *            the number of steps
	 * @param adaptation
	 *            ratio between origin and adapted corpus, only emission
	 * @param adaptStructure
	 *            whether to adapt transition (a) and initial (pi) probability
	 */
	public void trainAdaptationLog(int[][] oId, double[][] oCount, int steps, double adaptation,
			boolean adaptStructure) {
		int T = oId.length;
		double[][] fwd;
		double[][] bwd;
		double[][] gamma;
		double[][][] xi;
		double[][] distributionAdapt;

		double pi1[] = new double[numStates];
		double a1[][] = new double[numStates][numStates];
		double[][] distributionOrigin = new double[distribution.length][];
		for (int i = 0; i < distribution.length; i++)
			distributionOrigin[i] = distribution[i].clone();

		for (int s = 0; s < steps; s++) {
			/*
			 * calculation of Forward-Backward Variables from the current model
			 */
			b = calBjOtLog(oId, oCount);
			fwd = forwardProcLog();
			bwd = backwardProcLog();
			gamma = calGammaLog(fwd, bwd);
			xi = calXiLog(fwd, bwd);

			if (adaptStructure) {
				/* re-estimation of initial state probabilities */
				for (int i = 0; i < numStates; i++)
					pi1[i] = Math.exp(gamma[i][0]);
								
				/* re-estimation of transition probabilities */
				for (int i = 0; i < numStates; i++) {
					for (int j = 0; j < numStates; j++) {
						double num = LZERO;
						double denom = LZERO;
						for (int t = 0; t <= T - 1; t++) {
							num = LAdd(xi[i][j][t], num);
							denom = LAdd(gamma[i][t], denom);
						}
						a1[i][j] = Math.exp(num - denom);
					}
				}				
			}
			for (int t = 0; t < T; t++) {
				for (int i = 0; i < numStates; i++) {
					gamma[i][t] = Math.exp(gamma[i][t]);
				}
			}

			/* re-estimation of emission probabilities */
			distributionAdapt = sparseOToDistribution(gamma, oId, oCount);
			for (int i = 0; i < numStates; i++) {
				for (int d = 0; d < sigmaSize; d++) {
					distribution[i][d] = adaptation * distributionAdapt[i][d]
							+ (1 - adaptation) * distributionOrigin[i][d];
				}
			}
			computelogLM();
			pi = pi1;
			a = a1;
		}
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
	 */
	public void trainAdaptationLogHard(int[][] oId, double[][] oCount, int steps, double adaptation) {
		int T = oId.length;
		int state_seq[] = new int[T];
		double[][] delta = new double[numStates][T];
		int[][] path = new int[numStates][T];
		double[][] gamma = new double[numStates][T];
		double bestScore;
		double[][] distributionAdapt;
		
		double[][] distributionOrigin = new double[distribution.length][];
		for (int i = 0; i < distribution.length; i++)
			distributionOrigin[i] = distribution[i].clone();

		for (int s = 0; s < steps; s++) {
			b = calBjOtLog(oId, oCount);
			/* initialization (time 0) */
			for (int i = 0; i < numStates; i++)
				delta[i][0] = Math.log(pi[i]) + b[i][0];

			/* induction */
			for (int t = 0; t <= T - 2; t++) {
				for (int j = 0; j < numStates; j++) {
					delta[j][t + 1] = LZERO;
					path[j][t + 1] = -1;
					for (int i = 0; i < numStates; i++) {
						bestScore = delta[i][t] + Math.log(a[i][j]);
						if (bestScore > delta[j][t + 1]) {
							delta[j][t + 1] = bestScore;
							path[j][t + 1] = i;
						}
					}
					delta[j][t + 1] += b[j][t + 1];
				}
			}

			/* backtrack */
			state_seq[T - 1] = numStates - 1;
			bestScore = delta[numStates - 1][T - 1];

			if (!toTail) {
				for (int i = numStates - 2; i >= 0; i--) {
					if (delta[i][T - 1] > bestScore) {
						bestScore = delta[i][T - 1];
						state_seq[T - 1] = i;
					}
				}
			}
			for (int t = T - 2; t >= 0; t--) {
				state_seq[t] = path[state_seq[t + 1]][t + 1];
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

	/** Viterbi */
	public Map<Double, int[]> viterbi(int[][] oId, double[][] oCount) {
		int T = oId.length;
		int state_seq[] = new int[T];
		double[][] delta = new double[numStates][T];
		int[][] path = new int[numStates][T];
		double bestScore;

		b = calBjOt(oId, oCount);

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			delta[i][0] = pi[i] * b[i][0];

		/* induction */
		for (int t = 0; t <= T - 2; t++) {
			for (int j = 0; j < numStates; j++) {
				delta[j][t + 1] = 0;
				path[j][t + 1] = -1;
				for (int i = 0; i < numStates; i++) {
					bestScore = delta[i][t] * a[i][j];
					if (bestScore > delta[j][t + 1]) {
						delta[j][t + 1] = bestScore;
						path[j][t + 1] = i;
					}

				}
				delta[j][t + 1] *= b[j][t + 1];
			}
		}
		/* backtrack */
		state_seq[T - 1] = numStates - 1;
		bestScore = delta[numStates - 1][T - 1];

		if (!toTail) {
			for (int i = numStates - 2; i >= 0; i--) {
				if (delta[i][T - 1] > bestScore) {
					bestScore = delta[i][T - 1];
					state_seq[T - 1] = i;
				}
			}
		}
		for (int t = T - 2; t >= 0; t--)
			state_seq[t] = path[state_seq[t + 1]][t + 1];

		Map<Double, int[]> resultMap = new HashMap<Double, int[]>();
		resultMap.put(bestScore, state_seq);
		return resultMap;
	}

	/** Viterbi log version */
	public Map<Double, int[]> viterbiLog(int[][] oId, double[][] oCount) {
		int T = oId.length;
		int state_seq[] = new int[T];
		double[][] delta = new double[numStates][T];
		int[][] path = new int[numStates][T];
		double bestScore;

		b = calBjOtLog(oId, oCount); /* log version TBD */

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			delta[i][0] = Math.log(pi[i]) + b[i][0];

		/* induction */
		for (int t = 0; t <= T - 2; t++) {
			for (int j = 0; j < numStates; j++) {
				delta[j][t + 1] = LZERO;
				path[j][t + 1] = -1;
				for (int i = 0; i < numStates; i++) {
					bestScore = delta[i][t] + Math.log(a[i][j]);
					if (bestScore > delta[j][t + 1]) {
						delta[j][t + 1] = bestScore;
						path[j][t + 1] = i;
					}

				}
				delta[j][t + 1] += b[j][t + 1];
			}
		}
		/* backtrack */
		state_seq[T - 1] = numStates - 1;
		bestScore = delta[numStates - 1][T - 1];
		if (!toTail) {
			for (int i = numStates - 2; i >= 0; i--) {
				if (delta[i][T - 1] > bestScore) {
					bestScore = delta[i][T - 1];
					state_seq[T - 1] = i;
				}
			}
		}
		for (int t = T - 2; t >= 0; t--)
			state_seq[t] = path[state_seq[t + 1]][t + 1];

		Map<Double, int[]> resultMap = new HashMap<Double, int[]>();
		resultMap.put(bestScore, state_seq);
		return resultMap;
	}

	/** Viterbi log version with constraint */
	public Map<Double, int[]> viterbiLogConstraint(int[][] oId, double[][] oCount,
			List<ArrayList<Integer>> validState) {
		int T = oId.length;
		assert T == validState.size();
		
		int state_seq[] = new int[T];
		double[][] delta = new double[numStates][T];
		int[][] path = new int[numStates][T];
		double bestScore;

		b = calBjOtLog(oId, oCount); /* log version TBD */
		
		for (int i = 0; i < numStates; i++) {
			for (int t = 0; t < T; t++) {
				delta[i][t] = LZERO;
				path[i][t] = -1;
			}
		}
		/* initialization (time 0) */
		assert validState.get(0).size() != 0;
		for (int i = 0; i < validState.get(0).size(); i++)
			delta[validState.get(0).get(i)][0] = Math.log(pi[validState.get(0)
					.get(i)]) + b[validState.get(0).get(i)][0];

		/* induction */
		for (int t = 0; t <= T - 2; t++) {
			assert validState.get(t+1).size() != 0;
			for (int j = 0; j < validState.get(t+1).size(); j++) {
				int toState = validState.get(t+1).get(j);
				delta[toState][t + 1] = LZERO;
				path[toState][t + 1] = -1;
				for (int i = 0; i < validState.get(t).size(); i++) {
					int fromState = validState.get(t).get(i);
					bestScore = delta[fromState][t] + Math.log(a[fromState][toState]);
					if (bestScore > delta[toState][t + 1]) {
						delta[toState][t + 1] = bestScore;
						path[toState][t + 1] = fromState;
					}
				}
				delta[toState][t + 1] += b[toState][t + 1];
				if (path[toState][t + 1] == -1) {
					System.out.println("Computation Error: cannot find valid state");
					path[toState][t + 1] = validState.get(t).get(0);
				}
			}
		}
		/* backtrack */
		state_seq[T - 1] = numStates - 1;
		bestScore = delta[numStates - 1][T - 1];
		if (!toTail) {
			for (int i = numStates - 2; i >= 0; i--) {
				if (delta[i][T - 1] > bestScore) {
					bestScore = delta[i][T - 1];
					state_seq[T - 1] = i;
				}
			}
		}
		for (int t = T - 2; t >= 0; t--)
			state_seq[t] = path[state_seq[t + 1]][t + 1];

		Map<Double, int[]> resultMap = new HashMap<Double, int[]>();
		resultMap.put(bestScore, state_seq);
		return resultMap;
	}
	
	public double[][] computeStatePosterior(int[][] oId, double[][] oCount) {
		double[][] fwd;
		double[][] bwd;

		b = calBjOtLog(oId, oCount);
		fwd = forwardProcLog();
		bwd = backwardProcLog();
		return calGammaLog(fwd, bwd);		
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

	/**
	 * calculation of Forward-Variables f(i,t) for state i at time t for output
	 * sequence O with the current HMM parameters
	 * 
	 * @return an array f(i,t) over states and times, containing the
	 *         Forward-variables.
	 */
	public double[][] forwardProc() {
		int T = b[0].length;
		double[][] fwd = new double[numStates][T];

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			fwd[i][0] = pi[i] * b[i][0];

		/* induction */
		for (int t = 0; t <= T - 2; t++) {
			for (int j = 0; j < numStates; j++) {
				fwd[j][t + 1] = 0;
				for (int i = 0; i < numStates; i++)
					fwd[j][t + 1] += (fwd[i][t] * a[i][j]);
				fwd[j][t + 1] *= b[j][t + 1];
			}
		}

		return fwd;
	}

	/**
	 * calculation of Forward-Variables f(i,t) for state i at time t for output
	 * sequence O with the current HMM parameters log version
	 * 
	 * @return an array f(i,t) over states and times, containing the
	 *         Forward-variables.
	 */
	public double[][] forwardProcLog() {
		int T = b[0].length;
		double[][] fwd = new double[numStates][T];

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			fwd[i][0] = Math.log(pi[i]) + b[i][0];

		/* induction */
		for (int t = 0; t <= T - 2; t++) {
			for (int j = 0; j < numStates; j++) {
				fwd[j][t + 1] = LZERO;
				for (int i = 0; i < numStates; i++)
					fwd[j][t + 1] = LAdd(fwd[i][t] + Math.log(a[i][j]),
							fwd[j][t + 1]);
				fwd[j][t + 1] += b[j][t + 1];
			}
		}

		return fwd;
	}

	/**
	 * calculation of Backward-Variables b(i,t) for state i at time t for output
	 * sequence O with the current HMM parameters
	 * 
	 * @return an array b(i,t) over states and times, containing the
	 *         Backward-Variables.
	 */
	public double[][] backwardProc() {
		int T = b[0].length;
		double[][] bwd = new double[numStates][T];

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			bwd[i][T - 1] = 1;

		/* induction */
		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < numStates; i++) {
				bwd[i][t] = 0;
				for (int j = 0; j < numStates; j++)
					bwd[i][t] += (bwd[j][t + 1] * a[i][j] * b[j][t + 1]);
			}
		}

		return bwd;
	}

	/**
	 * calculation of Backward-Variables b(i,t) for state i at time t for output
	 * sequence O with the current HMM parameters Log version
	 * 
	 * @return an array b(i,t) over states and times, containing the
	 *         Backward-Variables.
	 */
	public double[][] backwardProcLog() {
		int T = b[0].length;
		double[][] bwd = new double[numStates][T];

		/* initialization (time 0) */
		for (int i = 0; i < numStates; i++)
			bwd[i][T - 1] = 0;

		/* induction */
		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < numStates; i++) {
				bwd[i][t] = LZERO;
				for (int j = 0; j < numStates; j++)
					bwd[i][t] = LAdd(bwd[j][t + 1] + Math.log(a[i][j])
							+ b[j][t + 1], bwd[i][t]);
			}
		}

		return bwd;
	}

	/** prints all the parameters of an HMM */
	public void print() {
		DecimalFormat fmt = new DecimalFormat();
		fmt.setMinimumFractionDigits(5);
		fmt.setMaximumFractionDigits(5);

		for (int i = 0; i < numStates; i++)
			System.out.println("pi(" + i + ") = " + fmt.format(pi[i]));
		System.out.println();

		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++)
				System.out.print("a(" + i + "," + j + ") = "
						+ fmt.format(a[i][j]) + "  ");
			System.out.println();
		}

		System.out.println();
		for (int i = 0; i < numStates; i++) {
			for (int k = 0; k < b[0].length; k++)
				System.out.print("b(" + i + "," + k + ") = "
						+ fmt.format(b[i][k]) + "  ");
			System.out.println();
		}

		System.out.println();
		for (int i = 0; i < numStates; i++) {
			for (int k = 0; k < distribution[0].length; k++)
				System.out.print("distribution(" + i + "," + k + ") = "
						+ fmt.format(distribution[i][k]) + "  ");
			System.out.println();
		}
	}

	/** divides two doubles. 0 / 0 = 0! */
	public double divide(double n, double d) {
		if (n == 0)
			return 0;
		else
			return n / d;
	}

	/**
	 * calculation of probability P(X_t = s_i, X_t+1 = s_j | O, m).
	 * 
	 * @param fwd
	 *            the Forward-Variables for o
	 * @param bwd
	 *            the Backward-Variables for o
	 * @return xi (i.e., P)
	 */
	public double[][][] calXi(double[][] fwd, double[][] bwd) {
		int T = b[0].length;
		double[][][] xi = new double[numStates][numStates][T];
		double num;
		double denom;

		for (int t = 0; t <= T - 2; t++) {
			denom = 0;
			for (int k = 0; k < numStates; k++)
				denom += (fwd[k][t] * bwd[k][t]);

			for (int i = 0; i < numStates; i++) {
				for (int j = 0; j < numStates; j++) {
					num = fwd[i][t] * a[i][j] * b[j][t + 1] * bwd[j][t + 1];
					xi[i][j][t] = divide(num, denom);
				}
			}
		}

		denom = 0;

		for (int k = 0; k < numStates; k++)
			denom += (fwd[k][T - 1] * bwd[k][T - 1]);

		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				num = fwd[i][T - 1] * a[i][j];
				xi[i][j][T - 1] = divide(num, denom);
			}
		}

		return xi;
	}

	/**
	 * calculation of probability P(X_t = s_i, X_t+1 = s_j | O, m). Log version
	 * 
	 * @param fwd
	 *            the Forward-Variables for o
	 * @param bwd
	 *            the Backward-Variables for o
	 * @return xi (i.e., P)
	 */
	public double[][][] calXiLog(double[][] fwd, double[][] bwd) {
		int T = b[0].length;
		double[][][] xi = new double[numStates][numStates][T];
		double num;
		double denom;

		for (int t = 0; t <= T - 2; t++) {
			denom = LZERO;
			for (int k = 0; k < numStates; k++)
				denom = LAdd(fwd[k][t] + bwd[k][t], denom);

			for (int i = 0; i < numStates; i++) {
				for (int j = 0; j < numStates; j++) {
					num = fwd[i][t] + Math.log(a[i][j]) + b[j][t + 1]
							+ bwd[j][t + 1];
					xi[i][j][t] = num - denom;
				}
			}
		}

		denom = LZERO;

		for (int k = 0; k < numStates; k++)
			denom = LAdd(fwd[k][T - 1] * bwd[k][T - 1], denom);

		for (int i = 0; i < numStates; i++) {
			for (int j = 0; j < numStates; j++) {
				num = fwd[i][T - 1] + Math.log(a[i][j]);
				xi[i][j][T - 1] = num - denom;
			}
		}

		return xi;
	}

	/** computes gamma(i, t) */
	public double[][] calGamma(double[][] fwd, double[][] bwd) {
		int T = b[0].length;
		double[][] gamma = new double[numStates][T];
		double denom;
		double num;

		for (int t = 0; t < T; t++) {
			denom = 0;
			for (int j = 0; j < numStates; j++)
				denom += fwd[j][t] * bwd[j][t];

			for (int i = 0; i < numStates; i++) {
				num = fwd[i][t] * bwd[i][t];
				gamma[i][t] = divide(num, denom);
			}
		}

		return gamma;
	}

	/** computes gamma(i, t); log version */
	public double[][] calGammaLog(double[][] fwd, double[][] bwd) {
		int T = b[0].length;
		double[][] gamma = new double[numStates][T];
		double denom;
		double num;

		for (int t = 0; t < T; t++) {
			denom = LZERO;
			for (int j = 0; j < numStates; j++)
				denom = LAdd(fwd[j][t] + bwd[j][t], denom);

			for (int i = 0; i < numStates; i++) {
				num = fwd[i][t] + bwd[i][t];
				gamma[i][t] = num - denom;
			}
		}
		return gamma;
	}
}