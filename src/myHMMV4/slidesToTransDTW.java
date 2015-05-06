package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class slidesToTransDTW {
	/**
	 * run alignment; future, training.
	 * 
	 * @param transcription
	 *            (directory) o[t], sentence or segment based, line split,
	 *            sparse
	 * @param slides
	 *            (directory) state, distribution; format, 1: wordId... (page
	 *            title) 2: wordId... (bullet) 2: wordId... 1: wordId... (next
	 *            page)
	 * @param traininglist
	 *            training file list and mapping to ID


	 * @param vocabSize
	 *            size of dictionary
	 * @param verbose
	 *            >= 1: output linking result 
	 *            >= 2: output aligning result
	 * @param adaptVersion
	 *            1: hard decide segment, adapt emission probability 
	 *            2: soft decide segment, adapt emission probability (bug: initializing pi1, a1)
	 *            3: soft decide segment, adapt emission, transition, and initial probability
	 *            4: all to one state -> No Segment
	 * @param isUsedHMM
	 *            1: with transition probability
	 *            0: without transition probability
	 * @param feaSelType
	 *            0: no feature selection
	 *            1: selecting feature by keyword list (tuning the weight), vocabSize should be a_b, 
	 *               where a is total dict size (ordinary + appendix), b is ordinary dict size
	 *            2: using the first n sentences in each section of textbook
	 * @param smoothing
	 *            language model smoothing in computing emission probability
	 * @param trainingStep
	 *            step of training
	 * @param slidesTransRatio
	 *            ratio between slides and transcription when updating emission
	 *            probability
	 * @param transitionDecayForw
	 *            coefficient of (forward) exponential decay in transition
	 *            probability
	 * @param transitionDecayBackw
	 *            coefficient of (backward) exponential decay in transition
	 *            probability
	 * @throws IOException
	 */
	
	public static void main(String[] args) throws IOException {
		int verbose = Integer.parseInt(args[4]);
		//int adaptVersion = Integer.parseInt(args[5]);
		//int isUseHMM = Integer.parseInt(args[6]);
		int feaSelType = Integer.parseInt(args[5]);
		int cvFold = Integer.parseInt(args[6]);
		int vocabSize = 0;
		int ordinaryVocabSize = 0;
		
		String[] vocabs = args[3].split("_");
		vocabSize = Integer.parseInt(vocabs[0]);
		if (feaSelType == 1) {			
			ordinaryVocabSize = Integer.parseInt(vocabs[1]);
		}
		else {						
			ordinaryVocabSize = vocabSize;
		}

		int trainingStep = 5;
		BufferedReader br = null;

		//List<String> idList = new ArrayList<String>();
		List<TranscriptionClass> transObjArray = new ArrayList<TranscriptionClass>();
		List<SlidesClass> slideDistributionSparseArray = new ArrayList<SlidesClass>();
		// CV set * chapter

		/* Load Transcription and slides */
		try {
			br = new BufferedReader(new FileReader(args[2]));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			String line = br.readLine();
			while (line != null) {
				String[] slots = line.split(".pdf\t");
				if (verbose >= 2)
					System.out.println(slots[1]);
				//idList.add(slots[1]);
				
				SlidesClass slideDistributionSparse = new SlidesClass(args[1] + "/"
						+ slots[0]);				
				TranscriptionClass transObj = new TranscriptionClass(args[0]
						+ "/" + slots[0]);
				
				transObjArray.add(transObj);				
				slideDistributionSparseArray.add(slideDistributionSparse);
				line = br.readLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		double[] result = { 0.0 };
		
		for (int cvIndex = 0; cvIndex < cvFold; cvIndex++) {
			tuneConfiguration tuneResult;
			tuneResult = new tuneConfiguration(1);

			List<TranscriptionClass> transObjTrain = new ArrayList<TranscriptionClass>();
			List<TranscriptionClass> transObjTest = new ArrayList<TranscriptionClass>();
			List<SlidesClass> slideDistributionSparseTrain = new ArrayList<SlidesClass>();
			List<SlidesClass> slideDistributionSparseTest = new ArrayList<SlidesClass>();
			/* evaluation -> 1 accuracy; 0 likelihood */
			int evaluation = 1;
			if (cvFold == 1) { 
				evaluation = 0; /* don't tune parameter on train? */
				transObjTrain = transObjArray;
				transObjTest = transObjArray;
				slideDistributionSparseTrain = slideDistributionSparseArray;
				slideDistributionSparseTest = slideDistributionSparseArray;
			}
			else {
				for (int i = 0; i < transObjArray.size(); i++) {
					if ((i % cvFold) == cvIndex) {
						transObjTest.add(transObjArray.get(i));
						slideDistributionSparseTest
								.add(slideDistributionSparseArray.get(i));
					} else {
						transObjTrain.add(transObjArray.get(i));
						slideDistributionSparseTrain
								.add(slideDistributionSparseArray.get(i));
					}
				}
			}
			
			for (double smoothing : new double[] {10, 1, 0.1, 0.01}) { /* original 10, 1, 0.1, 0.01, 0.001 */
				for (double slidesTransRatio : new double[] { 0.9, 0.7, 0.5, 0.3, 0.1 }) {
					System.out.println("smoothing: " + smoothing
							+ " slidesTransRatio: " + slidesTransRatio);
					tuneResult = tuneFeaSel(tuneResult, smoothing,
							slidesTransRatio, transObjTrain,
							slideDistributionSparseTrain, vocabSize,
							trainingStep, ordinaryVocabSize, feaSelType,
							evaluation);
				}
			}

			System.out.println("========== Cross Validation Set " + cvIndex
					+ " Result ==========");
			/*System.out.println(tuneResult.nConfigSet);
			System.out.println(tuneResult.bestConfiguration.length);
			System.out.println(tuneResult.bestConfiguration[0].length);*/
			for (int i = 0; i < tuneResult.nConfigSet; i++) {
				System.out.println("smoothing: " + tuneResult.bestConfiguration[i][0]);
				System.out.println("slidesTransRatio: " + tuneResult.bestConfiguration[i][1]);
				System.out.println("keywordWeight: " + tuneResult.bestConfiguration[i][2]);				
				System.out.println("slidesWindowing: " + tuneResult.bestConfiguration[i][3]);
				System.out.println("dtwWindowSize: " + tuneResult.bestConfiguration[i][4]);
				
				double res[] = paramTune(transObjTest,
						slideDistributionSparseTest, vocabSize,
						tuneResult.bestConfiguration[i][0], trainingStep,
						tuneResult.bestConfiguration[i][1], verbose,
						ordinaryVocabSize, tuneResult.bestConfiguration[i][2],
						1, tuneResult.bestConfiguration[i][3],
						(int) tuneResult.bestConfiguration[i][4]);
				result[i] += res[i];
				if (verbose >= 1)
                    System.out.println("========== End Transcription linked result ==========");
			}
		}
		System.out
				.println("========== Averaged Validation Set Result ==========");
		for (int i = 0; i < result.length; i++)
			System.out.println(result[i] / cvFold);
	}

	public static tuneConfiguration tuneFeaSel(tuneConfiguration tuneResult,
			double smoothing, double slidesTransRatio,
			List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			int trainingStep, int ordinaryVocabSize, int feaSelType,
			int evaluation) {
		switch (feaSelType) {
		case 1:
			for (double keywordWeight : new double[] { 10, 5, 2, 1, 0.5, 0.2,
					0.1 })
				tuneResult = tuneTransition(tuneResult, smoothing,
						slidesTransRatio, transObjArray,
						slideDistributionSparseArray, vocabSize, trainingStep,
						ordinaryVocabSize, keywordWeight, evaluation);
			break;
		default:
			tuneResult = tuneTransition(tuneResult, smoothing,
					slidesTransRatio, transObjArray,
					slideDistributionSparseArray, vocabSize, trainingStep,
					ordinaryVocabSize, 0.0, evaluation);
		}
		return tuneResult;
	}

	public static tuneConfiguration tuneTransition(
			tuneConfiguration tuneResult, double smoothing,
			double slidesTransRatio, List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			int trainingStep, int ordinaryVocabSize, double keywordWeight,
			int evaluation) {
		for (double slidesWindowing : new double[] { 10, 7, 5, 3, 1, 0.1 }) {
			for (int dtwWindowSize : new int[] { 150, 125, 100, 90, 80, 70, 60, 50}) {
				double result[] = paramTune(transObjArray,
						slideDistributionSparseArray, vocabSize, smoothing,
						trainingStep, slidesTransRatio, -1, ordinaryVocabSize,
						keywordWeight, evaluation, slidesWindowing,
						dtwWindowSize);
				tuneResult.updateConfiguration(result, new double[] {
						smoothing, slidesTransRatio, keywordWeight,
						slidesWindowing, dtwWindowSize });
			}
		}
		return tuneResult;
	}
	
	public static double[] paramTune(List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			double smoothing, int trainingStep, double slidesTransRatio,
			int verbose, int ordinaryVocabSize, double keywordWeight,
			int evaluation, double slidesWindowing, int dtwWindowSize) {
		/* Training set */
		if (verbose >= 1)
            System.out.println("========== Start trans seg result ==========");
		
		double [] res = {0};
		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			double accSum = 0.0;
			double sentNum = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);
				
				double result[] = null;
				result = alignSlides(transObjArray.get(i),
						slideDistributionSparseArray.get(i), vocabSize,
						smoothing, trainingStep, slidesTransRatio, verbose,
						ordinaryVocabSize, keywordWeight, evaluation,
						slidesWindowing, dtwWindowSize);
				
				accSum += result[0];
				sentNum += result[1];
			}
			
			if (verbose >= 0) {
				System.out.println("Number of Sents: " + sentNum);
				System.out.println("acc: " + accSum / sentNum);
			}
			res[0] = accSum / sentNum;
		}
		else {
			double ll = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);
				
				double result[] = null;
				result = alignSlides(transObjArray.get(i),
						slideDistributionSparseArray.get(i), vocabSize,
						smoothing, trainingStep, slidesTransRatio, verbose,
						ordinaryVocabSize, keywordWeight, evaluation,
						slidesWindowing, dtwWindowSize);
				
				ll += result[0];
			}
			
			if (verbose >= 0) {
				System.out.println("Likelihood: " + ll);
			}
			res[0] = ll;			
		}
		
		if (verbose >= 1)
            System.out.println("========== End trans seg result ==========");		
		return res;
	}
	
	public static double[] alignSlides(TranscriptionClass transObj,
			SlidesClass slideDistributionSparse, int vocabSize,
			double smoothing, int trainingStep, double slidesTransRatio,
			int verbose, int ordinaryVocabSize, double keywordWeight,
			int evaluation, double slidesWindowing, int dtwWindowSize) {
		int segWordsMax = 0; /* use all transcription for training */
		DTW dtwModel = new DTW(slideDistributionSparse.observationId.length, vocabSize,
				ordinaryVocabSize, keywordWeight, segWordsMax, smoothing);

		/* word distribution from slides */
		dtwModel.setDistributionSparseWindowing(slideDistributionSparse.observationId,
				slideDistributionSparse.observationCount, slidesWindowing);

		dtwModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, slidesTransRatio, dtwWindowSize);
		
		Map<Double, int[]> resultMap = dtwModel.viterbiLog(transObj.observationId,
				transObj.observationCount, dtwWindowSize);
		double [] score = {0};
		int[] stateSeq = null;
		for (Double key: resultMap.keySet()) {
			score[0] = key;
			stateSeq = resultMap.get(key);
		}
		if (verbose >= 1) {
			double [][] alignedPosterior = dtwModel.computeStatePosterior(stateSeq);			
			
			String[] stateSeqS = new String [stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = Integer.toString(stateSeq[i] + 1);  /* IMPORTANT */
			transObj.printResult(stateSeqS, alignedPosterior);
		}
		
		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			String[] stateSeqS = new String [stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = Integer.toString(stateSeq[i] + 1);  /* IMPORTANT */
			return transObj.acc(stateSeqS);
		}
		else { return score; }
	}
}