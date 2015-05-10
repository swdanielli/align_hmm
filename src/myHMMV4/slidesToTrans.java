package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import myHMMV4.TranscriptionClass;

public class slidesToTrans {
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
	 * 
	 * 
	 * @param vocabSize
	 *            size of dictionary
	 * @param verbose
	 *            >= 1: output linking result >= 2: output aligning result
	 * @param adaptVersion
	 *            1: hard decide segment, adapt emission probability 2: soft
	 *            decide segment, adapt emission probability (bug: initializing
	 *            pi1, a1) 3: soft decide segment, adapt emission, transition,
	 *            and initial probability 4: all to one state -> No Segment
	 * @param isUsedHMM
	 *            1: with transition probability 0: without transition
	 *            probability
	 * @param feaSelType
	 *            0: no feature selection 1: selecting feature by keyword list
	 *            (tuning the weight), vocabSize should be a_b, where a is total
	 *            dict size (ordinary + appendix), b is ordinary dict size 2:
	 *            using the first n sentences in each section of textbook
	 * @param segmentType
	 *            0: sentences 1: overlapped window 2: non-overlapped chunks
	 *            with fixed length
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
		int adaptVersion = Integer.parseInt(args[5]);
		if (adaptVersion != 4) {
			adaptVersion = 1;
		}

		int isUseHMM = Integer.parseInt(args[6]);
		int feaSelType = Integer.parseInt(args[7]);
		int cvFold = Integer.parseInt(args[8]);
		int segmentType = Integer.parseInt(args[9]);
		int vocabSize = 0;
		int ordinaryVocabSize = 0;

		String[] vocabs = args[3].split("_");
		vocabSize = Integer.parseInt(vocabs[0]);
		if (feaSelType == 1) {
			ordinaryVocabSize = Integer.parseInt(vocabs[1]);
		} else {
			ordinaryVocabSize = vocabSize;
		}

		int trainingStep = 5;
		BufferedReader br = null;

		// List<String> idList = new ArrayList<String>();
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
				// idList.add(slots[1]);

				SlidesClass slideDistributionSparse = new SlidesClass(args[1]
						+ "/" + slots[0]);
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
			} else {
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

			double smoothing = 0.1;
			for (double slidesTransRatio : new double[] { 0.3, 0.2, 0.1 }) {
				System.out.println("slidesTransRatio: " + slidesTransRatio);
				tuneResult = tuneSegment(tuneResult, smoothing,
						slidesTransRatio, adaptVersion, isUseHMM,
						transObjTrain, slideDistributionSparseTrain, vocabSize,
						trainingStep, ordinaryVocabSize, feaSelType,
						evaluation, segmentType);
			}

			System.out.println("========== Cross Validation Set " + cvIndex
					+ " Result ==========");
			/*
			 * System.out.println(tuneResult.nConfigSet);
			 * System.out.println(tuneResult.bestConfiguration.length);
			 * System.out.println(tuneResult.bestConfiguration[0].length);
			 */
			for (int i = 0; i < tuneResult.nConfigSet; i++) {
				System.out.println("smoothing: "
						+ tuneResult.bestConfiguration[i][0]);
				System.out.println("slidesTransRatio: "
						+ tuneResult.bestConfiguration[i][1]);
				System.out.println("transition: "
						+ tuneResult.bestConfiguration[i][2]);
				System.out.println("keywordWeight: "
						+ tuneResult.bestConfiguration[i][3]);
				System.out.println("segment_length: "
						+ tuneResult.bestConfiguration[i][4]);

				List<TranscriptionClass> transObjArray_seg = new ArrayList<TranscriptionClass>();
				for (int j = 0; j < transObjTest.size(); j++)
					transObjArray_seg.add(transObjTest.get(j).segment(
							segmentType,
							(int) tuneResult.bestConfiguration[i][4]));

				@SuppressWarnings("unchecked")
				double res[] = paramTune(transObjArray_seg,
						slideDistributionSparseTest, vocabSize,
						tuneResult.bestConfiguration[i][0], trainingStep,
						tuneResult.bestConfiguration[i][1], adaptVersion,
						tuneResult.bestConfiguration[i][2], verbose,
						ordinaryVocabSize, tuneResult.bestConfiguration[i][3],
						1).result;
				result[i] += res[i];
				if (verbose >= 1)
					System.out
							.println("========== End Transcription linked result ==========");
			}
		}
		System.out
				.println("========== Averaged Validation Set Result ==========");
		for (int i = 0; i < result.length; i++)
			System.out.println(result[i] / cvFold);
	}

	public static tuneConfiguration tuneSegment(tuneConfiguration tuneResult,
			double smoothing, double slidesTransRatio, int adaptVersion,
			int isUseHMM, List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			int trainingStep, int ordinaryVocabSize, int feaSelType,
			int evaluation, int segmentType) {
		/* length shouldn't be too long such that n_segment < n_states */
		int win_length[] = { 1, 2, 3, 4, 5, 6 }; /* left + right + self */
		int chunk_length[] = { 1, 2, 3, 4, 5, 6 }; /* left + right + self */
		int sent_length[] = { 0 };
		int length[] = null;

		switch (segmentType) {
		case 1:
			length = win_length;
			break;
		case 2:
			length = chunk_length;
			break;
		default:
			length = sent_length;
		}

		for (int seg_len : length) {
			List<TranscriptionClass> transObjArray_seg = new ArrayList<TranscriptionClass>();
			for (int i = 0; i < transObjArray.size(); i++)
				transObjArray_seg.add(transObjArray.get(i).segment(segmentType,
						seg_len));

			tuneResult = tuneFeaSel(tuneResult, smoothing, slidesTransRatio,
					adaptVersion, isUseHMM, transObjArray_seg,
					slideDistributionSparseArray, vocabSize, trainingStep,
					ordinaryVocabSize, feaSelType, evaluation, seg_len);
		}
		return tuneResult;
	}

	public static tuneConfiguration tuneFeaSel(tuneConfiguration tuneResult,
			double smoothing, double slidesTransRatio, int adaptVersion,
			int isUseHMM, List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			int trainingStep, int ordinaryVocabSize, int feaSelType,
			int evaluation, int seg_len) {
		switch (feaSelType) {
		case 1:
			double keywordWeight = 5;
			tuneResult = tuneTransition(tuneResult, smoothing,
					slidesTransRatio, adaptVersion, isUseHMM, transObjArray,
					slideDistributionSparseArray, vocabSize, trainingStep,
					ordinaryVocabSize, keywordWeight, evaluation, seg_len);
			break;
		default:
			tuneResult = tuneTransition(tuneResult, smoothing,
					slidesTransRatio, adaptVersion, isUseHMM, transObjArray,
					slideDistributionSparseArray, vocabSize, trainingStep,
					ordinaryVocabSize, 0.0, evaluation, seg_len);
		}
		return tuneResult;
	}

	@SuppressWarnings("unchecked")
	public static tuneConfiguration tuneTransition(
			tuneConfiguration tuneResult, double smoothing,
			double slidesTransRatio, int adaptVersion, int isUseHMM,
			List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			int trainingStep, int ordinaryVocabSize, double keywordWeight,
			int evaluation, int seg_len) {
		if (isUseHMM == 0) {
			double transition = 0;
			double result[] = paramTune(transObjArray,
					slideDistributionSparseArray, vocabSize, smoothing,
					trainingStep, slidesTransRatio, adaptVersion, transition,
					-1, ordinaryVocabSize, keywordWeight, evaluation).result;
			tuneResult.updateConfiguration(result, new double[] { smoothing,
					slidesTransRatio, transition, keywordWeight, seg_len });
		} else {
			for (double transition : new double[] { 8.1, 2.7, 0.9, 0.3, 0.1,
					0.03 }) {
				double result[] = paramTune(transObjArray,
						slideDistributionSparseArray, vocabSize, smoothing,
						trainingStep, slidesTransRatio, adaptVersion,
						transition, -1, ordinaryVocabSize, keywordWeight,
						evaluation).result;
				tuneResult.updateConfiguration(result, new double[] {
						smoothing, slidesTransRatio, transition, keywordWeight,
						seg_len });
			}
		}

		return tuneResult;
	}

	public static Interpolate_result paramTune(List<TranscriptionClass> transObjArray,
			List<SlidesClass> slideDistributionSparseArray, int vocabSize,
			double smoothing, int trainingStep, double slidesTransRatio,
			int adaptVersion, double transition, int verbose, int ordinaryVocabSize,
			double keywordWeight, int evaluation, List<TranscriptionClass>... ref_trans_objs) {
		List<TranscriptionClass> interpolated_trans = null;
		boolean is_interpolation = false;
		if (ref_trans_objs.length == 2) {
			is_interpolation = true;
			interpolated_trans = new ArrayList<TranscriptionClass>();
		}
		/* Training set */
		if (verbose >= 1)
			System.out.println("========== Start trans seg result ==========");

		double[] res = { 0 };
		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			double accSum = 0.0;
			double sentNum = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);

				TranscriptionClass obs_trans_obj = null;
				TranscriptionClass label_trans_obj = null;
				if (is_interpolation) {
					obs_trans_obj = ref_trans_objs[0].get(i);
					label_trans_obj = ref_trans_objs[1].get(i);
				}
				Interpolate_result hmm_result = alignSlides(transObjArray.get(i),
						slideDistributionSparseArray.get(i), vocabSize,
						smoothing, trainingStep, slidesTransRatio,
						adaptVersion, transition, verbose, ordinaryVocabSize,
						keywordWeight, evaluation, obs_trans_obj, label_trans_obj);

				accSum += hmm_result.result[0];
				sentNum += hmm_result.result[1];
				if (is_interpolation) interpolated_trans.add(hmm_result.transObjs.get(0));
			}

			if (verbose >= 0) {
				System.out.println("Number of Sents: " + sentNum);
				System.out.println("acc: " + accSum / sentNum);
			}
			res[0] = accSum / sentNum;
		} else {
			double ll = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);

				TranscriptionClass obs_trans_obj = null;
				TranscriptionClass label_trans_obj = null;
				if (is_interpolation) {
					obs_trans_obj = ref_trans_objs[0].get(i);
					label_trans_obj = ref_trans_objs[1].get(i);
				}
				Interpolate_result hmm_result = alignSlides(transObjArray.get(i),
						slideDistributionSparseArray.get(i), vocabSize,
						smoothing, trainingStep, slidesTransRatio,
						adaptVersion, transition, verbose, ordinaryVocabSize,
						keywordWeight, evaluation, obs_trans_obj, label_trans_obj);

				ll += hmm_result.result[0];
				if (is_interpolation) interpolated_trans.add(hmm_result.transObjs.get(0));
			}

			if (verbose >= 0) {
				System.out.println("Likelihood: " + ll);
			}
			res[0] = ll;
		}

		if (verbose >= 1)
			System.out.println("========== End trans seg result ==========");
		return new Interpolate_result(interpolated_trans, res);
	}

	public static Interpolate_result alignSlides(TranscriptionClass transObj, 
			SlidesClass slideDistributionSparse, int vocabSize, double smoothing,
			int trainingStep, double slidesTransRatio, int adaptVersion, double transition,
			int verbose, int ordinaryVocabSize, double keywordWeight, int evaluation,
			TranscriptionClass... ref_trans_obj) {
		List<TranscriptionClass> interpolated_trans = null;
		boolean is_interpolation = false;
		if (ref_trans_obj.length == 2 && ref_trans_obj[0] != null && ref_trans_obj[1] != null) {
			is_interpolation = true;
			interpolated_trans = new ArrayList<TranscriptionClass>();
		}
		
		int segWordsMax = 0; /* use all transcription for training */
		HMM hmmModel = new HMM(slideDistributionSparse.observationId.length,
				vocabSize, ordinaryVocabSize, keywordWeight, segWordsMax,
				smoothing);
		/* To tail ?? */

		/* Compute model */
		// double normalizationTerm = 0.0;
		// int idxGap = 0;

		/* transition probability */
		double[][] a = new double[slideDistributionSparse.observationId.length][];
		for (int i = 0; i < a.length; i++)
			a[i] = new double[slideDistributionSparse.observationId.length];

		if (transition == 0) {
			for (int i = 0; i < a.length; i++) {
				for (int j = 0; j < a[i].length; j++) {
					a[i][j] = 1.0 / a[i].length;
				}
			}
		} else {
			for (int i = 0; i < a.length; i++) {
				for (int j = 0; j < a[i].length; j++) {
					if (i == a.length - 1) {
						if (i == j) {
							a[i][j] = 1;
						} else {
							a[i][j] = 0;
						}
					} else {
						if (i == j) {
							a[i][j] = 1 - transition;
						} else if ((i + 1) == j) {
							a[i][j] = transition;
						} else {
							a[i][j] = 0;
						}
					}
				}
			}
		}

		hmmModel.setA(a);

		/*
		 * for (int i = 0; i < a.length; i++) { normalizationTerm = 0.0; for
		 * (int j = 0; j < a[i].length; j++) { idxGap = i - j; if (i >= j) {
		 * a[i][j] = Math.exp(-transitionDecayBackw * idxGap); } else { a[i][j]
		 * = Math.exp(transitionDecayForw * idxGap); } normalizationTerm +=
		 * a[i][j]; } for (int j = 0; j < a[i].length; j++) a[i][j] /=
		 * normalizationTerm; }
		 */

		/* initial probability */
		double[] pi = new double[slideDistributionSparse.observationId.length];
		pi[0] = 1;
		for (int i = 1; i < pi.length; i++)
			pi[i] = 0;
		hmmModel.setPi(pi);

		hmmModel.setDistributionSparse(slideDistributionSparse.observationId,
				slideDistributionSparse.observationCount);
		/*
		 * hmmModel.setDistributionSparseWindowing(slideDistributionSparse.
		 * observationId, slideDistributionSparse.observationCount,
		 * slidesWindowing);
		 */

		switch (adaptVersion) {
		case 1:
			hmmModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, slidesTransRatio);
			break;
		case 2:
			hmmModel.trainAdaptationLog(transObj.observationId,
					transObj.observationCount, trainingStep, slidesTransRatio,
					false);
			break;
		case 3:
			hmmModel.trainAdaptationLog(transObj.observationId,
					transObj.observationCount, trainingStep, slidesTransRatio,
					true);
			break;
		case 4:
			hmmModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, slidesTransRatio);
			break;
		}

		Map<Double, int[]> resultMap = hmmModel.viterbiLog(
				transObj.observationId, transObj.observationCount);
		double[] score = { 0 };
		int[] stateSeq = null;
		for (Double key : resultMap.keySet()) {
			score[0] = key;
			stateSeq = resultMap.get(key);
		}

		if (verbose >= 1) {
			double[][] alignedPosterior = hmmModel.computeStatePosterior(
					transObj.observationId, transObj.observationCount);

			String[] stateSeqS = new String[stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = Integer.toString(stateSeq[i] + 1); /* IMPORTANT */
			transObj.printResult(stateSeqS, alignedPosterior);
		}
		
		if (is_interpolation) {
			interpolated_trans.add(transObj.interpolate(stateSeq,
					slideDistributionSparse.observationId,
					slideDistributionSparse.observationCount, ref_trans_obj[0],
					ref_trans_obj[1]));
		}

		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			String[] stateSeqS = new String[stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = Integer.toString(stateSeq[i] + 1); /* IMPORTANT */
			return new Interpolate_result(interpolated_trans, transObj.acc(stateSeqS));
		} else {
			return new Interpolate_result(interpolated_trans, score);
		}
	}
}
