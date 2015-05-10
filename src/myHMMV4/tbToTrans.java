package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import myHMMV4.TranscriptionClass;

public class tbToTrans {
	/**
	 * run alignment; future, training.
	 * 
	 * @param transcription
	 *            (directory) o[t], sentence or segment based, line split,
	 *            sparse
	 * @param tbList
	 * @param transList
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
	 * @param cvFold
	 *            # cv
	 * @param isUseSlides
	 *            interpolated with slides (1) or not (0)
	 * @param decayStyle
	 *            (\d)_(\d) $1: gamma decay style: from the beginning of chapter
	 *            (0), from the beginning of the first chapter all (1) $2:
	 *            transition decay style: section (0), chapter (1)
	 * 
	 * 
	 * 
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

		int isUseSlides = Integer.parseInt(args[9]);
		int segmentType = Integer.parseInt(args[10]);

		String gammaDecayStyle = args[11];

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

		/* Chapter * section */
		List<TextbookClass> textbookArray = new ArrayList<TextbookClass>();

		/* Load Textbook */
		try {
			br = new BufferedReader(new FileReader(args[1]));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			String line = br.readLine();
			while (line != null) {
				textbookArray.add(new TextbookClass(line));
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

		/* Load Transcription and slides */
		br = null;
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

				TranscriptionClass transObj = null;
				if (isUseSlides == 1) {
					transObj = new TranscriptionClass(args[0] + "/" + slots[0],
							args[0].replace("trans", "slides") + "/" + slots[0]);
				} else {
					transObj = new TranscriptionClass(args[0] + "/" + slots[0]);
				}

				transObjArray.add(transObj);
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

			List<TextbookClass> textbookModels = textbookArray;
			/* evaluation -> 1 accuracy; 0 likelihood */
			int evaluation = 1;
			if (cvFold == 1) {
				evaluation = 0; /* don't tune parameter on train? */
				transObjTrain = transObjArray;
				transObjTest = transObjArray;
			} else {
				for (int i = 0; i < transObjArray.size(); i++) {
					if ((i % cvFold) == cvIndex) {
						transObjTest.add(transObjArray.get(i));
					} else {
						transObjTrain.add(transObjArray.get(i));
					}
				}
			}

			double wSlides[] = { 0.3, 0.2, 0.1 };
			double woSlides[] = { 0.0 };
			double slidesTransRatios[] = null;
			if (isUseSlides == 1) {
				slidesTransRatios = wSlides;
			} else {
				slidesTransRatios = woSlides;
			}

			double smoothing = 0.1;
			String tDecayStyle = "1";
			double tbTransRatio = 0.1;
			for (double slidesTransRatio : slidesTransRatios) {
				System.out.println("slidesTransRatio: " + slidesTransRatio);
				tuneResult = tuneSegment(tuneResult, smoothing,
						slidesTransRatio, tbTransRatio, adaptVersion, isUseHMM,
						transObjTrain, textbookModels, vocabSize, trainingStep,
						ordinaryVocabSize, feaSelType, evaluation,
						gammaDecayStyle + tDecayStyle, segmentType);
			}

			System.out.println("========== Cross Validation Set " + cvIndex
					+ " Result ==========");
			/*
			 * System.out.println(tuneResult.nConfigSet);
			 * System.out.println(tuneResult.bestConfiguration.length);
			 * System.out.println(tuneResult.bestConfiguration[0].length);
			 */
			for (int i = 0; i < tuneResult.nConfigSet; i++) {
				System.out.println("smoothing: " + tuneResult.bestConfiguration[i][0]);
				System.out.println("slidesTransRatio: " + tuneResult.bestConfiguration[i][1]);
				System.out.println("transitionDecayBackw: " + tuneResult.bestConfiguration[i][2]);
				System.out.println("transitionDecayForw: " + tuneResult.bestConfiguration[i][3]);
				System.out.println("initDecay: " + tuneResult.bestConfiguration[i][4]);
				System.out.println("keywordWeight: " + tuneResult.bestConfiguration[i][5]);
				System.out.println("decayStyle: "
								+ Integer.toString((int) tuneResult.bestConfiguration[i][6])
								+ "_"
								+ Integer.toString((int) tuneResult.bestConfiguration[i][7]));
				System.out.println("segment_length: "
						+ tuneResult.bestConfiguration[i][8]);

				List<TranscriptionClass> transObjArray_seg = new ArrayList<TranscriptionClass>();
				for (int j = 0; j < transObjTest.size(); j++)
					transObjArray_seg.add(transObjTest.get(j).segment(
							segmentType,
							(int) tuneResult.bestConfiguration[i][8]));

				@SuppressWarnings("unchecked")
				double res[] = paramTune(
						transObjArray_seg, textbookModels, vocabSize, tuneResult.bestConfiguration[i][0], 
						trainingStep, tuneResult.bestConfiguration[i][1], tbTransRatio, adaptVersion,
						tuneResult.bestConfiguration[i][2],	tuneResult.bestConfiguration[i][3],
						tuneResult.bestConfiguration[i][4], verbose, ordinaryVocabSize,
						tuneResult.bestConfiguration[i][5],	1,
						Integer.toString((int) tuneResult.bestConfiguration[i][6])
								+ "_"
								+ Integer.toString((int) tuneResult.bestConfiguration[i][7])).result;

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
			double smoothing, double slidesTransRatio, double tbTransRatio,
			int adaptVersion, int isUseHMM,
			List<TranscriptionClass> transObjArray,
			List<TextbookClass> textbookArray, int vocabSize, int trainingStep,
			int ordinaryVocabSize, int feaSelType, int evaluation,
			String decayStyle, int segmentType) {
		/* length shouldn't be too long such that n_segment < n_states */
		int win_length[] = { 1, 2, 4, 6 }; /* left + right + self */
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
					tbTransRatio, adaptVersion, isUseHMM, transObjArray_seg,
					textbookArray, vocabSize, trainingStep, ordinaryVocabSize,
					feaSelType, evaluation, decayStyle, seg_len);
		}
		return tuneResult;
	}

	public static tuneConfiguration tuneFeaSel(tuneConfiguration tuneResult,
			double smoothing, double slidesTransRatio, double tbTransRatio,
			int adaptVersion, int isUseHMM,
			List<TranscriptionClass> transObjArray,
			List<TextbookClass> textbookArray, int vocabSize, int trainingStep,
			int ordinaryVocabSize, int feaSelType, int evaluation,
			String decayStyle, int seg_len) {
		switch (feaSelType) {
		case 1:
			double keywordWeight = 5;
			tuneResult = tuneTransition(tuneResult, smoothing,
					slidesTransRatio, tbTransRatio, adaptVersion, isUseHMM,
					transObjArray, textbookArray, vocabSize, trainingStep,
					ordinaryVocabSize, keywordWeight, evaluation, decayStyle,
					seg_len);
			break;
		default:
			tuneResult = tuneTransition(tuneResult, smoothing,
					slidesTransRatio, tbTransRatio, adaptVersion, isUseHMM,
					transObjArray, textbookArray, vocabSize, trainingStep,
					ordinaryVocabSize, 0.0, evaluation, decayStyle, seg_len);
		}
		return tuneResult;
	}

	@SuppressWarnings("unchecked")
	public static tuneConfiguration tuneTransition(
			tuneConfiguration tuneResult, double smoothing,
			double slidesTransRatio, double tbTransRatio, double adaptVersion,
			int isUseHMM, List<TranscriptionClass> transObjArray,
			List<TextbookClass> textbookArray, int vocabSize, int trainingStep,
			int ordinaryVocabSize, double keywordWeight, int evaluation,
			String decayStyle, int seg_len) {
		if (isUseHMM == 0) {
			double transitionDecayForw = 0;
			double transitionDecayBackw = 0;
			double initDecay = 0;

			double result[] = paramTune(transObjArray, textbookArray,
					vocabSize, smoothing, trainingStep, slidesTransRatio,
					tbTransRatio, (int) adaptVersion, transitionDecayBackw,
					transitionDecayForw, initDecay, -1, ordinaryVocabSize,
					keywordWeight, evaluation, decayStyle).result;
			String tags[] = decayStyle.split("_");
			tuneResult.updateConfiguration(
					result,
					new double[] { smoothing, slidesTransRatio,
							transitionDecayBackw, transitionDecayForw,
							initDecay, keywordWeight,
							Float.parseFloat(tags[0]),
							Float.parseFloat(tags[1]), seg_len });
		} else {
			for (double transitionDecayForw : new double[] { 5, 2, 1, 0.3, 0.1 }) {
				for (double transitionDecayBackw : new double[] { 5, 2, 1, 0.3, 0.1 }) {
					for (double initDecay : new double[] { 0.1, 0.03, 0.01 }) {
						double result[] = paramTune(transObjArray,
								textbookArray, vocabSize, smoothing,
								trainingStep, slidesTransRatio, tbTransRatio,
								(int) adaptVersion, transitionDecayBackw,
								transitionDecayForw, initDecay, -1,
								ordinaryVocabSize, keywordWeight, evaluation,
								decayStyle).result;
						String tags[] = decayStyle.split("_");
						tuneResult.updateConfiguration(
								result,
								new double[] { smoothing, slidesTransRatio,
										transitionDecayBackw,
										transitionDecayForw, initDecay,
										keywordWeight,
										Float.parseFloat(tags[0]),
										Float.parseFloat(tags[1]), seg_len });
					}
				}
			}
		}

		return tuneResult;
	}

	public static Interpolate_result paramTune(List<TranscriptionClass> transObjArray,
			List<TextbookClass> textbookArray, int vocabSize, double smoothing,
			int trainingStep, double slidesTransRatio, double tbTransRatio,
			int adaptVersion, double transitionDecayBackw,
			double transitionDecayForw, double initDecay, int verbose,
			int ordinaryVocabSize, double keywordWeight, int evaluation,
			String decayStyle, List<TranscriptionClass>... ref_trans_objs) {
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
				Interpolate_result hmm_result = alignTB(transObjArray.get(i), textbookArray, 
						vocabSize, smoothing, trainingStep, tbTransRatio, adaptVersion,
						transitionDecayBackw, transitionDecayForw, initDecay,
						verbose, ordinaryVocabSize, keywordWeight, evaluation,
						decayStyle, obs_trans_obj, label_trans_obj);

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
				Interpolate_result hmm_result = alignTB(transObjArray.get(i), textbookArray, 
						vocabSize, smoothing, trainingStep, tbTransRatio, adaptVersion,
						transitionDecayBackw, transitionDecayForw, initDecay,
						verbose, ordinaryVocabSize, keywordWeight, evaluation,
						decayStyle, obs_trans_obj, label_trans_obj);

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

	public static Interpolate_result alignTB(TranscriptionClass transObj,
			List<TextbookClass> textbookArray, int vocabSize, double smoothing,
			int trainingStep, double tbTransRatio, int adaptVersion,
			double transitionDecayBackw, double transitionDecayForw,
			double initDecay, int verbose, int ordinaryVocabSize,
			double keywordWeight, int evaluation, String decayStyle, 
			TranscriptionClass... ref_trans_obj) {
		List<TranscriptionClass> interpolated_trans = null;
		boolean is_interpolation = false;
		if (ref_trans_obj.length == 2 && ref_trans_obj[0] != null && ref_trans_obj[1] != null) {
			is_interpolation = true;
			interpolated_trans = new ArrayList<TranscriptionClass>();
		}
		
		int segWordsMax = 0; /* use all transcription for training */
		String tags[] = decayStyle.split("_");
		int initDecayVersion = Integer.parseInt(tags[0]);
		int transitionDecayVersion = Integer.parseInt(tags[1]);

		List<int[]> arrayArrayBuffer = new ArrayList<int[]>();
		List<double[]> arrayArrayCntBuffer = new ArrayList<double[]>();
		List<String> sectionIdBuffer = new ArrayList<String>();
		List<Integer> chapterList = new ArrayList<Integer>();

		for (int i = 0; i < textbookArray.size(); i++) {
			for (int j = 0; j < textbookArray.get(i).observationId.length; j++) {
				arrayArrayBuffer.add(textbookArray.get(i).observationId[j]);
				arrayArrayCntBuffer
						.add(textbookArray.get(i).observationCount[j]);
				sectionIdBuffer.add(textbookArray.get(i).sectionId[j]);
				chapterList.add(i);
			}
		}

		int observationId[][] = new int[arrayArrayBuffer.size()][];
		for (int i = 0; i < arrayArrayBuffer.size(); i++)
			observationId[i] = arrayArrayBuffer.get(i);

		double observationCount[][] = new double[arrayArrayCntBuffer.size()][];
		for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
			observationCount[i] = arrayArrayCntBuffer.get(i);

		/* Compute model */
		HMM hmmModel = new HMM(observationId.length, vocabSize,
				ordinaryVocabSize, keywordWeight, segWordsMax, smoothing);

		/* Emission probability */
		hmmModel.setDistributionSparse(observationId, observationCount);
		// hmmModel.setDistributionSparseWindowing(observationId,
		// observationCount, slidesWindowing);

		/* transition probability */
		double normalizationTerm = 0.0;
		int idxGap = 0;

		double[][] a = new double[observationId.length][];
		for (int i = 0; i < a.length; i++)
			a[i] = new double[observationId.length];
		for (int i = 0; i < a.length; i++) {
			normalizationTerm = 0.0;
			for (int j = 0; j < a[i].length; j++) {
				if (transitionDecayVersion == 0) {
					idxGap = i - j;
				} else if (transitionDecayVersion == 1) {
					idxGap = chapterList.get(i) - chapterList.get(j);
					if (i > j)
						idxGap += 1;
					else if (i < j)
						idxGap -= 1;
				}

				if (i >= j) {
					a[i][j] = Math.exp(-transitionDecayBackw * idxGap);
				} else {
					a[i][j] = Math.exp(transitionDecayForw * idxGap);
				}
				normalizationTerm += a[i][j];
			}
			for (int j = 0; j < a[i].length; j++)
				a[i][j] /= normalizationTerm;
		}
		hmmModel.setA(a);

		/* initial probability */
		int idxStart = 0;
		double[] pi = new double[observationId.length];
		normalizationTerm = 0.0;
		for (int i = 0; i < pi.length; i++) {
			if (initDecayVersion == 1) {
				idxStart = 0;
			} else if (initDecayVersion == 0) {
				if (i != 0 && chapterList.get(i) != chapterList.get(i - 1))
					idxStart = i;
			}

			idxGap = i - idxStart;

			pi[i] = Math.exp(-initDecay * idxGap);
			normalizationTerm += pi[i];
		}
		for (int i = 0; i < pi.length; i++)
			pi[i] /= normalizationTerm;
		hmmModel.setPi(pi);

		switch (adaptVersion) {
		case 1:
			hmmModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, tbTransRatio);
			break;
		case 2:
			hmmModel.trainAdaptationLog(transObj.observationId,
					transObj.observationCount, trainingStep, tbTransRatio,
					false);
			break;
		case 3:
			hmmModel.trainAdaptationLog(transObj.observationId,
					transObj.observationCount, trainingStep, tbTransRatio, true);
			break;
		case 4:
			hmmModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, tbTransRatio);
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
				stateSeqS[i] = sectionIdBuffer.get(stateSeq[i]); /* IMPORTANT */
			transObj.printResult(stateSeqS, alignedPosterior);
		}

		if (is_interpolation) {
			interpolated_trans.add(transObj.interpolate(stateSeq,
					observationId, observationCount, ref_trans_obj[0],
					ref_trans_obj[1]));
		}

		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			String[] stateSeqS = new String[stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = sectionIdBuffer.get(stateSeq[i]); /* IMPORTANT */
			return new Interpolate_result(interpolated_trans, transObj.acc(stateSeqS));
		} else {
			return new Interpolate_result(interpolated_trans, score);
		}
	}
}
