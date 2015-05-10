package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import myHMMV4.TranscriptionClass;

public class tb_slides_2_trans {
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
		int verbose = Integer.parseInt(args[6]);
		/********************** TODO *************************/
		verbose = 0;
		/*****************************************************/
		int adaptVersion = Integer.parseInt(args[7]);
		if (adaptVersion != 4) {
			adaptVersion = 1;
		}
		int isUseHMM = Integer.parseInt(args[8]);
		int feaSelType = Integer.parseInt(args[9]);
		int cvFold = Integer.parseInt(args[10]);
		int segmentType = Integer.parseInt(args[11]);
		int n_co_train_iter = Integer.parseInt(args[12]);
		String gammaDecayStyle = args[13];
		
		int vocabSize = 0;
		int ordinaryVocabSize = 0;

		String[] vocabs = args[5].split("_");
		vocabSize = Integer.parseInt(vocabs[0]);
		if (feaSelType == 1) {
			ordinaryVocabSize = Integer.parseInt(vocabs[1]);
		} else {
			ordinaryVocabSize = vocabSize;
		}

		int trainingStep = 5;
		BufferedReader br = null;

		/* List of lectures */
		List<TranscriptionClass> transObjs_label_slides = new ArrayList<TranscriptionClass>();
		List<SlidesClass> slideDistributionSparseArray = new ArrayList<SlidesClass>();

		List<TranscriptionClass> transObjs_label_tb = new ArrayList<TranscriptionClass>();
		/* List of chapters */
		List<TextbookClass> textbookArray = new ArrayList<TextbookClass>();

		/* Load Textbook */
		try {
			br = new BufferedReader(new FileReader(args[3]));
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
			br = new BufferedReader(new FileReader(args[4]));
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

				SlidesClass slideDistributionSparse = new SlidesClass(args[2]
						+ "/" + slots[0]);
				TranscriptionClass transObj_label_slides = new TranscriptionClass(args[1]
						+ "/" + slots[0]);
				TranscriptionClass transObj_label_tb = new TranscriptionClass(args[0]
						+ "/" + slots[0]);

				transObjs_label_slides.add(transObj_label_slides);
				transObjs_label_tb.add(transObj_label_tb);
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

		double[] result_slides = { 0.0 };
		double[] result_tb = { 0.0 };

		for (int cvIndex = 0; cvIndex < cvFold; cvIndex++) {
			List<tuneConfiguration> tuneResults_slides = new ArrayList<tuneConfiguration>();
			List<tuneConfiguration> tuneResults_tb = new ArrayList<tuneConfiguration>();

			List<TranscriptionClass> transObjs_train_slides = new ArrayList<TranscriptionClass>();
			List<TranscriptionClass> transObjs_test_slides = new ArrayList<TranscriptionClass>();
			List<SlidesClass> slideDistributionSparseTrain = new ArrayList<SlidesClass>();
			List<SlidesClass> slideDistributionSparseTest = new ArrayList<SlidesClass>();

			List<TranscriptionClass> transObjs_train_tb = new ArrayList<TranscriptionClass>();
			List<TranscriptionClass> transObjs_test_tb = new ArrayList<TranscriptionClass>();
			List<TextbookClass> textbookModels = textbookArray;
			
			/* evaluation -> 1 accuracy; 0 likelihood */
			int evaluation = 1;
			if (cvFold == 1) {
				evaluation = 0; /* don't tune parameter on train? */
				transObjs_train_slides = transObjs_label_slides;
				transObjs_test_slides = transObjs_label_slides;
				
				transObjs_train_tb = transObjs_label_tb;
				transObjs_test_tb = transObjs_label_tb;

				slideDistributionSparseTrain = slideDistributionSparseArray;
				slideDistributionSparseTest = slideDistributionSparseArray;
			} else {
				for (int i = 0; i < transObjs_label_slides.size(); i++) {
					if ((i % cvFold) == cvIndex) {
						transObjs_test_slides.add(transObjs_label_slides.get(i));
						transObjs_test_tb.add(transObjs_label_tb.get(i));
						slideDistributionSparseTest
								.add(slideDistributionSparseArray.get(i));
					} else {
						transObjs_train_slides.add(transObjs_label_slides.get(i));
						transObjs_train_tb.add(transObjs_label_tb.get(i));
						slideDistributionSparseTrain
								.add(slideDistributionSparseArray.get(i));
					}
				}
			}

			double smoothing = 0.1;
			String tDecayStyle = "1";
			double tbTransRatio = 0.1;
			double interpolation_weight = 0.1;

			for (int co_train_iter = 0; co_train_iter < n_co_train_iter; co_train_iter++) {
				tuneConfiguration tuneResult_slides = new tuneConfiguration(1);
				tuneConfiguration tuneResult_tb = new tuneConfiguration(1);

				for (double slidesTransRatio : new double[] { 0.3, 0.2, 0.1 }) {
					System.out.println("slidesTransRatio: " + slidesTransRatio);
					tuneResult_slides = slidesToTrans.tuneSegment(
							tuneResult_slides, smoothing, slidesTransRatio,
							adaptVersion, isUseHMM, transObjs_train_slides,
							slideDistributionSparseTrain, vocabSize,
							trainingStep, ordinaryVocabSize, feaSelType,
							evaluation, segmentType);
				}
				transObjs_train_tb = interpolate_slides_to_trans(tuneResult_slides,
						transObjs_train_slides, slideDistributionSparseTrain,
						interpolation_weight, transObjs_train_tb, adaptVersion,
						vocabSize, trainingStep, -1, 0, ordinaryVocabSize,
						feaSelType, evaluation, segmentType).transObjs;

				tuneResult_tb = tbToTrans.tuneSegment(tuneResult_tb, smoothing,
						0.0, tbTransRatio, adaptVersion, isUseHMM,
						transObjs_train_tb, textbookModels, vocabSize, trainingStep,
						ordinaryVocabSize, feaSelType, evaluation,
						gammaDecayStyle + tDecayStyle, segmentType);
				
				tuneResults_slides.add(tuneResult_slides);
				tuneResults_tb.add(tuneResult_tb);
			}
			
			System.out.println("========== Cross Validation Set " + cvIndex
					+ " Result ==========");
			for (int i = 0; i < tuneResults_slides.get(0).nConfigSet; i++) {
				for (int co_train_iter = 0; co_train_iter < n_co_train_iter-1; co_train_iter++) {
					System.out.println("========== Co-train iter: " + (co_train_iter+1) + " ==========");
					System.out.println("smoothing: " + tuneResults_slides.get(co_train_iter).bestConfiguration[i][0]);
					System.out.println("slidesTransRatio: "	+ tuneResults_slides.get(co_train_iter).bestConfiguration[i][1]);
					System.out.println("transition: " + tuneResults_slides.get(co_train_iter).bestConfiguration[i][2]);
					System.out.println("keywordWeight: " + tuneResults_slides.get(co_train_iter).bestConfiguration[i][3]);
					System.out.println("segment_length: " + tuneResults_slides.get(co_train_iter).bestConfiguration[i][4]);
					
					
					System.out.println("smoothing: " + tuneResults_tb.get(co_train_iter).bestConfiguration[i][0]);
					System.out.println("slidesTransRatio: "	+ tuneResults_tb.get(co_train_iter).bestConfiguration[i][1]);
					System.out.println("transitionDecayBackw: "	+ tuneResults_tb.get(co_train_iter).bestConfiguration[i][2]);
					System.out.println("transitionDecayForw: " + tuneResults_tb.get(co_train_iter).bestConfiguration[i][3]);
					System.out.println("initDecay: " + tuneResults_tb.get(co_train_iter).bestConfiguration[i][4]);
					System.out.println("keywordWeight: " + tuneResults_tb.get(co_train_iter).bestConfiguration[i][5]);
					System.out.println("decayStyle: "
									+ Integer.toString((int) tuneResults_tb.get(co_train_iter).bestConfiguration[i][6])
									+ "_"
									+ Integer.toString((int) tuneResults_tb.get(co_train_iter).bestConfiguration[i][7]));
					System.out.println("segment_length: " + tuneResults_tb.get(co_train_iter).bestConfiguration[i][8]);

					transObjs_test_tb = interpolate_slides_to_trans(
							tuneResults_slides.get(co_train_iter), transObjs_test_slides,
							slideDistributionSparseTest, interpolation_weight, transObjs_test_tb,
							adaptVersion, vocabSize, trainingStep, -1, i, ordinaryVocabSize,
							feaSelType, evaluation, segmentType).transObjs;
				}

				System.out.println("========== Co-train iter: " + (n_co_train_iter) + " ==========");
				System.out.println("smoothing: " + tuneResults_slides.get(n_co_train_iter-1).bestConfiguration[i][0]);
				System.out.println("slidesTransRatio: "	+ tuneResults_slides.get(n_co_train_iter-1).bestConfiguration[i][1]);
				System.out.println("transition: " + tuneResults_slides.get(n_co_train_iter-1).bestConfiguration[i][2]);
				System.out.println("keywordWeight: " + tuneResults_slides.get(n_co_train_iter-1).bestConfiguration[i][3]);
				System.out.println("segment_length: " + tuneResults_slides.get(n_co_train_iter-1).bestConfiguration[i][4]);
				
				
				System.out.println("smoothing: " + tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][0]);
				System.out.println("slidesTransRatio: "	+ tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][1]);
				System.out.println("transitionDecayBackw: "	+ tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][2]);
				System.out.println("transitionDecayForw: " + tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][3]);
				System.out.println("initDecay: " + tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][4]);
				System.out.println("keywordWeight: " + tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][5]);
				System.out.println("decayStyle: "
								+ Integer.toString((int) tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][6])
								+ "_"
								+ Integer.toString((int) tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][7]));
				System.out.println("segment_length: " + tuneResults_tb.get(n_co_train_iter-1).bestConfiguration[i][8]);

				result_slides[i] += interpolate_slides_to_trans(
						tuneResults_slides.get(n_co_train_iter-1),
						transObjs_test_slides, slideDistributionSparseTest,
						interpolation_weight, transObjs_test_tb, adaptVersion,
						vocabSize, trainingStep, verbose, i, ordinaryVocabSize,
						feaSelType, evaluation, segmentType).result[i];
				result_tb[i] += interpolate_tb_to_trans(
						tuneResults_tb.get(n_co_train_iter-1),
						transObjs_test_tb, textbookModels, interpolation_weight,
						transObjs_test_slides, tbTransRatio, adaptVersion, vocabSize,
						trainingStep, verbose, i, ordinaryVocabSize, feaSelType,
						evaluation,	segmentType).result[i];

				if (verbose >= 1)
					System.out
							.println("========== End Transcription linked result ==========");
			}
		}
		System.out
				.println("========== Averaged Validation Set Result ==========");
		for (int i = 0; i < result_slides.length; i++)
			System.out.println(result_slides[i] / cvFold);
		for (int i = 0; i < result_tb.length; i++)
			System.out.println(result_tb[i] / cvFold);
	}
	
	@SuppressWarnings("unchecked")
	public static Interpolate_result interpolate_slides_to_trans(
			tuneConfiguration tuneResult, List<TranscriptionClass> transObjs_slides,
			List<SlidesClass> slideDistributionSparse, double interpolation_weight, 
			List<TranscriptionClass> transObjs_tb, int adaptVersion, int vocabSize, 
			int trainingStep, int verbose, int config_set_id,
			int ordinaryVocabSize, int feaSelType, int evaluation,
			int segmentType) {
		List<TranscriptionClass> transObjs_seg = new ArrayList<TranscriptionClass>();
		for (int j = 0; j < transObjs_slides.size(); j++)
			transObjs_seg.add(transObjs_slides.get(j).segment(
					segmentType,
					(int) tuneResult.bestConfiguration[config_set_id][4]));

		return slidesToTrans.paramTune(
				transObjs_seg, slideDistributionSparse, vocabSize,
				tuneResult.bestConfiguration[config_set_id][0], trainingStep,
				tuneResult.bestConfiguration[config_set_id][1], adaptVersion,
				tuneResult.bestConfiguration[config_set_id][2], verbose,
				ordinaryVocabSize, tuneResult.bestConfiguration[config_set_id][3],
				evaluation, transObjs_slides, transObjs_tb);
	}
	
	@SuppressWarnings("unchecked")
	public static Interpolate_result interpolate_tb_to_trans(
			tuneConfiguration tuneResult, List<TranscriptionClass> transObjs_tb,
			List<TextbookClass> textbookModels, double interpolation_weight,
			List<TranscriptionClass> transObjs_slides, double tbTransRatio,
			int adaptVersion, int vocabSize, int trainingStep, int verbose,
			int config_set_id, int ordinaryVocabSize, int feaSelType,
			int evaluation,	int segmentType) {
		List<TranscriptionClass> transObjs_seg = new ArrayList<TranscriptionClass>();
		for (int j = 0; j < transObjs_tb.size(); j++)
			transObjs_seg.add(transObjs_tb.get(j).segment(
					segmentType,
					(int) tuneResult.bestConfiguration[config_set_id][8]));

		return tbToTrans.paramTune(
				transObjs_seg, textbookModels,
				vocabSize, tuneResult.bestConfiguration[config_set_id][0],
				trainingStep, tuneResult.bestConfiguration[config_set_id][1],
				tbTransRatio, adaptVersion,
				tuneResult.bestConfiguration[config_set_id][2],
				tuneResult.bestConfiguration[config_set_id][3],
				tuneResult.bestConfiguration[config_set_id][4],
				verbose, ordinaryVocabSize,
				tuneResult.bestConfiguration[config_set_id][5],	evaluation,
				Integer.toString((int) tuneResult.bestConfiguration[config_set_id][6])
						+ "_"
						+ Integer.toString((int) tuneResult.bestConfiguration[config_set_id][7]),
				transObjs_tb, transObjs_slides);
	}
}
