package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import myHMMV4.TranscriptionClass;

public class tbToTransNoParamTuneCV {
	/**
	 * run alignment; future, training.
	 * 
	 * @param transcription
	 *            (directory) o[t], sentence or segment based, line split,
	 *            sparse
	 * @param tbList
	 * @param transList
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
	 * @param cvFold
	 *            # cv
	 * @param isUseSlides
	 *            interpolated with slides (1) or not (0)
	 *            
	 *            
	 * @param decayStyle
	 *            (\d)_(\d)
	 *            $1: gamma decay style: chapter (0), all (1)
	 *            $2: transition decay style: section (0), chapter (1)

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
		int feaSelType = Integer.parseInt(args[5]);
		int isUseSlides = Integer.parseInt(args[6]);		
		String gammaDecayStyle = args[7];
		
		int vocabSize = 0;
		int ordinaryVocabSize = 0;
		
		String[] vocabs = args[3].split("_");
		vocabSize = Integer.parseInt(vocabs[0]);
		if (feaSelType == 1) { ordinaryVocabSize = Integer.parseInt(vocabs[1]); }
		else { ordinaryVocabSize = vocabSize; }

		int trainingStep = 5;
		BufferedReader br = null;

		//List<String> idList = new ArrayList<String>();
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
			
		List<TextbookClass> textbookModels = textbookArray;
	        List<TranscriptionClass> transObjTrain = new ArrayList<TranscriptionClass>();
         	int evaluation = 1;
		transObjTrain = transObjArray;

		double smoothing = Float.parseFloat(args[8]);
		double slidesTransRatio = Float.parseFloat(args[9]);
		double tbTransRatio = Float.parseFloat(args[10]);
		String tDecayStyle = args[11];
		double keywordWeight = Float.parseFloat(args[12]);
		double adaptVer = Float.parseFloat(args[13]);
		double transitionDecayForw = Float.parseFloat(args[14]);
		double transitionDecayBackw = Float.parseFloat(args[15]);
		double initDecay = Float.parseFloat(args[16]);
		double slidesWindowing = 10;
		System.out.println("smoothing: " + smoothing);
		System.out.println("slidesTransRatio: " + slidesTransRatio);
		System.out.println("tbTransRatio: " + tbTransRatio);
		System.out.println("adaptVersion: " + adaptVer);
		System.out.println("transitionDecayBackw: " + transitionDecayBackw);
		System.out.println("transitionDecayForw: " + transitionDecayForw);
		System.out.println("initDecay: " + initDecay);
		System.out.println("keywordWeight: " + keywordWeight);
		System.out.println("segWordsMax: " + 0);
		System.out.println("slidesWindowing: " + slidesWindowing);
		System.out.println("decayStyle: " + gammaDecayStyle + tDecayStyle);

		paramTune(transObjTrain, textbookModels, vocabSize,
					smoothing, trainingStep, slidesTransRatio, tbTransRatio,
					(int) adaptVer, transitionDecayBackw, transitionDecayForw,
					initDecay, verbose, ordinaryVocabSize, keywordWeight, 0,
					evaluation, slidesWindowing, gammaDecayStyle + tDecayStyle);
		
                        	
	}
	
	public static double[] paramTune(List<TranscriptionClass> transObjArray,
			List<TextbookClass> textbookArray, int vocabSize, double smoothing,
			int trainingStep, double slidesTransRatio, double tbTransRatio,
			int adaptVersion, double transitionDecayBackw,
			double transitionDecayForw, double initDecay, int verbose,
			int ordinaryVocabSize, double keywordWeight, int segWordsMax,
			int evaluation, double slidesWindowing, String decayStyle) {
		/* Training set */
		if (verbose >= 1)
            System.out.println("========== Start trans seg result ==========");
		
		double [] res = {0};
		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			double[] accSum = new double[transObjArray.size()+1];
			double[] sentNum = new double[transObjArray.size()+1];
                        for (int i = 0; i < accSum.length; i++) {
				accSum[i] = 0.0;
				sentNum[i] = 0.0;
			}
			//double accSum = 0.0;
			//double sentNum = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);
				
				TranscriptionClass interpolatedTrans = null;
				if (slidesTransRatio == 0.0)
					interpolatedTrans = transObjArray.get(i);
				else
					interpolatedTrans = transObjArray.get(i)
							.slidesInterpolation(slidesTransRatio);
				
				double result[] = null;
				result = alignTB(interpolatedTrans, textbookArray,
						vocabSize, smoothing, trainingStep, tbTransRatio,
						adaptVersion, transitionDecayBackw,
						transitionDecayForw, initDecay, verbose,
						ordinaryVocabSize, keywordWeight, evaluation,
						slidesWindowing, decayStyle);
				
				accSum[i] = result[0];
				sentNum[i] = result[1];
				accSum[transObjArray.size()] += result[0];
				sentNum[transObjArray.size()] += result[1];
			}
			if (verbose >= 1) {
				for (int i = 0; i < transObjArray.size(); i++) {
					System.out.println("leave course " + i + " out acc: " + (accSum[transObjArray.size()]-accSum[i]) / (sentNum[transObjArray.size()]-sentNum[i]) );
				}
			}
			if (verbose >= 0) {
				System.out.println("Number of Sents: " + sentNum[transObjArray.size()]);
				System.out.println("acc: " + accSum[transObjArray.size()] / sentNum[transObjArray.size()]);
			}
			res[0] = accSum[transObjArray.size()] / sentNum[transObjArray.size()];
		}
		else {
			double ll = 0.0;
			for (int i = 0; i < transObjArray.size(); i++) {
				if (verbose >= 1)
					System.out.println("transcription/slides : " + i);
				
				TranscriptionClass interpolatedTrans = transObjArray.get(i)
						.slidesInterpolation(slidesTransRatio);
				
				double result[] = null;
				result = alignTB(interpolatedTrans, textbookArray,
						vocabSize, smoothing, trainingStep, tbTransRatio,
						adaptVersion, transitionDecayBackw,
						transitionDecayForw, initDecay, verbose,
						ordinaryVocabSize, keywordWeight, evaluation,
						slidesWindowing, decayStyle);
				
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
	
	public static double[] alignTB(TranscriptionClass transObj,
			List<TextbookClass> textbookArray, int vocabSize, double smoothing,
			int trainingStep, double tbTransRatio, int adaptVersion,
			double transitionDecayBackw, double transitionDecayForw,
			double initDecay, int verbose, int ordinaryVocabSize,
			double keywordWeight, int evaluation, double slidesWindowing,
			String decayStyle) {
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
		
	    int observationId [][]= new int[arrayArrayBuffer.size()][];
		for (int i = 0; i < arrayArrayBuffer.size(); i++)
			observationId[i] = arrayArrayBuffer.get(i);

		double observationCount [][] = new double[arrayArrayCntBuffer.size()][];
		for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
			observationCount[i] = arrayArrayCntBuffer.get(i);

		/* Compute model */
		HMM hmmModel = new HMM(observationId.length, vocabSize,
				ordinaryVocabSize, keywordWeight, segWordsMax, smoothing);
		
		/* Emission probability */
		hmmModel.setDistributionSparseWindowing(observationId,
				observationCount, slidesWindowing);
		
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
					if (i > j) idxGap += 1;
					else if (i < j) idxGap -= 1;
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
					transObj.observationCount, trainingStep, tbTransRatio,
					true);
			break;
		case 4:
			hmmModel.trainAdaptationLogHard(transObj.observationId,
					transObj.observationCount, trainingStep, tbTransRatio);
			break;
		}
		
		Map<Double, int[]> resultMap = hmmModel.viterbiLog(transObj.observationId,
				transObj.observationCount);
		double [] score = {0};
		int[] stateSeq = null;
		for (Double key: resultMap.keySet()) {
			score[0] = key;
			stateSeq = resultMap.get(key);
		}
		if (verbose >= 1) {
			double [][] alignedPosterior = hmmModel.computeStatePosterior(transObj.observationId,
					transObj.observationCount);			
			
			String[] stateSeqS = new String [stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = sectionIdBuffer.get(stateSeq[i]);  /* IMPORTANT */
			transObj.printResult(stateSeqS, alignedPosterior);
		}
		
		/* evaluation -> 1 accuracy; 0 likelihood */
		if (evaluation == 1) {
			String[] stateSeqS = new String [stateSeq.length];
			for (int i = 0; i < stateSeq.length; i++)
				stateSeqS[i] = sectionIdBuffer.get(stateSeq[i]);  /* IMPORTANT */
			return transObj.acc(stateSeqS);
		}
		else { return score; }
	}	
}
