package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class TranscriptionClass {
	public List<ArrayList<String>> truth = new ArrayList<ArrayList<String>>();
	public int observationId[][];
	public double observationCount[][];
	public int observationSlidesId[][];
	public double observationSlidesCount[][];
	public String courseId;

	public TranscriptionClass(List<ArrayList<String>> truth,
			int observationId[][], double observationCount[][], String courseId) {
		this.truth = truth;
		this.observationId = observationId;
		this.observationCount = observationCount;
		this.courseId = courseId;
	}

	public TranscriptionClass(List<ArrayList<String>> truth,
			int observationId[][], double observationCount[][],
			int observationSlidesId[][], double observationSlidesCount[][],
			String courseId) {
		this.truth = truth;
		this.observationId = observationId;
		this.observationCount = observationCount;
		this.observationSlidesId = observationSlidesId;
		this.observationSlidesCount = observationSlidesCount;
		this.courseId = courseId;
	}

	public TranscriptionClass(String transcriptionFileName, boolean... params) {
		assert params.length <= 1;
		boolean isAggregated = params.length > 0 ? params[0] : false;

		// System.out.println(transcriptionFileName);
		List<int[]> arrayArrayBuffer = null;
		List<double[]> arrayArrayCntBuffer = null;

		String temps[] = transcriptionFileName.split("/");
		courseId = temps[temps.length - 1];

		BufferedReader br = null;
		/* Load transcription, o */
		try {
			br = new BufferedReader(new FileReader(transcriptionFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			String line = br.readLine();
			String sentBuffer = "";

			List<String> truthBuffer = new ArrayList<String>();

			arrayArrayBuffer = new ArrayList<int[]>();
			arrayArrayCntBuffer = new ArrayList<double[]>();

			int breakPointNum = 0;
			while (line != null) {
				if (line.contains("++++++++++")) {
					if (!isAggregated) {
						if (breakPointNum > 0) {
							truth.add((ArrayList<String>) truthBuffer);

							sentBuffer = sentBuffer.trim();
							int[] arrayBuffer = null;
							double[] arrayCntBuffer = null;

							if (sentBuffer.length() > 0) {
								String wordId[] = sentBuffer.split("\\s+");
								arrayBuffer = new int[wordId.length];
								arrayCntBuffer = new double[wordId.length];

								for (int i = 0; i < wordId.length; i++) {
									String slots[] = wordId[i].split(":");
									arrayBuffer[i] = Integer.parseInt(slots[0]);
									arrayCntBuffer[i] = Float
											.parseFloat(slots[1]);
								}
							} else {
								arrayBuffer = new int[0];
								arrayCntBuffer = new double[0];
							}

							arrayArrayBuffer.add(arrayBuffer);
							arrayArrayCntBuffer.add(arrayCntBuffer);

							sentBuffer = "";
						}
						breakPointNum++;
						truthBuffer = new ArrayList<String>();
					}
				} else if (line.contains("==========")) {
					line = br.readLine();
					continue;
				} else {
					String chunks[] = line.split("\t");

					/* change format, label separated with _ */
					if (sentBuffer.length() > 0) {
						sentBuffer += " " + chunks[0].trim();
					} else {
						sentBuffer = chunks[0].trim();
					}

					if (chunks.length > 1) {
						truthBuffer.add(chunks[1]);
					}
				}
				line = br.readLine();
			}

			if (isAggregated) {
				truth.add((ArrayList<String>) truthBuffer);
				sentBuffer = sentBuffer.trim();

				int[] arrayBuffer = null;
				double[] arrayCntBuffer = null;

				if (sentBuffer.length() > 0) {
					String wordId[] = sentBuffer.split("\\s+");
					arrayBuffer = new int[wordId.length];
					arrayCntBuffer = new double[wordId.length];

					for (int i = 0; i < wordId.length; i++) {
						String slots[] = wordId[i].split(":");
						arrayBuffer[i] = Integer.parseInt(slots[0]);
						arrayCntBuffer[i] = Float.parseFloat(slots[1]);
					}
				} else {
					arrayBuffer = new int[0];
					arrayCntBuffer = new double[0];
				}

				arrayArrayBuffer.add(arrayBuffer);
				arrayArrayCntBuffer.add(arrayCntBuffer);

				sentBuffer = "";
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

		observationId = new int[arrayArrayBuffer.size()][];
		for (int i = 0; i < arrayArrayBuffer.size(); i++)
			observationId[i] = arrayArrayBuffer.get(i);

		observationCount = new double[arrayArrayCntBuffer.size()][];
		for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
			observationCount[i] = arrayArrayCntBuffer.get(i);
	}

	public TranscriptionClass(String transcriptionFileName,
			String slidesFileName, boolean... params) {
		assert params.length <= 1;
		boolean isAggregated = params.length > 0 ? params[0] : false;

		// System.out.println(transcriptionFileName);
		List<int[]> arrayArrayBuffer = null;
		List<double[]> arrayArrayCntBuffer = null;
		List<int[]> arrayArrayBufferSlides = null;
		List<double[]> arrayArrayCntBufferSlides = null;

		String temps[] = transcriptionFileName.split("/");
		courseId = temps[temps.length - 1];

		BufferedReader br = null;
		BufferedReader brSlides = null;

		/* Load transcription, o */
		try {
			br = new BufferedReader(new FileReader(transcriptionFileName));
			brSlides = new BufferedReader(new FileReader(slidesFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			String line = br.readLine();
			String slidesBuffer = "";
			String sentBuffer = "";

			List<String> truthBuffer = new ArrayList<String>();

			arrayArrayBuffer = new ArrayList<int[]>();
			arrayArrayCntBuffer = new ArrayList<double[]>();
			arrayArrayBufferSlides = new ArrayList<int[]>();
			arrayArrayCntBufferSlides = new ArrayList<double[]>();

			int breakPointNum = 0;
			while (line != null) {
				if (line.contains("++++++++++")) {
					if (!isAggregated) {
						if (breakPointNum > 0) {
							truth.add((ArrayList<String>) truthBuffer);

							sentBuffer = sentBuffer.trim();
							slidesBuffer = slidesBuffer.trim();
							int[] arrayBuffer = null;
							double[] arrayCntBuffer = null;
							int[] arrayBufferSlides = null;
							double[] arrayCntBufferSlides = null;

							if (sentBuffer.length() > 0) {
								String wordId[] = sentBuffer.split("\\s+");
								arrayBuffer = new int[wordId.length];
								arrayCntBuffer = new double[wordId.length];

								for (int i = 0; i < wordId.length; i++) {
									String slots[] = wordId[i].split(":");
									arrayBuffer[i] = Integer.parseInt(slots[0]);
									arrayCntBuffer[i] = Float
											.parseFloat(slots[1]);
								}
							} else {
								arrayBuffer = new int[0];
								arrayCntBuffer = new double[0];
							}

							if (slidesBuffer.length() > 0) {
								String wordId[] = slidesBuffer.split("\\s+");
								arrayBufferSlides = new int[wordId.length];
								arrayCntBufferSlides = new double[wordId.length];

								for (int i = 0; i < wordId.length; i++) {
									String slots[] = wordId[i].split(":");
									arrayBufferSlides[i] = Integer
											.parseInt(slots[0]);
									arrayCntBufferSlides[i] = Float
											.parseFloat(slots[1])
											/ truthBuffer.size();
								}
							} else {
								arrayBufferSlides = new int[0];
								arrayCntBufferSlides = new double[0];
							}

							arrayArrayBuffer.add(arrayBuffer);
							arrayArrayCntBuffer.add(arrayCntBuffer);
							arrayArrayBufferSlides.add(arrayBufferSlides);
							arrayArrayCntBufferSlides.add(arrayCntBufferSlides);

							sentBuffer = "";
							slidesBuffer = "";
						}
						breakPointNum++;
						truthBuffer = new ArrayList<String>();
					}
				} else if (line.contains("==========")) {
					line = br.readLine();
					continue;
				} else {
					String chunks[] = line.split("\t"); /* NEW Separation tag */

					if (sentBuffer.length() > 0) {
						slidesBuffer += " " + brSlides.readLine().trim();
						sentBuffer += " " + chunks[0].trim();
					} else {
						slidesBuffer = brSlides.readLine().trim();
						sentBuffer = chunks[0].trim();
					}

					if (chunks.length > 1) {
						truthBuffer.add(chunks[1]);
					}
				}
				line = br.readLine();
			}

			if (isAggregated) {
				truth.add((ArrayList<String>) truthBuffer);
				sentBuffer = sentBuffer.trim();
				slidesBuffer = slidesBuffer.trim();
				int[] arrayBuffer = null;
				double[] arrayCntBuffer = null;
				int[] arrayBufferSlides = null;
				double[] arrayCntBufferSlides = null;

				if (sentBuffer.length() > 0) {
					String wordId[] = sentBuffer.split("\\s+");
					arrayBuffer = new int[wordId.length];
					arrayCntBuffer = new double[wordId.length];

					for (int i = 0; i < wordId.length; i++) {
						String slots[] = wordId[i].split(":");
						arrayBuffer[i] = Integer.parseInt(slots[0]);
						arrayCntBuffer[i] = Float.parseFloat(slots[1]);
					}
				} else {
					arrayBuffer = new int[0];
					arrayCntBuffer = new double[0];
				}

				if (slidesBuffer.length() > 0) {
					String wordId[] = slidesBuffer.split("\\s+");
					arrayBufferSlides = new int[wordId.length];
					arrayCntBufferSlides = new double[wordId.length];

					for (int i = 0; i < wordId.length; i++) {
						String slots[] = wordId[i].split(":");
						arrayBufferSlides[i] = Integer.parseInt(slots[0]);
						arrayCntBufferSlides[i] = Float.parseFloat(slots[1])
								/ truthBuffer.size();
					}
				} else {
					arrayBufferSlides = new int[0];
					arrayCntBufferSlides = new double[0];
				}

				arrayArrayBuffer.add(arrayBuffer);
				arrayArrayCntBuffer.add(arrayCntBuffer);
				arrayArrayBufferSlides.add(arrayBufferSlides);
				arrayArrayCntBufferSlides.add(arrayCntBufferSlides);

				sentBuffer = "";
				slidesBuffer = "";
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

		observationId = new int[arrayArrayBuffer.size()][];
		for (int i = 0; i < arrayArrayBuffer.size(); i++)
			observationId[i] = arrayArrayBuffer.get(i);

		observationCount = new double[arrayArrayCntBuffer.size()][];
		for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
			observationCount[i] = arrayArrayCntBuffer.get(i);

		observationSlidesId = new int[arrayArrayBufferSlides.size()][];
		for (int i = 0; i < arrayArrayBufferSlides.size(); i++)
			observationSlidesId[i] = arrayArrayBufferSlides.get(i);

		observationSlidesCount = new double[arrayArrayCntBufferSlides.size()][];
		for (int i = 0; i < arrayArrayCntBufferSlides.size(); i++)
			observationSlidesCount[i] = arrayArrayCntBufferSlides.get(i);
	}

	public TranscriptionClass compute_interpolation(double ratio) {
		int newObservationId[][] = new int[observationId.length][];
		double newObservationCount[][] = new double[observationCount.length][];

		for (int i = 0; i < observationId.length; i++) {
			HashMap<Integer, Double> count = new HashMap<Integer, Double>();
			for (int j = 0; j < observationId[i].length; j++) {
				if (count.containsKey(observationId[i][j])) {
					count.put(observationId[i][j],
							count.get(observationId[i][j])
									+ observationCount[i][j]
									* (1 - ratio));
				} else {
					count.put(observationId[i][j], observationCount[i][j]
							* (1 - ratio));
				}
			}

			for (int j = 0; j < observationSlidesId[i].length; j++) {
				if (count.containsKey(observationSlidesId[i][j])) {
					count.put(observationSlidesId[i][j],
							count.get(observationSlidesId[i][j])
									+ observationSlidesCount[i][j]
									* ratio);
				} else {
					count.put(observationSlidesId[i][j],
							observationSlidesCount[i][j] * ratio);
				}
			}

			int[] arrayBuffer = new int[count.size()];
			double[] arrayCntBuffer = new double[count.size()];
			int counter = 0;
			for (int key : count.keySet()) {
				arrayBuffer[counter] = key;
				arrayCntBuffer[counter] = count.get(key);
				counter++;
			}
			newObservationId[i] = arrayBuffer;
			newObservationCount[i] = arrayCntBuffer;
		}

		return new TranscriptionClass(truth, newObservationId,
				newObservationCount, courseId);
	}

	/**
	 * compute accuracy
	 * 
	 * @param hypo
	 *            hypothesis
	 */
	/*
	 * public double [] acc (String [] hypo) { double numAcc = 0.0; double
	 * numSent = 0.0;
	 * 
	 * for (int i = 0; i < hypo.length; i++) { for (int j = 0; j <
	 * truth.get(i).size() ; j++) { String wordId []=
	 * truth.get(i).get(j).split(" "); for (int k = 0; k < wordId.length; k++) {
	 * if (wordId[k].equals("0")) { continue; } if (hypo[i].equals(wordId[k])) {
	 * numAcc += 1; } numSent += 1; } } }
	 * 
	 * double[] result = {numAcc, numSent}; return result; }
	 */

	public double[] acc(String[] hypo) {
		double numAcc = 0.0;
		double numSent = 0.0;

		for (int i = 0; i < hypo.length; i++) {
			for (int j = 0; j < truth.get(i).size(); j++) {
				String wordId[] = truth.get(i).get(j).split(" ");
				double truthClass = 0.0;
				double hitClass = 0.0;
				for (int k = 0; k < wordId.length; k++) {
					if (wordId[k].equals("0")) {
						continue;
					}
					if (hypo[i].equals(wordId[k])) {
						hitClass += 1;
					}
					truthClass += 1;
				}
				if (truthClass != 0.0) {
					numAcc += hitClass / truthClass;
					numSent += 1;
				}
			}
		}

		double[] result = { numAcc, numSent };
		return result;
	}
	
	public TranscriptionClass interpolate(int[] stateSeq, int[][] oId,
			double[][] oCount, TranscriptionClass obs_trans_obj, TranscriptionClass label_trans_obj) {
		List<ArrayList<Pair>> pairs_list = new ArrayList<ArrayList<Pair>>();
		int top_n = 49; /* Pick the 50 most frequent words */

		for (int state_idx = 0; state_idx < oId.length; state_idx++) {
			Pair [] pairs_raw = new Pair[oId[state_idx].length];
			for (int word_idx = 0; word_idx < oId[state_idx].length; word_idx++) {
				pairs_raw[word_idx] = new Pair(oId[state_idx][word_idx], oCount[state_idx][word_idx]);
			}
			Arrays.sort(pairs_raw);

			List<Pair> pairs = new ArrayList<Pair>();
			double freq_th = 0.0;
			for (int word_idx = 0; word_idx < pairs_raw.length; word_idx++) {
				if (word_idx == top_n) freq_th = pairs_raw[word_idx].value;
				if (pairs_raw[word_idx].value < freq_th) break;
				pairs.add(pairs_raw[word_idx]);
			}
			pairs_list.add((ArrayList<Pair>)  pairs);
		}

		int interpolate_oId[][] = new int[label_trans_obj.truth.size()][];
		double interpolate_oCount[][] = new double[label_trans_obj.truth.size()][];
		int sent_idx = 0;
		for (int chunk_idx = 0; chunk_idx < stateSeq.length; chunk_idx++) {
			int n_words = pairs_list.get(stateSeq[chunk_idx]).size();
			int[] oId_buffer = new int[n_words];
			double[] oCount_buffer = new double[n_words];
			for (int word_idx = 0; word_idx < n_words; word_idx++) {
				oId_buffer[word_idx] = pairs_list.get(stateSeq[chunk_idx]).get(word_idx).id;
				oCount_buffer[word_idx] = pairs_list.get(stateSeq[chunk_idx]).get(word_idx).value;
			}
			
			for (int i = 0; i < truth.get(chunk_idx).size(); i++) {
				interpolate_oId[sent_idx] = new int[n_words];
				interpolate_oCount[sent_idx] = new double[n_words];
	        	System.arraycopy(oId_buffer, 0, interpolate_oId[sent_idx], 0, n_words);
	        	System.arraycopy(oCount_buffer, 0, interpolate_oCount[sent_idx], 0, n_words);
				sent_idx++;
			}
		}
 
		assert(label_trans_obj.truth.size() == obs_trans_obj.observationId.length);
		assert(label_trans_obj.truth.size() == sent_idx);
		assert(stateSeq.length == observationId.length);

		double ratio = 0.1;
		return new TranscriptionClass(label_trans_obj.truth,
				obs_trans_obj.observationId, obs_trans_obj.observationCount,
				interpolate_oId, interpolate_oCount, courseId)
				.compute_interpolation(ratio);
	}

	public TranscriptionClass segment(int segmentType, int segment_len) {
		if (segmentType == 0) {
			return new TranscriptionClass(truth, observationId,
					observationCount, observationSlidesId,
					observationSlidesCount, courseId);
		}

		int n_segment = observationId.length;
		if (segmentType == 2) n_segment = (int) Math.ceil((double) observationId.length/(2*segment_len+1));

		int newObservationId[][] = new int[n_segment][];
		int newObservationSlidesId[][] = null;
		double newObservationCount[][] = new double[n_segment][];
		double newObservationSlidesCount[][] = null;
		if (observationSlidesId != null && observationSlidesId.length > 0) {
			newObservationSlidesId = new int[n_segment][];
			newObservationSlidesCount = new double[n_segment][];
		}
		
		List<ArrayList<String>> new_truth = null;
		if (segmentType == 1) new_truth = truth;
		else if (segmentType == 2) new_truth = new ArrayList<ArrayList<String>>();

		int counter_segs = 0;
		int idx = 0;
		if (segmentType == 2) idx = segment_len;
		while(true) {
			HashMap<Integer, Double> count = new HashMap<Integer, Double>();
			HashMap<Integer, Double> count_slides = new HashMap<Integer, Double>();
			List<String> truthBuffer = new ArrayList<String>();
			for (int i = idx - segment_len; i <= idx + segment_len; i++) {
				if (i < 0) continue;
				else if (i >= observationId.length) break;

				if (segmentType == 2) {
					for (int j = 0; j < truth.get(i).size(); j++) truthBuffer.add(truth.get(i).get(j));
				}

				for (int j = 0; j < observationId[i].length; j++) {
					if (count.containsKey(observationId[i][j])) {
						count.put(observationId[i][j],
								count.get(observationId[i][j]) + observationCount[i][j]);
					} else {
						count.put(observationId[i][j], observationCount[i][j]);
					}
				}

				if (observationSlidesId != null && observationSlidesId.length > 0) {
					for (int j = 0; j < observationSlidesId[i].length; j++) {
						if (count_slides.containsKey(observationSlidesId[i][j])) {
							count_slides.put(observationSlidesId[i][j],
									count_slides.get(observationSlidesId[i][j]) + observationSlidesCount[i][j]);
						} else {
							count_slides.put(observationSlidesId[i][j], observationSlidesCount[i][j]);
						}
					}
				}
			}

			if (segmentType == 2) new_truth.add((ArrayList<String>) truthBuffer);

			int[] arrayBuffer = new int[count.size()];
			double[] arrayCntBuffer = new double[count.size()];
			int counter = 0;
			for (int key : count.keySet()) {
				arrayBuffer[counter] = key;
				arrayCntBuffer[counter] = count.get(key);
				counter++;
			}
			newObservationId[counter_segs] = arrayBuffer;
			newObservationCount[counter_segs] = arrayCntBuffer;

			if (observationSlidesId != null && observationSlidesId.length > 0) {
				int[] arraySlidesBuffer = new int[count_slides.size()];
				double[] arraySlidesCntBuffer = new double[count_slides.size()];
				int counter_slides = 0;
				for (int key : count_slides.keySet()) {
					arraySlidesBuffer[counter_slides] = key;
					arraySlidesCntBuffer[counter_slides] = count_slides.get(key);
					counter_slides++;
				}
				newObservationSlidesId[counter_segs] = arraySlidesBuffer;
				newObservationSlidesCount[counter_segs] = arraySlidesCntBuffer;
			}

			counter_segs += 1;
			if (segmentType == 1) {
				idx += 1;
				if (idx >= observationId.length) break;
			}
			else if (segmentType == 2) {
				idx += 2 * segment_len + 1;
				if ((idx - segment_len) >= observationId.length) break;
			}
		}

		if (counter_segs != n_segment) System.out.println("Error in TranscriptionClass Segment");
		
		return new TranscriptionClass(new_truth, newObservationId,
				newObservationCount, newObservationSlidesId,
				newObservationSlidesCount, courseId);
	}

	public void printResult(String[] hypo, double[][] alignedPosterior) {
		int sentCnt = 0;
		System.out.println("course : " + courseId);
		for (int i = 0; i < hypo.length; i++) {
			for (int j = 0; j < truth.get(i).size(); j++) {
				sentCnt++;
				System.out.print("sent id: " + sentCnt + ", slides: " + hypo[i]
						+ ", posterior:");
				for (int k = 0; k < alignedPosterior.length; k++) {
					System.out.print(" " + alignedPosterior[k][i]);
				}
				System.out.println("");
			}
		}
	}
}
