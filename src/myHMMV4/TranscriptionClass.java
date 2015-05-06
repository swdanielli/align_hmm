package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class TranscriptionClass {
	public List<ArrayList<String>> truth = new ArrayList<ArrayList<String>>();
	public int observationId [][];
	public double observationCount [][];
	public int observationSlidesId [][];
	public double observationSlidesCount [][];
	public String courseId;
	
	public TranscriptionClass(List<ArrayList<String>> truth,
			int observationId[][], double observationCount[][], String courseId) {
		this.truth = truth;
		this.observationId = observationId;
		this.observationCount = observationCount;
		this.courseId = courseId;
	}
	
	public TranscriptionClass(String transcriptionFileName, boolean... params) {
		assert params.length <= 1;
		boolean isAggregated = params.length > 0 ? params[0] : false;
		
//		System.out.println(transcriptionFileName);
		List<int[]> arrayArrayBuffer = null;
		List<double[]> arrayArrayCntBuffer = null;

		String temps []= transcriptionFileName.split("/");
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
	        
	        arrayArrayBuffer = new ArrayList<int []>();
	        arrayArrayCntBuffer = new ArrayList<double []>();

	        int breakPointNum = 0;
	        while (line != null) {
	        	if (line.contains("++++++++++")) {
	        		if (!isAggregated) {
		        		if (breakPointNum > 0) {
	        				truth.add((ArrayList<String>) truthBuffer);
		        			
				        	sentBuffer = sentBuffer.trim();	        			
		        			int [] arrayBuffer = null;
		        			double [] arrayCntBuffer = null;
		        			
		        			if (sentBuffer.length() > 0) {
		        				String wordId []= sentBuffer.split("\\s+");
		        				arrayBuffer = new int [wordId.length];
		        				arrayCntBuffer = new double [wordId.length];
		        				
		        				for (int i = 0; i < wordId.length; i++) {
		        					String slots [] = wordId[i].split(":");
		        					arrayBuffer[i] = Integer.parseInt(slots[0]);
		        					arrayCntBuffer[i] = Float.parseFloat(slots[1]);
		        				}
		        			}
		        			else {
		        				arrayBuffer = new int [0];
		        				arrayCntBuffer = new double [0];
		        			}
				        	
				        	arrayArrayBuffer.add(arrayBuffer);
				        	arrayArrayCntBuffer.add(arrayCntBuffer);
				        	
				        	sentBuffer = "";
			        	}
		        		breakPointNum++;
		        		truthBuffer = new ArrayList<String>();
	        		}
	        	}
	        	else if (line.contains("==========")) {
	        		line = br.readLine();
	        		continue;
	        	}
	        	else {
	        		String chunks [] = line.split("\t");
	        		
	        		/* change format, label separated with _ */	        		
	        		if (sentBuffer.length() > 0) { sentBuffer += " " + chunks[0].trim(); }
	        		else { sentBuffer = chunks[0].trim(); }
	        		
	        	    if(chunks.length > 1) {
	        	    	truthBuffer.add(chunks[1]);
	        	    }
	        	}
	        	line = br.readLine();
	        }
	        
	        if (isAggregated) {
		        truth.add((ArrayList<String>) truthBuffer);
		        sentBuffer = sentBuffer.trim();
		        
		        int [] arrayBuffer = null;
    			double [] arrayCntBuffer = null;
    			
    			if (sentBuffer.length() > 0) {
    				String wordId []= sentBuffer.split("\\s+");
    				arrayBuffer = new int [wordId.length];
    				arrayCntBuffer = new double [wordId.length];
    				
    				for (int i = 0; i < wordId.length; i++) {
    					String slots [] = wordId[i].split(":");
    					arrayBuffer[i] = Integer.parseInt(slots[0]);
    					arrayCntBuffer[i] = Float.parseFloat(slots[1]);
    				}
    			}
    			else {
    				arrayBuffer = new int [0];
    				arrayCntBuffer = new double [0];
    			}
	        	
	        	arrayArrayBuffer.add(arrayBuffer);
	        	arrayArrayCntBuffer.add(arrayCntBuffer);

	        	sentBuffer = "";
	        }
	    } catch (IOException e) {
			e.printStackTrace();
		} finally {
	        try { br.close(); }
	        catch (IOException e) { e.printStackTrace(); }
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
		
//		System.out.println(transcriptionFileName);
		List<int[]> arrayArrayBuffer = null;
		List<double[]> arrayArrayCntBuffer = null;
		List<int[]> arrayArrayBufferSlides = null;
		List<double[]> arrayArrayCntBufferSlides = null;

		String temps []= transcriptionFileName.split("/");
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
	        
	        arrayArrayBuffer = new ArrayList<int []>();
	        arrayArrayCntBuffer = new ArrayList<double []>();
	        arrayArrayBufferSlides = new ArrayList<int []>();
	        arrayArrayCntBufferSlides = new ArrayList<double []>();

	        int breakPointNum = 0;
	        while (line != null) {
	        	if (line.contains("++++++++++")) {
	        		if (!isAggregated) {
		        		if (breakPointNum > 0) {
	        				truth.add((ArrayList<String>) truthBuffer);
		        			
				        	sentBuffer = sentBuffer.trim();
				        	slidesBuffer = slidesBuffer.trim();
		        			int [] arrayBuffer = null;
		        			double [] arrayCntBuffer = null;
		        			int [] arrayBufferSlides = null;
		        			double [] arrayCntBufferSlides = null;
		        			
		        			if (sentBuffer.length() > 0) {
		        				String wordId []= sentBuffer.split("\\s+");
		        				arrayBuffer = new int [wordId.length];
		        				arrayCntBuffer = new double [wordId.length];
		        				
		        				for (int i = 0; i < wordId.length; i++) {
		        					String slots [] = wordId[i].split(":");
		        					arrayBuffer[i] = Integer.parseInt(slots[0]);
		        					arrayCntBuffer[i] = Float.parseFloat(slots[1]);
		        				}
		        			}
		        			else {
		        				arrayBuffer = new int [0];
		        				arrayCntBuffer = new double [0];
		        			}
				        	
		        			if (slidesBuffer.length() > 0) {
		        				String wordId []= slidesBuffer.split("\\s+");
		        				arrayBufferSlides = new int [wordId.length];
		        				arrayCntBufferSlides = new double [wordId.length];
		        				
		        				for (int i = 0; i < wordId.length; i++) {
		        					String slots [] = wordId[i].split(":");
		        					arrayBufferSlides[i] = Integer.parseInt(slots[0]);
									arrayCntBufferSlides[i] = Float
											.parseFloat(slots[1])
											/ truthBuffer.size();
		        				}
		        			}
		        			else {
		        				arrayBufferSlides = new int [0];
		        				arrayCntBufferSlides = new double [0];
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
	        	}
	        	else if (line.contains("==========")) {
	        		line = br.readLine();
	        		continue;
	        	}
	        	else {
	        		String chunks [] = line.split("\t"); /* NEW Separation tag */
	    	        	        		
	        		if (sentBuffer.length() > 0) {
	        			slidesBuffer += " " + brSlides.readLine().trim();
	        			sentBuffer += " " + chunks[0].trim();
	        		}
	        		else {
	        			slidesBuffer = brSlides.readLine().trim();
	        			sentBuffer = chunks[0].trim();
	        		}
	        		
	        	    if(chunks.length > 1) {
	        	    	truthBuffer.add(chunks[1]);
	        	    }
	        	}
	        	line = br.readLine();
	        }
	        
	        if (isAggregated) {
		        truth.add((ArrayList<String>) truthBuffer);
		        sentBuffer = sentBuffer.trim();
	        	slidesBuffer = slidesBuffer.trim();
    			int [] arrayBuffer = null;
    			double [] arrayCntBuffer = null;
    			int [] arrayBufferSlides = null;
    			double [] arrayCntBufferSlides = null;
    			
    			if (sentBuffer.length() > 0) {
    				String wordId []= sentBuffer.split("\\s+");
    				arrayBuffer = new int [wordId.length];
    				arrayCntBuffer = new double [wordId.length];
    				
    				for (int i = 0; i < wordId.length; i++) {
    					String slots [] = wordId[i].split(":");
    					arrayBuffer[i] = Integer.parseInt(slots[0]);
    					arrayCntBuffer[i] = Float.parseFloat(slots[1]);
    				}
    			}
    			else {
    				arrayBuffer = new int [0];
    				arrayCntBuffer = new double [0];
    			}

    			if (slidesBuffer.length() > 0) {
    				String wordId []= slidesBuffer.split("\\s+");
    				arrayBufferSlides = new int [wordId.length];
    				arrayCntBufferSlides = new double [wordId.length];
    				
    				for (int i = 0; i < wordId.length; i++) {
    					String slots [] = wordId[i].split(":");
    					arrayBufferSlides[i] = Integer.parseInt(slots[0]);
						arrayCntBufferSlides[i] = Float
								.parseFloat(slots[1])
								/ truthBuffer.size();
    				}
    			}
    			else {
    				arrayBufferSlides = new int [0];
    				arrayCntBufferSlides = new double [0];
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
	        try { br.close(); }
	        catch (IOException e) { e.printStackTrace(); }
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
	
	public TranscriptionClass slidesInterpolation (double slidesTransRatio) {		
		int newObservationId [][] = new int[observationId.length][];
		double newObservationCount [][] = new double[observationCount.length][];

		for (int i = 0; i < observationId.length; i++) {
			HashMap<Integer, Double> count = new HashMap<Integer, Double>();
			for (int j = 0; j < observationId[i].length; j++) {
				if (count.containsKey(observationId[i][j])) {
					count.put(observationId[i][j],
							count.get(observationId[i][j])
									+ observationCount[i][j] * (1 - slidesTransRatio));
				} else {
					count.put(observationId[i][j], observationCount[i][j] * (1 - slidesTransRatio));
				}
			}
			
			for (int j = 0; j < observationSlidesId[i].length; j++) {
				if (count.containsKey(observationSlidesId[i][j])) {
					count.put(observationSlidesId[i][j],
							count.get(observationSlidesId[i][j])
									+ observationSlidesCount[i][j] * slidesTransRatio);
				} else {
					count.put(observationSlidesId[i][j], observationSlidesCount[i][j] * slidesTransRatio);
				}
			}

			int [] arrayBuffer = new int [count.size()];
			double [] arrayCntBuffer = new double [count.size()];
			int counter = 0;
			for (int key : count.keySet()) {
				arrayBuffer[counter] = key;
				arrayCntBuffer[counter] = count.get(key);
				counter++;
			}
			newObservationId[i] = arrayBuffer;
			newObservationCount[i] = arrayCntBuffer;
		}

		TranscriptionClass interpolatedTrans = new TranscriptionClass(truth,
				newObservationId, newObservationCount, courseId);
		return interpolatedTrans;
	}
	
	/** compute accuracy
	 * 
	 * @param hypo
	 *            hypothesis
 	 */
	/*public double [] acc (String [] hypo) {
        double numAcc = 0.0;
        double numSent = 0.0;
        
        for (int i = 0; i < hypo.length; i++) {
        	for (int j = 0; j < truth.get(i).size() ; j++) {
        		String wordId []= truth.get(i).get(j).split(" ");
        		for (int k = 0; k < wordId.length; k++) {
        			if (wordId[k].equals("0")) { continue; }
        			if (hypo[i].equals(wordId[k])) { numAcc += 1; }
        			numSent += 1;
        		}
        	}
        }

        double[] result = {numAcc, numSent};
        return result;	    
	}*/
	
	public double [] acc (String [] hypo) {
        double numAcc = 0.0;
        double numSent = 0.0;
        
        for (int i = 0; i < hypo.length; i++) {
        	for (int j = 0; j < truth.get(i).size() ; j++) {
        		String wordId []= truth.get(i).get(j).split(" ");
        		double truthClass = 0.0;
        		double hitClass = 0.0;
        		for (int k = 0; k < wordId.length; k++) {
        			if (wordId[k].equals("0")) { continue; }
        			if (hypo[i].equals(wordId[k])) { hitClass += 1; }
        			truthClass += 1;
        		}
        		if (truthClass != 0.0) {
        			numAcc += hitClass/truthClass;
        			numSent += 1;
        		}
        	}
        }

        double[] result = {numAcc, numSent};
        return result;	    
	}
	
	public void printResult (String [] hypo, double [][] alignedPosterior) {
		int sentCnt = 0;
		System.out.println("course : " + courseId);
        for (int i = 0; i < hypo.length; i++) {
        	for (int j = 0; j < truth.get(i).size() ; j++) {
        		sentCnt++;
        		System.out.print("sent id: " + sentCnt + ", slides: " + hypo[i] + ", posterior:");
        		for (int k = 0; k < alignedPosterior.length ; k++) {
        			System.out.print(" " + alignedPosterior[k][i]);
        		}
        		System.out.println("");
        	}
        }
	}
}
