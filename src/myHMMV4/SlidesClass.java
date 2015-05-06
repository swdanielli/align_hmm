package myHMMV4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/* NOTICE: the filtered page (e.g., too short, no keywords) should be listed with its 
 * index and an empty feature line (or duplicate the line it merges with) */

public class SlidesClass {
	public int observationId [][];
	public double observationCount [][];
	
	public SlidesClass(String slidesFileName) {
		List<int[]> arrayArrayBuffer = null;
		List<double[]> arrayArrayCntBuffer = null;

		/* Load slides, distribution */
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(slidesFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			String line = br.readLine();
			arrayArrayBuffer = new ArrayList<int[]>();
			arrayArrayCntBuffer = new ArrayList<double []>();
			String slidesBuffer = "";
			while (line != null) {
				String tags[] = line.split("_");
				slidesBuffer = tags[1].trim();
				int[] arrayBuffer = null;
				double [] arrayCntBuffer = null;
	        	
				if (slidesBuffer.length() > 0) {
					String wordId[] = slidesBuffer.split("\\s+");
					arrayBuffer = new int[wordId.length];
    				arrayCntBuffer = new double [wordId.length];

    				for (int i = 0; i < wordId.length; i++) {
    					String slots [] = wordId[i].split(":");
    					arrayBuffer[i] = Integer.parseInt(slots[0]);
    					arrayCntBuffer[i] = Float.parseFloat(slots[1]);
    				}
				} else {
					arrayBuffer = new int [0];
    				arrayCntBuffer = new double [0];
				}

				arrayArrayBuffer.add(arrayBuffer);
	        	arrayArrayCntBuffer.add(arrayCntBuffer);
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
		
		observationId = new int[arrayArrayBuffer.size()][];
        for (int i = 0; i < arrayArrayBuffer.size(); i++)
        	observationId[i] = arrayArrayBuffer.get(i);
	    
        observationCount = new double[arrayArrayCntBuffer.size()][];
        for (int i = 0; i < arrayArrayCntBuffer.size(); i++)
        	observationCount[i] = arrayArrayCntBuffer.get(i);

	}
}
